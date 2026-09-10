function [dom, qf, p2q, split] = refine_leaves(dom0, p2q, marked, rmax)
%REFINE_LEAVES   Split marked leaves of a quadforest and rebuild the mesh.
%   [DOM, QF, P2Q, SPLIT] = REFINE_LEAVES(DOM0, P2Q, MARKED, RMAX) splits
%   each leaf listed in MARKED into its four children, restores the 2:1
%   balance of the resulting quadforest, and materializes the corresponding
%   SURFACEMESH.
%
%   Unlike ADAP_REF, no refinement criterion is applied here: the caller
%   decides what to mark. This makes the routine re-entrant, so a mesh can
%   be refined repeatedly by feeding the returned P2Q back in.
%
%   DOM0 is the base (level-0) mesh and is never modified. The mesh at any
%   stage is described entirely by the pair (DOM0, P2Q), where P2Q is an
%   npat x 3 matrix whose row i gives
%
%       [tree_root, level, morton]
%
%   for patch i of the current mesh. Row index equals patch index, so
%   DOM.x{i} corresponds with P2Q(i,:). Pass
%
%       p2q = [(1:length(dom0.x)).', zeros(length(dom0.x), 2)]
%
%   to start from the uniform mesh, or the P2Q returned by ADAP_REF to
%   start from a geometrically adapted one.
%
%   MARKED is a vector of patch indices into the CURRENT mesh (i.e. row
%   indices of the input P2Q). It may be empty, in which case the mesh is
%   simply re-materialized.
%
%   RMAX is the maximum quadtree depth, passed to QUADFOREST as L_max.
%
%   Neighbour queries here do NOT go through QUADFOREST/FOREST_COLLEAGUES or
%   QUADFOREST/GET_SPLIT. Those resolve cross-patch adjacency with
%   ROTATE_NODE, which rotates in a 2^(l+1) grid for a level-l code and so
%   returns the wrong neighbour once trees are more than one level deep;
%   meshes built that way carry level-1 cells adjacent to level-3 ones and
%   are rejected by SURFACEOP's HPS solve. Instead the adjacency of the base
%   mesh (which patch, which side, and whether the shared edge is traversed
%   in the opposite direction) is read off the corner geometry of DOM0 once,
%   and both the balance condition and the SPLIT flags are derived from it.
%
%   See also SURFACEMESH/ADAP_REF, QUADFOREST, SURFACEFUN/PROLONG_LEAVES.

arguments (Input)
    dom0
    p2q  (:,3) double
    marked
    rmax (1,1) double
end

arguments (Output)
    dom
    qf
    p2q
    split
end

npat0 = length(dom0.x);
n = size(dom0.x{1}, 1); % polynomial order

marked = marked(:).';
if ( any(marked < 1 | marked > size(p2q, 1)) )
    error('SURFACEMESH:refine_leaves:marked', ...
        'MARKED contains patch indices outside the current mesh.');
end

if ( any(p2q(marked, 2) >= rmax) )
    error('SURFACEMESH:refine_leaves:rmax', ...
        'MARKED contains leaves already at the maximum level RMAX = %d.', rmax);
end

% -- Split the marked leaves ------------------------------------------------
% Child 0 keeps the parent's row, so that unmarked patches keep their index;
% children 1-3 are appended. This matches the convention used by ADAP_REF.
p2q = split_leaves(p2q, marked);

% -- Enforce the 2:1 balance condition --------------------------------------
A = base_adjacency(dom0);
p2q = balance_leaves(p2q, A, rmax);

% -- Hanging-edge flags -----------------------------------------------------
split = compute_split(p2q, A);

% -- Quadforest -------------------------------------------------------------
% Returned for the benefit of callers that walk Morton codes within a tree
% (TaylorState.intacyc / intbcyc). Its own balancer is not trusted, so the
% Morton lists are restored to exactly the balanced leaf set computed above.
C = dom0.connectivity.elem2elem;
[Q, tree_roots] = morton_lists(p2q, npat0, rmax);
qf = quadforest(Q, rmax, C, tree_roots);
qf.morton = Q;

% -- Materialize the mesh ---------------------------------------------------
% Each leaf is interpolated straight out of its level-0 ancestor with a
% single pair of barycentric interpolation matrices. This is exact for the
% degree-(n-1) interpolant and identical to composing the BL/BR chain down
% the Morton path, but needs no traversal.
npat = size(p2q, 1);
x = cell(npat, 1);
y = cell(npat, 1);
z = cell(npat, 1);
xc = chebpts(n, [-1 1]);
for i = 1:npat
    t = p2q(i, 1);
    l = p2q(i, 2);
    if ( l == 0 )
        x{i} = dom0.x{t};
        y{i} = dom0.y{t};
        z{i} = dom0.z{t};
        continue
    end
    [mx, my] = quadforest.deinterleave(uint64(p2q(i, 3)), l);
    h = 2/2^l;
    % Morton x sits in the even bits and indexes columns (the u direction,
    % right factor); Morton y sits in the odd bits and indexes rows.
    Bu = barymat(chebpts(n, [-1 + h*double(mx), -1 + h*(double(mx)+1)]), xc);
    Bv = barymat(chebpts(n, [-1 + h*double(my), -1 + h*(double(my)+1)]), xc);
    x{i} = Bv * dom0.x{t} * Bu.';
    y{i} = Bv * dom0.y{t} * Bu.';
    z{i} = Bv * dom0.z{t} * Bu.';
end

dom = surfacemesh(x, y, z, split);

end

% ==========================================================================

function p2q = split_leaves(p2q, marked)
%SPLIT_LEAVES   Replace each marked leaf by its four children.

nadd = 3*numel(marked);
p2q = [p2q; nan(nadd, 3)];
k = size(p2q, 1) - nadd;
for p = marked
    t = p2q(p, 1);
    l = p2q(p, 2) + 1;
    m = 4*p2q(p, 3);
    p2q(p, :) = [t, l, m];
    for j = 1:3
        k = k + 1;
        p2q(k, :) = [t, l, m + j];
    end
end

end

function A = base_adjacency(dom0)
%BASE_ADJACENCY   Side-to-side adjacency of the conforming base mesh.
%   A is npat0 x 4 x 3. A(t,s,:) = [t', s', rev] says that side s of base
%   patch t is shared with side s' of base patch t', and that rev is 1 when
%   the two patches traverse the shared edge in opposite directions. A(t,s,1)
%   is 0 for a boundary edge.
%
%   Sides are ordered [Left Right Down Up] as elsewhere in SURFACEMESH:
%   Left/Right are the u = -1 / u = +1 edges (columns), Down/Up the
%   v = -1 / v = +1 edges (rows). Reading this off the geometry avoids every
%   Morton-rotation convention.

npat0 = length(dom0.x);

% Corners in the order (u-,v-), (u-,v+), (u+,v-), (u+,v+). Rows index v,
% columns index u, matching @surfacemesh/private/buildConnectivity.
cor = zeros(4*npat0, 3);
for t = 1:npat0
    X = dom0.x{t}; Y = dom0.y{t}; Z = dom0.z{t};
    cor(4*(t-1)+1, :) = [X(1,1)     Y(1,1)     Z(1,1)];
    cor(4*(t-1)+2, :) = [X(end,1)   Y(end,1)   Z(end,1)];
    cor(4*(t-1)+3, :) = [X(1,end)   Y(1,end)   Z(1,end)];
    cor(4*(t-1)+4, :) = [X(end,end) Y(end,end) Z(end,end)];
end
[~, ~, nodeid] = uniquetol(cor, 1e-10, 'ByRows', true, ...
    'DataScale', max(1, max(abs(cor(:)))));
nodeid = reshape(nodeid, 4, npat0).';

% Endpoints of each side, ordered along the edge: Left/Right by increasing
% v, Down/Up by increasing u.
sidecor = [1 2; 3 4; 1 3; 2 4];

ends = zeros(npat0, 4, 2);
for t = 1:npat0
    for s = 1:4
        ends(t, s, :) = nodeid(t, sidecor(s, :));
    end
end

A = zeros(npat0, 4, 3);
% Index sides by their unordered endpoint pair.
keys = containers.Map('KeyType', 'char', 'ValueType', 'any');
for t = 1:npat0
    for s = 1:4
        e = sort(squeeze(ends(t, s, :)).');
        k = sprintf('%d_%d', e(1), e(2));
        if ( isKey(keys, k) )
            keys(k) = [keys(k); t s];
        else
            keys(k) = [t s];
        end
    end
end
for t = 1:npat0
    for s = 1:4
        e = sort(squeeze(ends(t, s, :)).');
        lst = keys(sprintf('%d_%d', e(1), e(2)));
        other = lst(lst(:,1) ~= t | lst(:,2) ~= s, :);
        if ( isempty(other) )
            continue    % boundary edge
        end
        t2 = other(1,1); s2 = other(1,2);
        rev = double(ends(t, s, 1) ~= ends(t2, s2, 1));
        A(t, s, :) = [t2, s2, rev];
    end
end

end

function [t2, m2, ok] = neighbour_cell(t, l, m, s, A)
%NEIGHBOUR_CELL   The same-level cell across side s of cell (t,l,m).
%   ok is false when that side is on the boundary of the surface.

N = 2^l;
[mx, my] = quadforest.deinterleave(uint64(m), l);
mx = double(mx); my = double(my);

inside = [mx > 0, mx < N-1, my > 0, my < N-1];
if ( inside(s) )
    d = [-1 1 0 0; 0 0 -1 1];
    t2 = t;
    m2 = double(quadforest.interleave(mx + d(1,s), my + d(2,s), l));
    ok = true;
    return
end

t2 = A(t, s, 1);
if ( t2 == 0 )
    m2 = 0; ok = false; return
end
s2 = A(t, s, 2);
rev = A(t, s, 3);

% Index along the shared edge: v for Left/Right, u for Down/Up.
if ( s == 1 || s == 2 )
    a = my;
else
    a = mx;
end
if ( rev )
    a = N - 1 - a;
end

switch s2
    case 1, mx2 = 0;     my2 = a;
    case 2, mx2 = N-1;   my2 = a;
    case 3, mx2 = a;     my2 = 0;
    case 4, mx2 = a;     my2 = N-1;
end
m2 = double(quadforest.interleave(mx2, my2, l));
ok = true;

end

function key = leaf_map(p2q)
key = containers.Map('KeyType', 'char', 'ValueType', 'double');
for i = 1:size(p2q, 1)
    key(sprintf('%d_%d_%d', p2q(i,1), p2q(i,2), p2q(i,3))) = i;
end
end

function [lc, idx] = covering_leaf(key, t, l, m)
%COVERING_LEAF   The leaf containing cell (t,l,m), found by walking up.
%   Returns lc = -1 when no ancestor is a leaf, i.e. the region is refined
%   finer than level l.
lc = -1; idx = 0;
for c = l:-1:0
    mc = double(bitshift(uint64(m), -2*(l - c)));
    k = sprintf('%d_%d_%d', t, c, mc);
    if ( isKey(key, k) )
        lc = c; idx = key(k); return
    end
end
end

function p2q = balance_leaves(p2q, A, rmax)
%BALANCE_LEAVES   Split leaves until no two edge-neighbours differ by more
%   than one quadtree level.
%
%   Every violating pair is detected from its FINE side: for a leaf at level
%   l, the same-level cell across each edge is located, and the leaf actually
%   covering that cell is found by walking up the Morton path. A covering
%   leaf below level l-1 is too coarse and is split. Finding none means the
%   neighbourhood is finer than l, which that leaf handles itself.

for sweep = 1:(4*rmax + 8)
    key = leaf_map(p2q);
    tosplit = false(size(p2q, 1), 1);
    for i = 1:size(p2q, 1)
        l = p2q(i, 2);
        if ( l < 2 )
            continue    % nothing can be more than one level coarser
        end
        for s = 1:4
            [t2, m2, ok] = neighbour_cell(p2q(i,1), l, p2q(i,3), s, A);
            if ( ~ok ), continue, end
            [lc, idx] = covering_leaf(key, t2, l, m2);
            if ( lc >= 0 && lc < l - 1 )
                tosplit(idx) = true;
            end
        end
    end

    if ( ~any(tosplit) )
        return
    end

    idx = find(tosplit).';
    if ( any(p2q(idx, 2) >= rmax) )
        error('SURFACEMESH:refine_leaves:balance', ...
            ['Enforcing the 2:1 balance condition needs to refine past ' ...
             'RMAX = %d. Increase RMAX or mark fewer patches.'], rmax);
    end
    p2q = split_leaves(p2q, idx);
end

error('SURFACEMESH:refine_leaves:balance', ...
    'The 2:1 balance sweep did not converge.');

end

function split = compute_split(p2q, A)
%COMPUTE_SPLIT   Hanging-edge flags for SURFACEMESH.
%   split{i}(s) is true when side s of leaf i faces two finer neighbours.
%   On a 2:1 balanced mesh that is exactly the case where the same-level
%   neighbour cell has no covering leaf, i.e. that region is one level finer.

npat = size(p2q, 1);
key = leaf_map(p2q);
split = cell(npat, 1);
for i = 1:npat
    split{i} = false(1, 4);
    l = p2q(i, 2);
    for s = 1:4
        [t2, m2, ok] = neighbour_cell(p2q(i,1), l, p2q(i,3), s, A);
        if ( ~ok ), continue, end
        lc = covering_leaf(key, t2, l, m2);
        if ( lc < 0 )
            split{i}(s) = true;
        end
    end
end

end

function [Q, tree_roots] = morton_lists(p2q, npat0, rmax)
%MORTON_LISTS   Per-tree Morton code lists in QUADFOREST's layout.

Q = cell(npat0, 1);
for i = 1:npat0
    Q{i} = cell(rmax, 1);
    for j = 1:rmax
        Q{i}{j} = uint64([]);
    end
end
for i = 1:size(p2q, 1)
    level = p2q(i, 2);
    if ( level > 0 )
        Q{p2q(i, 1)}{level} = [Q{p2q(i, 1)}{level} uint64(p2q(i, 3))];
    end
end
tree_roots = unique(p2q(p2q(:, 2) > 0, 1)).';

end
