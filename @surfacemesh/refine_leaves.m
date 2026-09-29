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
%   The 2:1 balance and the SPLIT flags come from QUADFOREST/BALANCE_QUADFOREST
%   and QUADFOREST/GET_SPLIT. QF.MORTON lists exactly the leaves of P2Q
%   below level 0.
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

% -- Split the marked leaves and restore the 2:1 balance -------------------
% Child 0 keeps the parent's row, so that unmarked patches keep their index;
% children 1-3 are appended.
p2q = split_leaves(p2q, marked);

Q = cell(npat0, 1);
for t = 1:npat0
    Q{t} = repmat({uint64([])}, rmax, 1);
end
for i = find(p2q(:, 2) > 0).'
    t = p2q(i, 1);
    l = p2q(i, 2);
    Q{t}{l} = [Q{t}{l} uint64(p2q(i, 3))];
end
qf = quadforest(Q, rmax, dom0.connectivity.elem2elem);
qf = qf.balance_quadforest();

% Balancing only splits, so split every leaf of P2Q that is not a leaf of QF
% until the two agree. This keeps the row convention above.
roots = qf.tree_roots;
leaves = [setdiff(1:npat0, roots).', zeros(npat0 - numel(roots), 2)];
for t = roots
    for l = 1:rmax
        m = double(qf.morton{t}{l}(:));
        leaves = [leaves; repmat([t l], numel(m), 1), m]; %#ok<AGROW>
    end
end
while true
    idx = find(~ismember(p2q, leaves, 'rows')).';
    if ( isempty(idx) )
        break
    end
    if ( any(p2q(idx, 2) >= rmax) )
        error('SURFACEMESH:refine_leaves:balance', ...
            'The balanced quadforest does not refine P2Q.');
    end
    p2q = split_leaves(p2q, idx);
end

split = qf.get_split(p2q);

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
