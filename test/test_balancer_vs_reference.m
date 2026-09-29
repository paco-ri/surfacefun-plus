% Differential test: REFINE_LEAVES, which balances with QUADFOREST, against
% an independent reference that reads the base-mesh adjacency off the corner
% geometry (the functions below the script, from refine_leaves.m before the
% switch to QUADFOREST).
%
% On every refinement step of 30 seeded trials, checks
%
%   leaves    the balanced leaf set equals the reference's
%   rows      row i of the new P2Q descends from row i of the old one via
%             child 0 (unmarked, unbalanced leaves keep their row)
%   split     SPLIT equals the reference hanging-edge flags
%   qf        QF.MORTON lists exactly the leaves below level 0, and
%             QF.TREE_ROOTS the trees that have any
%   colls     FOREST_COLLEAGUES equals the reference neighbour, every leaf
%             and side
%
% Base meshes: the uniform stellarator (every shared edge aligned), the same
% surface with each patch's parametrization rotated by a random multiple of
% 90 degrees (all 16 side pairings), and an open square (boundary edges).
fprintf('=== test_balancer_vs_reference start\n');

% ===== Creation of base meshes =====
rmax = 4; n = 5;
domA = surfacemesh.stellarator(n, 9, 3);
np0 = length(domA.x);
% Random rotations of the stellarator patches
rng(12345);
krot = randi(4, np0, 1) - 1;
xr = cell(np0, 1); yr = xr; zr = xr;
for i = 1:np0
    xr{i} = rot90(domA.x{i}, krot(i));
    yr{i} = rot90(domA.y{i}, krot(i));
    zr{i} = rot90(domA.z{i}, krot(i));
end
bases = {domA, surfacemesh(xr, yr, zr), surfacemesh.square(n, 2)};
bname = {'aligned stellarator', 'rotated stellarator', 'open square'};
modes = {'random', 'clustered', 'partial'};
allpass = true;

for b = 1:numel(bases)
    dom0 = bases{b};
    np0 = length(dom0.x);
    % QUADFOREST sees only elem2elem; the reference sees only A, built from
    % the corner geometry. Their neighbour patches must agree first (a
    % boundary edge is -1 in C and 0 in A).
    C = dom0.connectivity.elem2elem;
    A = base_adjacency(dom0);
    nC = nnz(max(C, 0) ~= A(:, :, 1));

    bad = struct('leaves', 0, 'rows', 0, 'split', 0, 'qf', 0, 'colls', 0);
    nmesh = 0; nleaf = 0; nbnd = 0;
    % Each trial starts from the base mesh and refines up to RMAX times,
    % feeding each output P2Q back in; every step is one mesh checked.
    % MODE picks how leaves are marked (names in MODES):
    %   1 random     uniform first step, then 10% of eligible leaves at random
    %   2 clustered  uniform first step, then the 10% nearest one random leaf
    %   3 partial    30% of base patches first, then as random
    for mode = 1:3
        for trial = 1:10
            rng(trial);
            p2q = [(1:np0).', zeros(np0, 2)];
            for lev = 1:rmax
                % After the first step, mark 10% of the leaves that can
                % still be split.
                elig = find(p2q(:, 2) < rmax);
                if isempty(elig), break, end
                k = max(1, round(0.10*numel(elig)));
                if lev == 1 && mode ~= 3
                    % Uniform first step.
                    marked = (1:size(p2q, 1)).';
                elseif lev == 1
                    % Partial: leaves most base patches at level 0, so level-0
                    % neighbours and the planting of new trees are exercised.
                    marked = randperm(np0, round(0.3*np0)).';
                elseif mode == 2
                    % Clustered: the K leaves nearest a random one, by patch
                    % centroid. Dense local refinement forces balancing to
                    % cascade over several levels and across patches.
                    cx = cellfun(@(a) mean(a(:)), dom.x);
                    cy = cellfun(@(a) mean(a(:)), dom.y);
                    cz = cellfun(@(a) mean(a(:)), dom.z);
                    s = elig(randi(numel(elig)));
                    dd = (cx(elig)-cx(s)).^2 + (cy(elig)-cy(s)).^2 + (cz(elig)-cz(s)).^2;
                    [~, ord] = sort(dd);
                    marked = elig(ord(1:k));
                else
                    marked = elig(randperm(numel(elig), k));
                end

                % The same step both ways: the reference splits and
                % balances on P2Q directly, REFINE_LEAVES via QUADFOREST.
                ref = balance_leaves(split_leaves(p2q, marked(:).'), A, rmax);
                [dom, qf, p2qN, split] = surfacemesh.refine_leaves(dom0, p2q, marked, rmax);
                nmesh = nmesh + 1;
                nleaf = nleaf + size(p2qN, 1);

                % Row order of the appended leaves may differ, so compare
                % as sets.
                bad.leaves = bad.leaves + ~isequal(sortrows(p2qN), sortrows(ref));

                % REFINE_LEAVES' row convention: row i stays in the region
                % of old patch i, in the same tree, with the old Morton code
                % followed by DL levels of child 0. So an unmarked leaf that
                % balancing leaves alone keeps its row.
                old = p2q;
                new = p2qN(1:size(old, 1), :);
                dl = new(:, 2) - old(:, 2);
                bad.rows = bad.rows + ~(all(new(:, 1) == old(:, 1)) && all(dl >= 0) ...
                    && all(new(:, 3) == old(:, 3).*4.^dl));

                % SPLIT is per row, so it is compared in P2QN's order.
                bad.split = bad.split + ~isequal(split, compute_split(p2qN, A));

                % TaylorState.intacyc / intbcyc walk QF.MORTON and
                % QF.TREE_ROOTS in place of P2Q, so both must describe the
                % same leaves. Level-0 leaves have no Morton code.
                Lq = zeros(0, 3);
                for t = 1:np0
                    for l = 1:rmax
                        m = double(qf.morton{t}{l}(:));
                        Lq = [Lq; repmat([t l], numel(m), 1), m]; %#ok<AGROW>
                    end
                end
                fine = p2qN(p2qN(:, 2) > 0, :);
                bad.qf = bad.qf + ~(isequal(sortrows(Lq), sortrows(fine)) && ...
                    isequal(qf.tree_roots, unique(fine(:, 1)).'));

                % Same-level neighbour across every side of every leaf. A
                % boundary side must come back as [0 0].
                for i = 1:size(p2qN, 1)
                    t = p2qN(i, 1); l = p2qN(i, 2); m = p2qN(i, 3);
                    fc = qf.forest_colleagues(uint64(m), t, l);
                    for s = 1:4
                        [t2, m2, ok] = neighbour_cell(t, l, m, s, A);
                        if ok
                            bad.colls = bad.colls + ~isequal(fc(s, :), [m2 t2]);
                        else
                            nbnd = nbnd + 1;
                            bad.colls = bad.colls + ~isequal(fc(s, :), [0 0]);
                        end
                    end
                end
                p2q = p2qN;
            end
        end
    end

    pass = nC == 0 && all(cellfun(@(v) v == 0, struct2cell(bad)));
    allpass = allpass && pass;
    fprintf(['%-20s %3d meshes, %6d leaves, %5d boundary sides; mismatches: ' ...
        'elem2elem %d, leaves %d, rows %d, split %d, qf %d, colls %d  %s\n'], ...
        bname{b}, nmesh, nleaf, nbnd, nC, bad.leaves, bad.rows, bad.split, ...
        bad.qf, bad.colls, passfail(pass));
end

fprintf('VERDICT: %s\n', passfail(allpass));
fprintf('DONE_SENTINEL\n');

% ==========================================================================

function s = passfail(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end

% ---- reference, from refine_leaves.m before the switch to QUADFOREST -----
% It shares no code with QUADFOREST: adjacency comes from matching patch
% corners in space rather than elem2elem, and the orientation of a shared
% edge from comparing its endpoints rather than a side-pairing rule. A leaf
% (t,l,m) is looked up by walking up its Morton path in a map of P2Q.

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
%   is 0 for a boundary edge. Sides are ordered [Left Right Down Up].

npat0 = length(dom0.x);

% Corners in the order (u-,v-), (u-,v+), (u+,v-), (u+,v+).
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
%   than one quadtree level. Each violation is found from its fine side.

for sweep = 1:(4*rmax + 8)
    key = leaf_map(p2q);
    tosplit = false(size(p2q, 1), 1);
    for i = 1:size(p2q, 1)
        l = p2q(i, 2);
        if ( l < 2 )
            continue
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
    p2q = split_leaves(p2q, find(tosplit).');
end

error('balance_leaves: the 2:1 balance sweep did not converge.');

end

function split = compute_split(p2q, A)
%COMPUTE_SPLIT   Hanging-edge flags: side s of leaf i faces two finer
%   neighbours exactly when the same-level cell across it has no covering
%   leaf.

npat = size(p2q, 1);
key = leaf_map(p2q);
split = cell(npat, 1);
for i = 1:npat
    split{i} = false(1, 4);
    l = p2q(i, 2);
    for s = 1:4
        [t2, m2, ok] = neighbour_cell(p2q(i,1), l, p2q(i,3), s, A);
        if ( ~ok ), continue, end
        if ( covering_leaf(key, t2, l, m2) < 0 )
            split{i}(s) = true;
        end
    end
end

end
