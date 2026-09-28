% Differential test: quadforest's neighbour, split and balance machinery
% against the verified geometric reference in refine_leaves.
%
%   T0  elem2elem vs base_adjacency; shared-edge pairings and reversals
%   T1  forest_colleagues vs neighbour_cell, every leaf x side
%   T2  get_split vs compute_split
%   T3  quadforest's own balancer vs balance_leaves, from the pre-balance
%       leaf set of every refinement step
%
% Run on two base meshes: the stellarator as built (every shared edge
% aligned, L-R and D-U), and the same surface with each patch's
% parametrization rotated by a random multiple of 90 degrees, which produces
% every side pairing and so exercises the rotation branches.
%
% Meshes: the seeded random/clustered sequences of test_random_refinement.m.
fprintf('=== diff_balancer start\n');

rmax = 4; n = 5; nv = 3; nu = 3*nv;
d = prepare_stellarator(n,nu,nv,16,40); domA = d{1};
np0 = length(domA.x);
rng(12345);
krot = randi(4, np0, 1) - 1;
xr = cell(np0,1); yr = xr; zr = xr;
for i = 1:np0
    xr{i} = rot90(domA.x{i}, krot(i));
    yr{i} = rot90(domA.y{i}, krot(i));
    zr{i} = rot90(domA.z{i}, krot(i));
end
domB = surfacemesh(xr, yr, zr);
bases = {domA, domB};
bname = {'ALIGNED base mesh', 'ROTATED base mesh'};
sname = {'L','R','D','U'};
modes = {'random', 'clustered', 'partial'};
allpass = true;

for b = 1:2
dom0 = bases{b};
fprintf('\n########## %s ##########\n', bname{b});
C = dom0.connectivity.elem2elem;
A = base_adjacency(dom0);

% ---- T0 -------------------------------------------------------------------
cmis = 0; nrev = 0; pairs = zeros(4,4);
for t = 1:np0
    for s = 1:4
        if C(t,s) ~= A(t,s,1), cmis = cmis + 1; end
        nrev = nrev + A(t,s,3);
        pairs(s, A(t,s,2)) = pairs(s, A(t,s,2)) + 1;
    end
end
fprintf('T0 elem2elem vs base_adjacency: %d of %d entries differ; reversed edges %d\n', ...
    cmis, 4*np0, nrev);
fprintf('T0 side pairings (row = my side L R D U, col = neighbour side L R D U):\n');
disp(pairs)

% ---- meshes ---------------------------------------------------------------
cases = struct('name', {}, 'pre', {}, 'post', {});
for mode = 1:3
    mname = modes{mode};
    for trial = 1:10
        rng(trial);
        p2q = [(1:np0).', zeros(np0,2)];
        for lev = 1:rmax
            npc = size(p2q,1);
            elig = find(p2q(:,2) < rmax);
            if isempty(elig), break, end
            if lev == 1 && mode ~= 3
                marked = (1:npc).';
            elseif lev == 1
                % partial: leaves most base patches at level 0, so the
                % level-0 paths and the planting of new trees are exercised
                marked = randperm(npc, round(0.3*npc)).';
            elseif mode == 1 || mode == 3
                k = max(1, round(0.10*numel(elig)));
                marked = elig(randperm(numel(elig), k));
            else
                k = max(1, round(0.10*numel(elig)));
                cx = cellfun(@(a) mean(a(:)), domR.x);
                cy = cellfun(@(a) mean(a(:)), domR.y);
                cz = cellfun(@(a) mean(a(:)), domR.z);
                s = elig(randi(numel(elig)));
                dd = (cx(elig)-cx(s)).^2 + (cy(elig)-cy(s)).^2 + (cz(elig)-cz(s)).^2;
                [~, ord] = sort(dd);
                marked = elig(ord(1:k));
            end
            pre = split_leaves(p2q, marked(:).');
            evalc('[domR, ~, p2qR] = surfacemesh.refine_leaves(dom0, p2q, marked, rmax);');
            cases(end+1) = struct('name', sprintf('%s t%d L%d', mname, trial, lev), ...
                'pre', pre, 'post', p2qR); %#ok<SAGROW>
            p2q = p2qR;
        end
    end
end
fprintf('%d meshes (every refinement step of 30 trials)\n', numel(cases));

% ---- T1 / T2 / T3 ---------------------------------------------------------
% forest_colleagues' row r is side r: L R D U.
dirside = [1 2 3 4];
t1within = zeros(2, 4);          % [mismatches; queries] per row
t1cross = zeros(4, 4, 2);        % (row, neighbour side, [mismatches queries])
t1ex = {};
t2 = struct('rows', 0, 'bad', 0, 'byside', zeros(1,4), 'miss', 0, 'extra', 0);
t2ex = {};
t3 = struct('err', 0, 'notpart', 0, 'viol', 0, 'differs', 0, 'extra', 0, 'missing', 0);
t3ex = {};
for c = 1:numel(cases)
    post = cases(c).post; pre = cases(c).pre;
    [Q, roots] = morton_lists(post, np0, rmax);
    evalc('qf = quadforest(Q, rmax, C, roots);');
    qf.morton = Q;

    for i = 1:size(post,1)
        t = post(i,1); l = post(i,2); m = post(i,3);
        fc = qf.forest_colleagues(uint64(m), t, l);
        [mx, my] = quadforest.deinterleave(uint64(m), l);
        N = 2^l;
        onb = [mx == 0, mx == N-1, my == 0, my == N-1];
        for r = 1:4
            s = dirside(r);
            [t2_, m2_] = neighbour_cell(t, l, m, s, A);
            bad = ~(fc(r,1) == m2_ && fc(r,2) == t2_);
            if onb(s)
                s2 = A(t,s,2);
                t1cross(r, s2, :) = t1cross(r, s2, :) + reshape([bad 1], 1, 1, 2);
            else
                t1within(:, r) = t1within(:, r) + [bad; 1];
            end
            if bad && numel(t1ex) < 6
                t1ex{end+1} = sprintf(['  %s  leaf (t=%d,l=%d,m=%d,x=%d,y=%d) ' ...
                    'row %d: fc -> (t=%d,m=%d)  ref side %s -> (t=%d,m=%d)'], ...
                    cases(c).name, t, l, m, mx, my, r, fc(r,2), fc(r,1), ...
                    sname{s}, t2_, m2_); %#ok<SAGROW>
            end
        end
    end

    % T2
    evalc('sq = qf.get_split(post, C);');
    sr = compute_split(post, A);
    for i = 1:size(post,1)
        t2.rows = t2.rows + 1;
        dd = sq{i} ~= sr{i};
        if any(dd)
            t2.bad = t2.bad + 1;
            t2.byside = t2.byside + dd;
            t2.miss = t2.miss + sum(sr{i} & ~sq{i});
            t2.extra = t2.extra + sum(sq{i} & ~sr{i});
            if numel(t2ex) < 6
                t2ex{end+1} = sprintf('  %s  leaf %d (t=%d,l=%d,m=%d): get_split %s  ref %s', ...
                    cases(c).name, i, post(i,1), post(i,2), post(i,3), ...
                    mat2str(sq{i}), mat2str(sr{i})); %#ok<SAGROW>
            end
        end
    end

    % T3: the constructor's balancer on the PRE-balance leaves.
    [Qp, rootsp] = morton_lists(pre, np0, rmax);
    try
        evalc('qb = quadforest(Qp, rmax, C, rootsp);');
    catch ME
        t3.err = t3.err + 1;
        if numel(t3ex) < 6
            t3ex{end+1} = sprintf('  %s  balancer ERROR: %s', cases(c).name, ME.message); %#ok<SAGROW>
        end
        continue
    end
    L = leaves_from_morton(qb.morton, np0);
    area = sum(4.^(-L(:,2)));
    part = abs(area - np0) < 1e-12;
    v = count_violations(L, A);
    same = isequal(sortrows(L), sortrows(post));
    t3.notpart = t3.notpart + ~part;
    t3.viol = t3.viol + (v > 0);
    t3.differs = t3.differs + ~same;
    if ~same
        t3.extra = t3.extra + size(setdiff(L, post, 'rows'), 1);
        t3.missing = t3.missing + size(setdiff(post, L, 'rows'), 1);
        if numel(t3ex) < 6
            t3ex{end+1} = sprintf(['  %s  pre %d leaves -> quadforest %d (area %.4f of %d, ' ...
                '%d 2:1 violations)  vs balance_leaves %d'], cases(c).name, size(pre,1), ...
                size(L,1), area, np0, v, size(post,1)); %#ok<SAGROW>
        end
    end
end

fprintf('\nT1 forest_colleagues vs neighbour_cell, by direction (mismatches/queries)\n');
fprintf('   within tree:   L %d/%d  R %d/%d  D %d/%d  U %d/%d\n', t1within);
fprintf('   across trees, by the neighbour side the edge lands on:\n');
fprintf('                  nbr L          nbr R          nbr D          nbr U\n');
rl = {'side L', 'side R', 'side D', 'side U'};
for r = 1:4
    fprintf('   %-9s', rl{r});
    for s2 = 1:4
        fprintf('  %5d/%-7d', t1cross(r,s2,1), t1cross(r,s2,2));
    end
    fprintf('\n');
end
t1bad = sum(t1within(1,:)) + sum(sum(t1cross(:,:,1)));
fprintf('   T1 total mismatches: %d\n', t1bad);
fprintf('%s\n', t1ex{:});

fprintf('\nT2 get_split vs compute_split: %d of %d leaves differ; by side L R D U = %s; missing flags %d, spurious flags %d\n', ...
    t2.bad, t2.rows, mat2str(t2.byside), t2.miss, t2.extra);
fprintf('%s\n', t2ex{:});

fprintf('\nT3 quadforest balancer vs balance_leaves over %d meshes:\n', numel(cases));
fprintf('   errors %d; not a partition %d; with 2:1 violations %d; differs from balance_leaves %d (extra leaves %d, missing %d)\n', ...
    t3.err, t3.notpart, t3.viol, t3.differs, t3.extra, t3.missing);
fprintf('%s\n', t3ex{:});

pass = t1bad == 0 && t2.bad == 0 && t3.err == 0 && t3.notpart == 0 && t3.viol == 0;
fprintf('\n%s: T1 %s, T2 %s, T3 %s\n', bname{b}, ...
    ternary(t1bad == 0), ternary(t2.bad == 0), ...
    ternary(t3.err == 0 && t3.notpart == 0 && t3.viol == 0));
allpass = allpass && pass;
end

if allpass, fprintf('\nVERDICT: PASS\n'); else, fprintf('\nVERDICT: FAIL\n'); end
fprintf('DONE_SENTINEL\n');

% ==========================================================================

function s = ternary(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end


function L = leaves_from_morton(M, np0)
%LEAVES_FROM_MORTON Leaf set of a quadforest's Morton lists, whether they
%   hold leaves only or every node: a code is a leaf when no descendant of it
%   is present.
rows = zeros(0,3);
for t = 1:np0
    for l = 1:numel(M{t})
        c = double(M{t}{l}(:));
        rows = [rows; repmat([t l], numel(c), 1) c]; %#ok<AGROW>
    end
end
key = leaf_map(rows);
isleaf = true(size(rows,1), 1);
hastree = false(np0, 1);
for i = 1:size(rows,1)
    t = rows(i,1); l = rows(i,2); m = rows(i,3);
    hastree(t) = true;
    for c = l-1:-1:1
        k = sprintf('%d_%d_%d', t, c, bitshift(m, -2*(l-c)));
        if isKey(key, k), isleaf(key(k)) = false; end
    end
end
L = [rows(isleaf,:); [find(~hastree), zeros(nnz(~hastree), 2)]];
L = unique(L, 'rows');
end

function v = count_violations(p2q, A)
key = leaf_map(p2q);
v = 0;
for i = 1:size(p2q,1)
    l = p2q(i,2);
    if l < 2, continue, end
    for s = 1:4
        [t2, m2, ok] = neighbour_cell(p2q(i,1), l, p2q(i,3), s, A);
        if ~ok, continue, end
        lc = covering_leaf(key, t2, l, m2);
        if lc >= 0 && lc < l - 1, v = v + 1; end
    end
end
end

% ---- reference implementation, copied verbatim from refine_leaves.m ------
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
