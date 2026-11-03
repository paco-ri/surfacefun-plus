% Explain the error floor that TEST_PCONV shows on meshes with hanging edges
% (the split flags themselves pass TEST_BALANCER_VS_REFERENCE).
%   (a) does a CONFORMING mesh show the same floor?  (uniform levels 0, 1)
%   (b) on a mesh with hanging edges, is the error concentrated on the
%       patches next to a hanging edge (coarse side, fine side, or both),
%       or simply on the coarsest ones?
fprintf('=== test_floor start ===\n');

ufun = @(x,y,z) sin(0.7*x).*cos(0.5*y) + 0.3*z.^2 - 0.2*x.*z;
nv = 3; nu = 3*nv; rmax = 4;
orders = [6 8 10 12];

% random trial 7, level 3 of TEST_RANDOM_REFINEMENT
dom0 = taylor_stellarator(orders(1), nu, nv);
steps = random_refinement_steps(dom0, 'random', 7, rmax);
p2qr = steps(3).p2q;

% ---- (a) conforming control ----
fprintf('\n(a) p-convergence, all CONFORMING (no hanging edges at all)\n');
fprintf('%-24s', 'mesh \ n'); fprintf('%12d', orders); fprintf('\n');
for lev = 0:1
    errs = nan(size(orders)); npat = 0;
    for j = 1:numel(orders)
        dom0 = taylor_stellarator(orders(j), nu, nv);
        np0 = length(dom0.x);
        if lev == 0, marked = []; else, marked = 1:np0; end
        dm = surfacemesh.refine_leaves(dom0, [(1:np0).' zeros(np0,2)], ...
            marked, rmax);
        npat = length(dm.x);
        errs(j) = manufactured_error(dm, ufun);
    end
    fprintf('%-24s', sprintf('uniform level %d (%d)', lev, npat));
    fprintf('%12.3e', errs); fprintf('\n');
end

% ---- (b) where does the error live on a mesh with hanging edges? ----
nord = orders(end);
dom0 = taylor_stellarator(nord, nu, nv);
np0 = length(dom0.x);
p0 = [(1:np0).', zeros(np0,2)];
meshes = {'every other patch', 'random_t7_l3'};
% children of a cell that touch its side [Left Right Down Up]; children of
% cell m are 4*m + k with k laid out as  2 3 / 0 1  (see QUADFOREST.GET_SPLIT)
kids = [0 2; 1 3; 0 1; 2 3];
classes = {'no hanging edge', 'coarse side', 'fine side', 'both sides'};
for m = 1:numel(meshes)
    if m == 1
        [dm, qf, p2q] = surfacemesh.refine_leaves(dom0, p0, 1:2:np0, rmax);
    else
        [dm, qf, p2q] = surfacemesh.refine_leaves(dom0, p2qr, [], rmax);
    end
    [~, perr] = manufactured_error(dm, ufun);
    npat = length(dm.x);

    % SPLIT marks the coarse side of a hanging edge. Its fine side is the
    % two children, on the touching side, of the same-level neighbour.
    sp = dm.split;
    coarse = cellfun(@(s) any(logical(s)), sp(:));
    fine = false(npat, 1);
    nside = 0; nfound = 0;
    for i = find(coarse).'
        [col, back] = qf.forest_colleagues(uint64(p2q(i,3)), p2q(i,1), p2q(i,2));
        for s = find(sp{i})
            nb = [col(s,2)*[1; 1], (p2q(i,2) + 1)*[1; 1], ...
                  4*double(col(s,1)) + kids(back(s),:).'];
            hit = ismember(p2q, nb, 'rows');
            fine(hit) = true;
            nside = nside + 1; nfound = nfound + nnz(hit);
        end
    end
    if nfound ~= 2*nside
        fprintf('  WARNING: %d fine neighbours found for %d split sides (expect %d)\n', ...
            nfound, nside, 2*nside);
    end
    cls = 1 + coarse + 2*fine;

    fprintf('\n(b) error by patch class, %s, n = %d\n', meshes{m}, nord);
    fprintf('  %-26s %5s %12s %12s\n', 'class', 'npat', 'max err', 'median err');
    for lv = unique(p2q(:,2)).'
        for c = 1:numel(classes)
            sel = p2q(:,2) == lv & cls == c;
            if ~any(sel), continue, end
            fprintf('  level %d, %-17s %5d %12.3e %12.3e\n', lv, classes{c}, ...
                sum(sel), max(perr(sel)), median(perr(sel)));
        end
    end
end
fprintf('DONE_SENTINEL\n');
