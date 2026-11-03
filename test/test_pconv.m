% p-convergence of the hanging-edge handling on FIXED meshes.
%
% Hold each mesh fixed and raise the polynomial order. A wrong split flag
% makes the merge inconsistent, so the error stalls however large n gets.
% With correct flags every row decays, but a coarse patch is coupled to its
% fine neighbours at points between its own nodes, so how fast depends on how
% well the coarse patch resolves the geometry along the edge. Meshes whose
% coarse side is level 0 ('every other patch') lag the conforming row; on
% this surface that comes from the kink between level-0 patches, which
% refinement does not remove (leaves are restricted from the level-0 mesh).
fprintf('=== test_pconv start ===\n');

ufun = @(x,y,z) sin(0.7*x).*cos(0.5*y) + 0.3*z.^2 - 0.2*x.*z;
nv = 3; nu = 3*nv; rmax = 4; orders = [6 8 10 12];

% random trial 7, level 3 of TEST_RANDOM_REFINEMENT; random leaf sets do not
% depend on the polynomial order
dom0 = taylor_stellarator(orders(1), nu, nv);
steps = random_refinement_steps(dom0, 'random', 7, rmax);
p2qr = steps(3).p2q;

fprintf('\n%-28s', 'mesh \ n');
fprintf('%12d', orders); fprintf('   %8s\n', 'hang');

cases = {'uniform level 1', 'every other patch', 'random_t7_l3'};
for c = 1:numel(cases)
    errs = nan(size(orders)); nh = 0; npat = 0;
    for j = 1:numel(orders)
        nord = orders(j);
        dom0 = taylor_stellarator(nord, nu, nv);
        np0 = length(dom0.x);
        p0 = [(1:np0).', zeros(np0,2)];
        switch c
            case 1, dm = surfacemesh.refine_leaves(dom0, p0, 1:np0, rmax);
            case 2, dm = surfacemesh.refine_leaves(dom0, p0, 1:2:np0, rmax);
            case 3, dm = surfacemesh.refine_leaves(dom0, p2qr, [], rmax);
        end
        npat = length(dm.x);
        nh = 0; sp = dm.split;
        for k = 1:numel(sp), nh = nh + sum(logical(sp{k})); end
        try
            errs(j) = manufactured_error(dm, ufun);
        catch ME
            fprintf('\n  [%s n=%d FAILED: %s]\n', cases{c}, nord, ME.message);
        end
    end
    fprintf('%-28s', sprintf('%s (%d pat)', cases{c}, npat));
    fprintf('%12.3e', errs);
    fprintf('   %8d\n', nh);
end

fprintf(['\nEvery row should decay. A row that flattens means a wrong split\n' ...
         'flag; rows with level-0 coarse patches decay more slowly.\n']);
fprintf('DONE_SENTINEL\n');
