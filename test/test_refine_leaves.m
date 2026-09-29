% Test refine_leaves on a uniform base mesh.
%   TEST1  marking nothing reproduces the base mesh node-for-node
%   TEST2  marking everything reproduces refine(dom0, 1), up to patch order
fprintf('=== test_refine_leaves start ===\n');

% geometry parameters
n = 6; nu = 6; nv = 3;
dom0 = surfacemesh.sharptorus(n, nu, nv);
npat0 = length(dom0.x); % number of base patches
fprintf('base mesh: %d patches, n = %d\n', npat0, n);

p2q0 = [(1:npat0).', zeros(npat0, 2)]; % uniform base mesh
rmax = 4; % max quadtree depth

%% 1. identity: marked = [] reproduces dom0
d1 = surfacemesh.refine_leaves(dom0, p2q0, [], rmax);
% max node difference over all patches and coordinates
err = 0;
for k = 1:npat0
    err = max([err, max(abs(d1.x{k}-dom0.x{k}),[],'all'), ...
                    max(abs(d1.y{k}-dom0.y{k}),[],'all'), ...
                    max(abs(d1.z{k}-dom0.z{k}),[],'all')]);
end
fprintf('TEST1 identity: npat %d (want %d), max node err %.3e\n', ...
    length(d1.x), npat0, err);

%% 2. marking everything reproduces uniform refine(dom0,1)
d2 = surfacemesh.refine_leaves(dom0, p2q0, 1:npat0, rmax);
dref = refine(dom0, 1);
fprintf('TEST2 uniform: npat %d (want %d)\n', length(d2.x), length(dref.x));
% the two meshes order patches differently, so compare sorted centroids
c2 = sortrows(cellfun(@(a) mean(a(:)), [d2.x, d2.y, d2.z]));
cr = sortrows(cellfun(@(a) mean(a(:)), [dref.x, dref.y, dref.z]));
fprintf('TEST2 centroid match: %.3e\n', max(abs(c2-cr),[],'all'));

fprintf('DONE_SENTINEL\n');
