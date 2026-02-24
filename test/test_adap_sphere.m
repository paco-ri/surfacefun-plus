clear

% set up domain
n = 8;
dom = surfacemesh.sphere(n, 1);
C = dom.connectivity.elem2elem;
npat0 = 4 * 6;

amr_tol = 1e-8;
marked = 1;
rmax = 2;
mode = 1;

dom = surfacemesh.adap_ref(dom, amr_tol, rmax, mode, marked);
plot(dom)

% Import expected split data and compare to dom.split
raw = readmatrix('split_test_adap_sphere.txt');
true_split = num2cell(logical(raw), 2);
diffIdx = find(~cellfun(@(a,b) isequal(a,b), true_split, dom.split));
if ~isempty(diffIdx)
    fprintf('dom.split does not match true_split at %d entries:\n', numel(diffIdx));
    fprintf('  %5s  %-20s  %-20s\n', 'idx', 'true_split', 'dom.split');
    for k = 1:numel(diffIdx)
        i = diffIdx(k);
        fprintf('  %5d  [%d %d %d %d]          [%d %d %d %d]\n', ...
            i, true_split{i}, dom.split{i});
    end
    error('dom.split does not match true_split at entries: %s', num2str(diffIdx(:)'));
end
% dom.split = true_split; % for now just set dom.split to true_split to test balance_quadforest

% l = 3; m = 2;
% sol = spherefun.sphharm(l, m);
% sol = surfacefun(@(x,y,z) sol(x,y,z), dom);
% f = -l*(l+1)*sol;
f = surfacefun(@(x, y, z) sin(x) + sin(y) + sin(z), dom);
pdo = [];
pdo.lap = 1;
L = surfaceop(dom, pdo, f);
L.rankdef = true;
u = L.solve();
plot(u)
hold on
plot(dom)