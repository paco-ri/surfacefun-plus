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