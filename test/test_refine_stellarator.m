clear

% set up domain
n = 8;
nu = 5;
nv = 3*nu;
dom = surfacemesh.stellarator(n,nu,nv);
C = dom.connectivity.elem2elem; % patch connectivity

% possible call signature
% dom = adap_ref(dom, amr_tol, rmax, mode, marked)

amr_tol = 1e-8;
rmax = 5;
mode = 1;

fname = 'test_refine_stellarator.vtk';
% surfacemesh_to_vtk(dom, fname, 'Unrefined stellarator geometry')

% 75 patches initially
marked = 1:10;
[dom, qf] = surfacemesh.adap_ref(dom, amr_tol, rmax, mode, marked);
% for i = 6:10
%     qf.plot_quadtree(i, i);
% end
% set(gcf, 'Color', 'w')
% ctr_ix = idivide(n, int32(2)) + 1;
% for i = 1:75
%     text(1.1 * dom.x{i}(ctr_ix, ctr_ix), 1.1 * dom.y{i}(ctr_ix, ctr_ix), 1.1 * dom.z{i}(ctr_ix, ctr_ix), int2str(i), 'Color', 'red');
% end
f = surfacefun(@(x, y, z) sin(x) + sin(y) + sin(z), dom);
pdo = [];
pdo.lap = 1;
L = surfaceop(dom, pdo, f);
L.rankdef = true;
u = L.solve();
surfacemesh_to_vtk(dom, fname, u, 'Title', 'Refined stellarator geometry');
% plot(u)
% hold on
% plot(dom)
