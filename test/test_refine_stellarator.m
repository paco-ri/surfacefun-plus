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

marked = 1:75;
dom = surfacemesh.adap_ref(dom, amr_tol, rmax, mode, marked);
surfacemesh_to_vtk(dom, fname, 'Refined stellarator geometry');

% TODO go back to quadtree generation