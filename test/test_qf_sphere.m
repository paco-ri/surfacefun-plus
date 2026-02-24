clear

dom = surfacemesh.sphere(8, 1);
C = dom.connectivity.elem2elem;
npat0 = 4 * 6;
Q = cell(npat0, 1);
Q{1} = cell(2, 1);
Q{1}{2} = uint64(0:15);

%quadforest.plot_quadtree(Q{1});
qf = quadforest(Q, 2, C, 1);
%qf.plot_quadtree(9)
qf.plot_quadforest();