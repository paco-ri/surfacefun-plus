n = 8;
dom = surfacemesh.sphere(n, 2);
surfacemesh_to_vtk(dom, 'test_vtk.vtk', 'Test VTK file')
plot(dom)