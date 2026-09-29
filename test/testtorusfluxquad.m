n = 4; 
nu = 6;
nv = nu*3;
dom = surfacemesh.torus(n,nu,nv);

f = surfacefunv(@(x,y,z) x.^2, ...
    @(x,y,z) x+y.^2, ...
    @(x,y,z) 1./(1 + z.^2), dom);

nr = 10;
nt = 40;
[qnodes, qweights] = torusfluxquad(nr,nt);
plot3(qnodes(1,:), qnodes(2,:), qnodes(3,:), '.')
hold on
plot(dom)
disp(sum(qweights))