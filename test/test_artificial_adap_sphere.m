clear
close all

% set up domain
n = 16; % polynomial order
dom = surfacemesh.sphere(n, 1);
npat0 = length(dom.x);

% divide patch 1 into four patches
[x_fin, y_fin, z_fin] = refine_patch(n, dom.x{1}, dom.y{1}, dom.z{1});
dom.x{1} = x_fin{1};
dom.y{1} = y_fin{1};
dom.z{1} = z_fin{1};
for i = 1:3
    dom.x{npat0 + i} = x_fin{i + 1};
    dom.y{npat0 + i} = y_fin{i + 1};
    dom.z{npat0 + i} = z_fin{i + 1};
end

elem2elem = dom.connectivity.elem2elem;

% get split info (L, R, D, U)
split = cell(length(dom.x), 1);
for k = 1:length(dom.x)
    split{k} = false(1, 4);
end
% neighbors of patch 1 each have one edge with a hanging node

% The elem2elem connectivity is ordered within each row in the order (L, R, D, U),
% i.e. elem2elem(j,:) = [L R D U]. Here, "left" means left in the local
% coordinate system of the element j-th element. If element j has been
% refined and thus has hanging nodes on each of its four sides, we can
% compute the split flags for each of its neighbors by finding which column of
% elem2elem(k,:) contains j, where k is L, R, D, U.
j = 1;
L = elem2elem(j, 1);
R = elem2elem(j, 2);
D = elem2elem(j, 3);
U = elem2elem(j, 4);
split{L} = (elem2elem(L,:) == j);
split{R} = (elem2elem(R,:) == j);
split{D} = (elem2elem(D,:) == j);
split{U} = (elem2elem(U,:) == j);

% generate new surfacemesh with updated x, y, z
dom = surfacemesh(dom.x, dom.y, dom.z, split);

% solve Laplace-Beltrami problem with spherical harmonic RHS
l = 3; m = 2;
sol = spherefun.sphharm(l, m);
sol = surfacefun(@(x,y,z) sol(x,y,z), dom);
f = -l*(l+1)*sol;
pdo = [];
pdo.lap = 1;
L = surfaceop(dom, pdo, f);
L.rankdef = true;
u = L.solve();
plot(u)
hold on
plot(dom)
shg

function [x_fin,y_fin,z_fin] = refine_patch(n,x_patch,y_patch,z_patch)
%REFINE_PATCH refine a patch

% tools for refining a patch
x  = chebpts(n, [-1 1]);
xL = chebpts(n, [-1 0]);
xR = chebpts(n, [ 0 1]);
BL = barymat(xL, x);
BR = barymat(xR, x);

% get nodes of fine patches
x_fin = cell(4,1);
y_fin = cell(4,1);
z_fin = cell(4,1);
x_fin{1} = BL * x_patch * BL.';
y_fin{1} = BL * y_patch * BL.';
z_fin{1} = BL * z_patch * BL.';
x_fin{2} = BL * x_patch * BR.';
y_fin{2} = BL * y_patch * BR.';
z_fin{2} = BL * z_patch * BR.';
x_fin{3} = BR * x_patch * BL.';
y_fin{3} = BR * y_patch * BL.';
z_fin{3} = BR * z_patch * BL.';
x_fin{4} = BR * x_patch * BR.';
y_fin{4} = BR * y_patch * BR.';
z_fin{4} = BR * z_patch * BR.';

end