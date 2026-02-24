clear
close all

% set up domain
n = 16; % polynomial order
dom = surfacemesh.sphere(n, 1);
npat0 = length(dom.x);

% divide patch 1 into four patches
[x_fin, y_fin, z_fin] = refine_patch(n, dom.x{1}, dom.y{1}, dom.z{1});
[x_fin2, y_fin2, z_fin2] = refine_patch(n, x_fin{2}, y_fin{2}, z_fin{2});
[x_fin21, y_fin21, z_fin21] = refine_patch(n, x_fin2{1}, y_fin2{1}, z_fin2{1});
dom.x{1} = x_fin{1};
dom.y{1} = y_fin{1};
dom.z{1} = z_fin{1};
dom.x{25} = x_fin2{1};
dom.y{25} = y_fin2{1};
dom.z{25} = z_fin2{1};
dom.x{26} = x_fin{3};
dom.y{26} = y_fin{3};
dom.z{26} = z_fin{3};
dom.x{27} = x_fin{4};
dom.y{27} = y_fin{4};
dom.z{27} = z_fin{4};
dom.x{28} = x_fin21{1};
dom.y{28} = y_fin21{1};
dom.z{28} = z_fin21{1};
dom.x{29} = x_fin2{3};
dom.y{29} = y_fin2{3};
dom.z{29} = z_fin2{3};
dom.x{30} = x_fin2{4};
dom.y{30} = y_fin2{4};
dom.z{30} = z_fin2{4};
for i = 1:3
    dom.x{30 + i} = x_fin21{i + 1};
    dom.y{30 + i} = y_fin21{i + 1};
    dom.z{30 + i} = z_fin21{i + 1};
end

npat = 33;
p2q = zeros(npat, 3);
p2q(1:24, 1) = 1:24;
p2q(1, 2) = 1;
p2q(25:33, 1) = 1;
p2q(25, 2) = 3;
p2q(26:27, 2) = 1;
p2q(28:30, 2) = 2;
p2q(31:33, 2) = 3;
p2q(25:33, 3) = [16, 2, 3, 5, 6, 7, 17, 18, 19];

rmax = 3;

% use p2q to populate Q
Q = cell(npat0, 1);
for i = 1:npat0
    Q{i} = cell(rmax, 1);
    for j = 1:rmax
        Q{i}{j} = uint64([]);
    end
end

for i = 1:npat
    level = p2q(i, 2);
    if level > 0
        Q{p2q(i, 1)}{level} = [Q{p2q(i, 1)}{level} p2q(i, 3)];
    end
end

% create quadforest object manually
Qact = cell(24, 1);
for i = 1:24
    Qact{i} = cell(rmax, 1);
    for j = 1:rmax
        Qact{i}{j} = uint64([]);
    end
end
Qact{1}{1} = [0, 2, 3];
Qact{1}{2} = [5, 6, 7];
Qact{1}{3} = [16, 17, 18, 19];

% for i in 1:npat0, if p2q(i, 2) is not zero, then there is a tree at Q{i}
tree_roots = p2q(:, 2) > 0;
qf = quadforest(Q, rmax, dom.connectivity.elem2elem, tree_roots);
% TODO get new patches required for balancing
% idea: at each level, collect list of patches to be refined and call
% refine routines below
% example: refine 1 --> 1 76 77 78; refine 76 --> 76 79 80 81
addl_p2q = nan([npat0*4^rmax 3]);
i = 1;
% columns of p2q are tree root index, quadtree level, Morton code
for t = 1:npat0
    for l = 1:rmax-1
        if ~isempty(qf.addl_patches{t}{l})
            p2q(t, 2) = p2q(t, 2) + 1; % patch is now one level deeper
            for p = qf.addl_patches{t}{l}
                addl_p2q(i, :) = [t l p];
                i = i + 1;
            end
        end
    end
end
addl_p2q = rmmissing(addl_p2q);
p2q = [p2q; addl_p2q];

% test equality of Q and Qact
disp(isequal(Q, Qact))

disp('addl_p2q')
disp(addl_p2q)

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