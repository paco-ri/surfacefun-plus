clear
close all

% set up domain
n = 16; % polynomial order
dom = surfacemesh.sphere(n, 1);
npat0 = length(dom.x);

% divide patch 1 into four patches
[x_fin, y_fin, z_fin] = refine_patch(n, dom.x{1}, dom.y{1}, dom.z{1});
[x_fin2, y_fin2, z_fin2] = refine_patch(n, x_fin{2}, y_fin{2}, z_fin{2});
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
for i = 1:3
    dom.x{27 + i} = x_fin2{i + 1};
    dom.y{27 + i} = y_fin2{i + 1};
    dom.z{27 + i} = z_fin2{i + 1};
end

p2q = zeros(30, 3);
p2q(1:24, 1) = 1:24;
p2q(1, 2) = 1;
p2q(25:30, 1) = 1;
p2q(25, 2) = 2;
p2q(26:27, 2) = 1;
p2q(28:30, 2) = 2;
p2q(25, 3) = 4;
p2q(26, 3) = 2;
p2q(27, 3) = 3;
p2q(28, 3) = 5;
p2q(29, 3) = 6;
p2q(30, 3) = 7;

rmax = 2;

% use p2q to populate Q
Q = cell(24, 1);
for i = 1:24
    Q{i} = cell(rmax, 1);
    for j = 1:rmax
        Q{i}{j} = uint64([]);
    end
end

for i = 1:30
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
Qact{1}{2} = [4, 5, 6, 7];

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