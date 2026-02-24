dom = surfacemesh.sphere(8, 1);
n_trees = length(dom);
% get element-to-element connectivity info
C = dom.connectivity.elem2elem;

L_max = 3;
% initialize quadforest
morton = cell(n_trees, 1);
for i = 1:n_trees(1)
    morton{i} = cell(L_max, 1);
    if i == 1
        morton{i}{1} = uint64([0b00 0b01]);
        morton{i}{2} = uint64([0b1000, 0b1001, 0b1010, 0b1100, 0b1101, 0b1110]);
        morton{i}{3} = uint64([0b101100, 0b101101, 0b101110, 0b101111, 0b111100, 0b111101, 0b111110, 0b111111]);
    else
        for j = 1:3
            morton{i}{j} = uint64([]);
        end
    end
end
Q = quadforest(morton, L_max, C, 1);
Q.plot_quadtree(Q.morton{5})
% TODO automatically generate shifts for quadforest plotting

% TODO call surfaceop constructor on dom + add'l small patches
% - for now, can start by calling surfacemesh(x, y, z) on manually-
%   calculated x, y, z

% TODO call buildWithQF(dom, Q)