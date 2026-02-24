function plot_quadforest(obj, varargin)
%PLOT_QUADFOREST Summary of this function goes here
%   Detailed explanation goes here

if length(varargin) >= 1
    fignum = varargin{1};
else
    fignum = 1;
end
loc_shift = [0 -1; 1 0; 0 1; -1 0];

% compute shifts for each tree
shifts = zeros(obj.n_trees, 2);

patches_done = zeros(obj.n_trees, 1);
patches_done(1) = 1;

prev_patches = zeros(obj.n_trees, 1);
prev_patches(1) = 1;

curr_patches = zeros(obj.n_trees, 1);
% curr_patches(1) = 1;

% patches_to_inspect = zeros(obj.n_trees, 1);
% patches_to_inspect(1) = 1; % start computing shifts for first patch
num_patches_done = 1;

while num_patches_done < obj.n_trees
    for i = 1:obj.n_trees
        if prev_patches(i)
            for j = 1:4
                if ~patches_done(obj.C(i, j))
                    shifts(obj.C(i, j), :) = shifts(i, :) + loc_shift(j, :);
                    curr_patches(obj.C(i, j)) = 1;
                    patches_done(obj.C(i, j)) = 1;
                    num_patches_done = num_patches_done + 1;
                end
            end
        end
    end
    prev_patches = curr_patches;
    curr_patches = zeros(obj.n_trees, 1);
end

figure(fignum)
for i = 1:obj.n_trees
    obj.plot_quadtree(i, fignum, shifts(i, :), obj.L_max);
    hold on
end

end