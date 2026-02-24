function plot_quadtree_merge(obj, elem, level)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Q = obj.morton{elem};
merge_idx = obj.merge_idx{elem};
merged_patches = obj.merged_patches{elem};
L_max = obj.L_max;

quadforest.plot_quadtree(Q, level)
hold on

% get coordinates of patches to merge together
scl = 2^(L_max - level);
% look at merge_idx(L_max - level + 1) to study patches at level in Q
num_merge = idivide(size(merge_idx{2 * (L_max - level) + 1}, 1), int32(2));
merge_x = zeros(1, 4);
merge_y = zeros(1, 4);
for i = 1:num_merge
    % top-left of 4 merged patches is in odd-numbered rows of merge_idx
    idx = merge_idx{2 * (L_max - level) + 1}(2 * i - 1, 1);
    if idx < 0
        [x, y] = quadforest.deinterleave(merged_patches(-idx), L_max);
    else 
        [x, y] = quadforest.deinterleave(Q{level}(idx), L_max);
    end
    
    merge_x(1) = scl * x;
    merge_x(2) = scl * x + 2^(L_max - level + 1);
    merge_x(3) = scl * x + 2^(L_max - level + 1);
    merge_x(4) = scl * x;

    merge_y(1) = scl * y;
    merge_y(2) = scl * y;
    merge_y(3) = scl * y + 2^(L_max - level + 1);
    merge_y(4) = scl * y + 2^(L_max - level + 1);

    fill(merge_x, merge_y, 'b')

    for j = (2 * i - 1):(2 * i)
        for k = 1:2
            dot_in_center(merge_idx{2 * (L_max - level) + 1}(j, k), level, Q, merged_patches, scl)
        end
    end
end


end

function dot_in_center(idx, level, Q, merged_patches, scl)

L_max = size(Q, 1);
if idx < 0
    [x, y] = quadforest.deinterleave(merged_patches(-idx), L_max);
else
    [x, y] = quadforest.deinterleave(Q{level}(idx), L_max);
end
% plot(x + 2^(L_max - level), y + 2^(L_max - level), 'ro', 'MarkerSize', 8)
plot(scl * double(x) + 2^(L_max - level - 1), scl * double(y) + 2^(L_max - level - 1), 'p', 'MarkerSize', 20, 'MarkerFaceColor', 'r')

end