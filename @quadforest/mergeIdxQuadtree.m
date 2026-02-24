function [merge_idx, merged_patches] = mergeIdxQuadtree(Q)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

L_max = size(Q, 1);
while L_max > 1
    if size(Q{L_max}) == 0
        L_max = L_max - 1;
    else
        break
    end
end
% L_max is the first level with leaves

merged_patches = []; % at each level, saves the patches that are the result of merging
if L_max == 1
    merge_idx = {};
    return
end
merge_idx = cell(2*L_max, 1);

% go to leaf level, which should have 4n leaves
% merge done in the obvious way
nrows_merge_idx = size(Q{L_max}, 2) / 2;
merge_idx{1} = zeros(nrows_merge_idx, 2);
merge_idx{1}(:, 2) = 2 * (1:size(Q{L_max}, 2)/2);
merge_idx{1}(:, 1) = merge_idx{1}(:, 2) - 1;
merge_idx{2} = merge_idx{1}(1:size(merge_idx{1}, 1)/2, :);
for i = 1:(size(Q{L_max}, 2) / 4)
    merged_patches = [merged_patches bitshift(Q{L_max}(4*(i-1)+1), -2)]; 
end
nnew_merged_patches = size(merged_patches, 2);
% disp("merge_idx{1}")
% disp(merge_idx{1})
% disp("merge_idx{2}")
% disp(merge_idx{2})

km = -1; k2 = 1; % counters
nprev_merged_patches = 0;
for l = 2:L_max
    nrows_merge_idx = (size(Q{L_max - l + 1}, 2) + nnew_merged_patches) / 2;
    merge_idx{2 * l - 1} = zeros(nrows_merge_idx, 2);
    kp = 1; % counter
    r = 1; c = 1; % indices

    % merge_idx{odd}
    i = 1; j = 1;
    while i <= size(Q{L_max - l + 1}, 2) & j <= nnew_merged_patches
        if Q{L_max - l + 1}(i) <= merged_patches(nprev_merged_patches + j)
            merge_idx{2 * l - 1}(r, c) = kp;
            kp = kp + 1;
            i = i + 1;
        else
            merge_idx{2 * l - 1}(r, c) = km;
            km = km - 1;
            j = j + 1;
        end
        if c == 1
            c = 2;
        else
            c = 1; r = r + 1;
        end
    end

    while i <= size(Q{L_max - l + 1}, 2)
        merge_idx{2 * l - 1}(r, c) = kp;
        kp = kp + 1;
        i = i + 1;

        if c == 1
            c = 2;
        else
            c = 1; r = r + 1;
        end
    end

    while j <= nnew_merged_patches
        merge_idx{2 * l - 1}(r, c) = km;
        km = km - 1;
        j = j + 1;

        if c == 1
            c = 2;
        else
            c = 1; r = r + 1;
        end
    end

    % merge_idx{even}
    merge_idx{2 * l} = zeros(nrows_merge_idx / 2, 2);
    merge_idx{2 * l}(:, 2) = 2 * (1:(nrows_merge_idx / 2));
    merge_idx{2 * l}(:, 1) = merge_idx{2 * l}(:, 2) - 1;

    % collect merged patches
    nprev_merged_patches = size(merged_patches, 2);
    nnew_merged_patches = size(merge_idx{2 * l}, 1);
    for i = 1:nnew_merged_patches
        merge_idx_val = merge_idx{2 * l - 1}(2 * i - 1, 1);
        if merge_idx_val > 0
            merged_patches = [merged_patches bitshift(Q{L_max - l + 1}(merge_idx_val), -2)]; 
        else
            merged_patches = [merged_patches bitshift(merged_patches(k2), -2)];
            k2 = k2 + 1;
        end
    end
    % fprintf("merge_idx{%d}\n", 2 * l - 1)
    % disp(merge_idx{2 * l - 1})
    % fprintf("merge_idx{%d}\n", 2 * l)
    % disp(merge_idx{2 * l})
end

end