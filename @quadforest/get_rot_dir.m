function dir = get_rot_dir(obj, root, root2nei_idx)
    %TODO docs
    % Get the rotation between root and nei from perspective of root

    C = obj.C;
    nei = C(root, root2nei_idx);
    nei2root_idx = find(C(nei, :) == root);
    if (nei2root_idx == root2nei_idx) % root and nei are both on same side from each other's POV
        dir = 0;
    elseif (root2nei_idx == 1 && nei2root_idx == 3 || root2nei_idx == 2 && nei2root_idx == 4 || root2nei_idx == 3 && nei2root_idx == 2 || root2nei_idx == 4 && nei2root_idx == 1)
        dir = 1; % clockwise
    elseif (root2nei_idx == 1 && nei2root_idx == 4 || root2nei_idx == 2 && nei2root_idx == 3 || root2nei_idx == 3 && nei2root_idx == 1 || root2nei_idx == 4 && nei2root_idx == 2)
        dir = -1; % counterclockwise
    else
        dir = 2; % don't rotate
    end
end