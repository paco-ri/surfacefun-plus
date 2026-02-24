function split = get_split(obj, p2q, C)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    obj 
    p2q % patch --> quadforest mapping
    C % connectivity
end

arguments (Output)
    split
end

npat = size(p2q, 1);
split = cell(npat, 1);
for i = 1:npat
    % p2q(i, :) = (tree, level, Morton)
    t = p2q(i, 1);
    l = p2q(i, 2);
    m = p2q(i, 3);
    split{i} = false(1, 4);
    if l == 0 % no tree on this patch
        nei_tree = C(t, 1);
        if ismember(1, obj.morton{nei_tree}{1})
            split{i}(1) = 1;
        end
        nei_tree = C(t, 2);
        if ismember(0, obj.morton{nei_tree}{1})
            split{i}(2) = 1;
        end
        nei_tree = C(t, 3);
        if ismember(0, obj.morton{nei_tree}{1})
            split{i}(3) = 1;
        end
        nei_tree = C(t, 4);
        if ismember(2, obj.morton{nei_tree}{1})
            split{i}(4) = 1;
        end
    else
        [x, y] = quadforest.deinterleave(m, l);

        % left neighbor
        if x > 0
            nei_tree = t; % left neighbor tree
            nei_child = 4 * quadforest.interleave(x - 1, y, l) + 1; % example Morton code one level deeper that borders m
        else % search to tree at left
            nei_tree = C(t, 1);
            nei_child = 4 * quadforest.interleave(2^l - 1, y, l) + 1;
        end
        rot_dir = obj.get_rot_dir(nei_tree, 1);
        if abs(rot_dir) < 2
            nei_child = obj.rotate_node(l + 1, nei_child, rot_dir);
        end
        if l < obj.L_max && ismember(nei_child, obj.morton{nei_tree}{l + 1})
            split{i}(1) = 1;
        end

        % right neighbor
        if x < 2^l - 1
            nei_tree = t;
            nei_child = 4 * quadforest.interleave(x + 1, y, l);
        else
            nei_tree = C(t, 2); 
            nei_child = 4 * quadforest.interleave(0, y, l);
        end
        rot_dir = obj.get_rot_dir(nei_tree, 2);
        if abs(rot_dir) < 2
            nei_child = obj.rotate_node(l + 1, nei_child, rot_dir);
        end
        if l < obj.L_max && ismember(nei_child, obj.morton{nei_tree}{l + 1})
            split{i}(2) = 1;
        end

        % neighbor below
        if y > 0
            nei_tree = t;
            nei_child = 4 * quadforest.interleave(x, y - 1, l);
        else
            if (i == 43)
                disp("hi")
            end
            nei_tree = C(t, 3);
            nei_child = 4 * quadforest.interleave(x, 2^l - 1, l);
        end
        rot_dir = obj.get_rot_dir(nei_tree, 3);
        if abs(rot_dir) < 2
            nei_child = obj.rotate_node(l + 1, nei_child, rot_dir);
        end
        if l < obj.L_max && ismember(nei_child, obj.morton{nei_tree}{l + 1})
            split{i}(3) = 1;
        end

        % neighbor above
        if y < 2^l - 1
            nei_tree = t;
            nei_child = 4 * quadforest.interleave(x, y + 1, l) + 2;
        else
            nei_tree = C(t, 4);
            nei_child = 4 * quadforest.interleave(x, 0, l) + 2;
        end
        rot_dir = obj.get_rot_dir(nei_tree, 4);
        if abs(rot_dir) < 2
            nei_child = obj.rotate_node(l + 1, nei_child, rot_dir);
        end
        if l < obj.L_max && ismember(nei_child, obj.morton{nei_tree}{l + 1})
            split{i}(4) = 1;
        end
    end
    
end

end