function [balanced_morton, addl_morton, remo_morton] = balance_quadforest(obj, morton, tree_roots)
%BALANCE_QUADFOREST
% Q = quadforest
% L_max = max level of quadforest
% C = edge connectivity

% utility when deciding to plant a new quadtree
temp = [1 3; 0 2; 0 1; 2 3];

L_max = obj.L_max;
balanced_morton = morton;
addl_morton = cell(obj.n_trees, 1);
remo_morton = cell(obj.n_trees, 1);
for i = 1:obj.n_trees
    addl_morton{i} = cell(L_max, 1);
    remo_morton{i} = cell(L_max, 1);
end
n_pat = size(obj.C, 1);
r = L_max; % current refinement level
while r > 0
    % for r-level leaves, get Morton codes of parents' colleagues
    coll_of_par = cell(n_pat,1);
    new_tree_roots = []; % patches that don't have a quadtree
    for k = 1:n_pat % loop over roots
        if length(balanced_morton{k}) >= r % quadtree has at least r levels
            n_nodes = length(balanced_morton{k}{r}); % number of nodes in tree k at level r
            if isempty(coll_of_par{k})
                coll_of_par{k} = nan(4*n_nodes, 1);
            end
            for i = 1:n_nodes % loop over these nodes
                parent = bitshift(balanced_morton{k}{r}(i),-2); % get parent of node
                colls = obj.forest_colleagues(parent, k, r-1); % get morton code and tree number of parent's colleagues
                for j = 1:4 % loop over colleages
                    if isempty(coll_of_par{colls(j,2)})
                        coll_of_par{colls(j,2)} = nan(4*n_nodes, 1);
                    end
                    coll_of_par{colls(j,2)}(4*(i-1)+j) = colls(j,1); % add colleages to coll_of_par in appropriate spot
                end
                new_tree_roots = union(new_tree_roots, colls(:,2).');
            end
        end
    end
    for k = 1:n_pat
        coll_of_par{k} = unique(coll_of_par{k});
        coll_of_par{k} = coll_of_par{k}(~isnan(coll_of_par{k}));
    end
    new_tree_roots = setdiff(new_tree_roots, tree_roots);
    new_tree_roots = reshape(new_tree_roots, 1, []);
    
    to_uproot = nan(1, n_pat);
    for root = new_tree_roots % only plant tree if neighbors have 2+ levels
        j = 1;
        uproot = true;
        while uproot && j <= 4
            nei = obj.C(root, j); % TODO only check relevant neighbors
            one_side_more_than_2 = false;
            nei_has_2_levels = false;
            for i = 2:L_max
                nei_has_2_levels = nei_has_2_levels || ~isempty(balanced_morton{nei}{i}); % true if tree at nei has level i >= 2
            end
            nei_idx = find(obj.C(nei, :) == root); % position of root from nei's POV
            if (nei_idx == j) % j and nei are both on same side from each other's POV
                rot_dir = 0; % 180 degrees
            elseif (j == 1 && nei_idx == 3 || j == 2 && nei_idx == 4 || j == 3 && nei_idx == 2 || j == 4 && nei_idx == 1)
                rot_dir = 1; % clockwise
            elseif (j == 1 && nei_idx == 4 || j == 2 && nei_idx == 3 || j == 3 && nei_idx == 1 || j == 4 && nei_idx == 2)
                rot_dir = -1; % counterclockwise
            else
                rot_dir = 2; % don't rotate
            end
            if rot_dir ~= 2
                levels_on_correct_side = ~ismember(obj.rotate_node(r, temp(j, 1), rot_dir), balanced_morton{nei}{1}); % when j = 1, consider left nei
                % true if nei missing level 1, morton 1, which is a node on nei's right side
                levels_on_correct_side = levels_on_correct_side || ~ismember(obj.rotate_node(r, temp(j, 2), rot_dir), balanced_morton{nei}{1}); 
                % when j = 1, same question for nei level 1, morton 3 
                % i.e., true if nei has 2+ levels on side bordering root
            else
                levels_on_correct_side = ~ismember(temp(j, 1), balanced_morton{nei}{1}); % when j = 1, consider left nei
                levels_on_correct_side = levels_on_correct_side || ~ismember(temp(j, 2), balanced_morton{nei}{1}); 
            end
            levels_border_root = nei_has_2_levels && levels_on_correct_side; % true if nei has level i >= 2 on side bordering root
            one_side_more_than_2 = one_side_more_than_2 || levels_border_root; % true if any level i >= 2 on side bordering root
            all_less_than_2 = ~one_side_more_than_2; % true if nei has < 2 levels on all sides bordering root
            uproot = uproot && all_less_than_2;
            j = j + 1;
        end
        if uproot
            to_uproot(root) = root;
        end
    end
    to_uproot = rmmissing(to_uproot);
    new_tree_roots = setdiff(new_tree_roots, to_uproot);

    % add trees at new tree roots
    if r > 1 && ~isempty(new_tree_roots)
        for root = new_tree_roots
            balanced_morton{root}{1} = uint64([0b00 0b01 0b10 0b11]);
            addl_morton{root}{1} = uint64([0b00 0b01 0b10 0b11]);
        end
    end
    tree_roots = union(tree_roots, new_tree_roots);

    % if nodes in par_of_coll{k} are not in quadforest{k},
    % refine until they are.
    % pass from level r-1 to level 1
    absent_parents = cell(n_pat,1);
    % ^ after for loop, the parents of colleagues not in forest
    for k = 1:n_pat
        if r > 1 && length(balanced_morton{k}) >= r-1
            absent_parents{k} = setdiff(coll_of_par{k}, balanced_morton{k}{r-1});
        else
            absent_parents{k} = coll_of_par{k};
        end

        % must also remove members of coll_of_par with descendants in
        % forest
        j = 1;
        for i = r:L_max
            if length(balanced_morton{k}) >= i
                j_ancestors = balanced_morton{k}{i};
                for l = 1:j
                    j_ancestors = quadforest.get_parents(j_ancestors);
                end
                absent_parents{k} = setdiff(absent_parents{k}, j_ancestors);
                j = j + 1;
            end
        end
    end
    if isempty(absent_parents)
        break;
    end

    absent_ancestors = cell(n_pat,r-1);
    % ^ after while loop, the ancestors of nodes in par_of_coll 
    %   absent from forest.
    for k = 1:n_pat
        level = r - 1;
        while ~isempty(absent_parents{k}) && level > 1
            absent_ancestors{k,level} = absent_parents{k};
            parents_of_absent_parents = quadforest.get_parents(absent_parents{k});
            absent_parents{k} = setdiff(parents_of_absent_parents, balanced_morton{k}{level-1});
            level = level - 1;
        end
        % parents_of_absent_parents = first present ancestors
        % level = level of these ancestors
    
        for i = level+1:r-1
            p = quadforest.get_parents(absent_ancestors{k,i});
            c = quadforest.get_children(p);
            addl_morton{k}{i-1} = unique([addl_morton{k}{i-1} quadforest.get_parents(absent_ancestors{k,i})]);
            addl_morton{k}{i} = unique([addl_morton{k}{i} c]);
            balanced_morton{k}{i-1} = setdiff(balanced_morton{k}{i-1}, quadforest.get_parents(absent_ancestors{k,i}));
            addl_morton{k}{i-1} = setdiff(addl_morton{k}{i-1}, quadforest.get_parents(absent_ancestors{k,i}));
            remo_morton{k}{i} = [remo_morton{k}{i} setdiff(c, balanced_morton{k}{i})];
            balanced_morton{k}{i} = unique([balanced_morton{k}{i} c]);
        end
    end
    r = r - 1;
end

end