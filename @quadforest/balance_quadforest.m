function obj = balance_quadforest(obj)
%BALANCE_QUADFOREST   Refine a quadforest until it is 2:1 balanced.
%   QF = BALANCE_QUADFOREST(QF) splits leaves of QF.MORTON until no two
%   leaves that share an edge, within a tree or across trees, differ by more
%   than one level. Leaves are only ever split, never merged. Also sets
%   QF.ADDL_PATCHES and QF.REMO_PATCHES.

L_max = obj.L_max;
balanced_morton = obj.morton;
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
    for k = 1:n_pat % loop over roots
        if length(balanced_morton{k}) >= r % quadtree has at least r levels
            n_nodes = length(balanced_morton{k}{r}); % number of nodes in tree k at level r
            for i = 1:n_nodes % loop over these nodes
                parent = bitshift(balanced_morton{k}{r}(i),-2); % get parent of node
                colls = obj.forest_colleagues(parent, k, r-1); % get morton code and tree number of parent's colleagues
                for j = 1:4 % loop over colleages
                    if colls(j,2) == 0
                        continue % surface boundary
                    end
                    % append morton code of colleague to coll_of_par{k}
                    coll_of_par{colls(j,2)}(end+1, 1) = colls(j,1);
                end
            end
        end
    end
    for k = 1:n_pat
        coll_of_par{k} = unique(coll_of_par{k});
    end
    % coll_of_par{k} lists the level-(r-1) cells that tree k must contain.
    % If tree k is still an unrefined base patch and r-1 >= 1, split it into
    % its four level-1 children.
    new_tree_roots = [];
    if r > 1
        for k = 1:n_pat
            % If tree k has required cells at level r-1 and is not yet refined, split it.
            if ~isempty(coll_of_par{k}) && all(cellfun(@isempty, balanced_morton{k}))
                new_tree_roots(end+1) = k; %#ok<AGROW>
            end
        end
    end

    % add trees at new tree roots
    if r > 1 && ~isempty(new_tree_roots)
        for root = new_tree_roots
            balanced_morton{root}{1} = uint64([0b00 0b01 0b10 0b11]);
            addl_morton{root}{1} = uint64([0b00 0b01 0b10 0b11]);
        end
    end

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

obj.morton = balanced_morton;
obj.addl_patches = addl_morton;
obj.remo_patches = remo_morton;

end
