classdef quadforest
    %QUADFOREST   A forest of quadtrees, one per patch of a base mesh.
    %   QF = QUADFOREST(MORTON, L_MAX, C) stores the leaves of each tree as
    %   Morton codes: MORTON{t}{l} lists the level-l leaves of tree t, and a
    %   tree with no codes is an unrefined base patch. C is the base mesh's
    %   element-to-element connectivity, ordered [Left Right Down Up]. At 
    %   construction, the forest is not necessarily 2:1 balanced. call 
    %   BALANCE_QUADFOREST to make it 2:1 balanced.

    properties ( Access = public )
        L_max
        morton % Morton encodings
        C % element-to-element connectivity
        n_trees = 0
        merge_idx
        merged_patches
        split
        addl_patches % patches added to level-restrict forest
        remo_patches % patches removed to level-restrict forest
    end

    properties ( Dependent )
        tree_roots % trees with at least one Morton code
    end

    methods
        function obj = quadforest(morton, L_max, C)
            obj.L_max = L_max;
            obj.n_trees = length(morton);
            obj.C = C;
            obj.morton = morton;
            obj.merge_idx = cell(obj.n_trees, 1);
            obj.merged_patches = cell(obj.n_trees, 1);
            % for i = 1:obj.n_trees
            %     [obj.merge_idx{i}, obj.merged_patches{i}] = quadforest.mergeIdxQuadtree(morton{i});
            % end
        end

        function r = get.tree_roots(obj)
            r = find(cellfun(@(q) iscell(q) && any(~cellfun(@isempty, q)), ...
                obj.morton)).';
        end

        obj = balance_quadforest(obj)
        plot_quadtree_merge(obj, elem, level)
        split = get_split(obj, p2q)
        plot_quadtree(obj, t, varargin)
        plot_quadforest(obj, varargin)
        rotated_node = rotate_node(obj, level, node, direction)
        [c, s2] = forest_colleagues(obj, xy, forest_id, n)
    end

    methods(Static)
        [x, y] = deinterleave(xy, n)
        result = interleave(x, y, n)
        parents = get_parents(children)
        children = get_children(parents)
        [merge_idx, merged_patches] = mergeIdxQuadtree(Q)
    end
end
