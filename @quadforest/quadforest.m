classdef quadforest
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

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
        tree_roots
    end

    methods
        function obj = quadforest(morton, L_max, C, tree_roots)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
            obj.L_max = L_max;
            obj.n_trees = length(morton);
            obj.C = C;
            [obj.morton, obj.addl_patches, obj.remo_patches] = obj.balance_quadforest(morton, tree_roots);
            obj.merge_idx = cell(obj.n_trees, 1);
            obj.merged_patches = cell(obj.n_trees, 1);
            obj.tree_roots = tree_roots;
            % for i = 1:obj.n_trees
            %     [obj.merge_idx{i}, obj.merged_patches{i}] = quadforest.mergeIdxQuadtree(morton{i});
            % end
        end

        [balanced_morton, addl_patches, remo_patches] = balance_quadforest(obj, morton, tree_roots)
        plot_quadtree_merge(obj, elem, level)
        split = get_split(obj, p2q)
        plot_quadtree(obj, t, varargin)
        plot_quadforest(obj, varargin)
        rotated_node = rotate_node(obj, level, node, direction)
        c = forest_colleagues(obj, xy, forest_id, n)
        dir = get_rot_dir(obj, root, root2nei_idx)
    end

    methods(Static)
        [x, y] = deinterleave(xy, n)
        result = interleave(x, y, n)
        parents = get_parents(children)
        children = get_children(parents)
        [merge_idx, merged_patches] = mergeIdxQuadtree(Q)
    end
end