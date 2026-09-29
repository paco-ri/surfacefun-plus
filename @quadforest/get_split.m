function split = get_split(obj, p2q)
%GET_SPLIT   Get hanging-edge flags for the leaves of a 2:1 balanced quadforest.
%   SPLIT = GET_SPLIT(OBJ, P2Q) returns a cell array with one logical
%   1 x 4 row per leaf, P2Q(i,:) = [tree, level, morton]. SPLIT{i}(s) is true
%   when side s of leaf i, ordered [Left Right Down Up], faces two finer
%   neighbours.

npat = size(p2q, 1);
split = cell(npat, 1);
% The children of cell m are 4*m + k, with k laid out as
%     2 3
%     0 1
% For side s = 1..4 (order L R D U), edgechild(s) = k, where k is an arbitrarily chosen
% child of the neighbor that touches side s.
edgechild = [0 1 0 2];
for i = 1:npat
    t = p2q(i, 1);
    l = p2q(i, 2);
    m = p2q(i, 3);
    split{i} = false(1, 4);
    if l >= obj.L_max
        continue
    end
    [col, back] = obj.forest_colleagues(uint64(m), t, l);
    for s = 1:4
        if col(s, 2) == 0
            continue % surface boundary
        end
        % Side s is a hanging edge iff the same-level neighbour across it is 
        % split. Take one of the neighbour's children that touches leaf i; under 2:1
        % balance the neighbour is split exactly when that child is a leaf.
        child = 4*col(s, 1) + edgechild(back(s));
        split{i}(s) = ismember(uint64(child), obj.morton{col(s, 2)}{l + 1});
    end
end

end
