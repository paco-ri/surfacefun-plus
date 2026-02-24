function rotated_node = rotate_node(obj, level, node, direction)
% ROTATE_NODE Rotate a node 90 degrees left or right in a quadtree
%
% rotated_node = ROTATE_NODE(obj, level, node, direction)
%   Returns the Morton encoding of the node that is rotated 90 degrees
%   from the given node in the specified direction.
%
% Inputs:
%   obj       - quadforest object
%   level     - The level of the node to rotate (1-indexed)
%   node      - Morton encoding of the node to rotate
%   direction - +1 for clockwise, -1 for counterclockwise, 0 for 180 degrees
%
% Output:
%   rotated_node - Morton encoding of the rotated node
%
% Rotation transformations:
%   Left (counterclockwise):  (x, y) -> (y, max - x)
%   Right (clockwise):        (x, y) -> (max - y, x)
%   180 degrees:              (x, y) -> (max - x, max - y)
%
% where max = 2^level - 1

% Validate level
if level < 1 || level > obj.L_max + 1
    error('Level must be between 1 and %d', obj.L_max);
end

% Validate direction
if direction ~= 1 && direction ~= -1 && direction ~= 0
    error('Direction must be +1 (clockwise), -1 (counterclockwise), or 0 (180 degrees)');
end

% Get the depth of the tree for this level (level 1 has depth 0)
n = level - 1;

% Deinterleave to get x and y coordinates
[x, y] = quadforest.deinterleave(node, n);

% Calculate max coordinate value at this level
max_coord = 2^n - 1;

% Apply rotation based on direction
if direction == -1
    % Left rotation (counterclockwise): (x, y) -> (y, max - x)
    new_x = y;
    new_y = max_coord - x;
elseif direction == 1
    % Right rotation (clockwise): (x, y) -> (max - y, x)
    new_x = max_coord - y;
    new_y = x;
else % direction == 0
    % 180 degree rotation: (x, y) -> (max - x, max - y)
    new_x = max_coord - x;
    new_y = max_coord - y;
end

% Interleave to get Morton encoding of rotated node
rotated_node = quadforest.interleave(new_x, new_y, n);

end
