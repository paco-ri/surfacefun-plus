function [c, s2] = forest_colleagues(obj, xy, forest_id, n)
%FOREST_COLLEAGUES   Same-level neighbours of a node, across tree boundaries.
%   C = FOREST_COLLEAGUES(OBJ, XY, FOREST_ID, N) returns a 4 x 2 array whose
%   row s is [morton, tree] of the level-N cell adjacent to side s of cell XY
%   in tree FOREST_ID. Sides are ordered [Left Right Down Up], as are the
%   columns of OBJ.C: Left/Right are x-1/x+1 (the u direction), Down/Up are
%   y-1/y+1 (the v direction).
%
%   [C, S2] = FOREST_COLLEAGUES(...) also returns S2, a 4-element vector 
%   that gives the side of the neighboring cell that touches cell XY.
%   Within a tree it is the opposite side, e.g. S2(1) = 2; across trees it 
%   depends on how the two patches are oriented.
%
%   A side on the boundary of an open surface (OBJ.C = -1) gets row [0 0]
%   and S2 = 0.

N = 2^n;
[x, y] = quadforest.deinterleave(xy, n);
x = double(x);
y = double(y);

step = [-1 0; 1 0; 0 -1; 0 1];
inside = [x > 0, x < N-1, y > 0, y < N-1];
opposite = [2 1 4 3];
sgn = [-1 1 1 -1];

c = zeros(4, 2);
s2 = zeros(4, 1);
for s = 1:4
    if inside(s)
        c(s,:) = [double(quadforest.interleave(x + step(s,1), y + step(s,2), n)), forest_id];
        s2(s) = opposite(s);
        continue
    end

    % Across a tree boundary, the neighbor's side comes from OBJ.C.
    t2 = obj.C(forest_id, s);
    if t2 < 1
        continue
    end
    back = find(obj.C(t2, :) == forest_id);
    if numel(back) ~= 1
        error('QUADFOREST:forest_colleagues:ambiguous', ...
            ['Patches %d and %d share %d edges, so the connectivity alone ' ...
             'cannot say which one side %d of patch %d is.'], ...
            forest_id, t2, numel(back), s, forest_id);
    end
    s2(s) = back;

    % The cell's x- or y-coordinate a along the shared edge either remains a
    % or becomes N-1-a. sgn(s) = +1 if a counter-clockwise walk around the
    % patch runs along side s in the direction its parameter increases
    % (Right, Down) and -1 otherwise (Left, Up). Consistently oriented
    % patches walk a shared edge in opposite directions, so a flips exactly
    % when sgn(s)*sgn(back) = 1. This one rule covers all sixteen side
    % pairings.
    if s <= 2
        a = y;
    else
        a = x;
    end
    if sgn(s)*sgn(back) == 1
        a = N - 1 - a;
    end
    switch back
        case 1, xn = 0;     yn = a;
        case 2, xn = N - 1; yn = a;
        case 3, xn = a;     yn = 0;
        case 4, xn = a;     yn = N - 1;
    end
    c(s,:) = [double(quadforest.interleave(xn, yn, n)), t2];
end

end
