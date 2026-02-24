function plot_quadtree(obj, t, varargin)

if length(varargin) >= 1
    fignum = varargin{1};
else
    fignum = 1;
end

if length(varargin) >= 2
    shift = varargin{2};
else
    shift = [0 0];
end

L_max = size(obj.morton{t},1);

if length(varargin) >= 3
    L_max_forest = varargin{3};
else
    L_max_forest = L_max;
end

figure(fignum)
set(gca, 'YDir','reverse')
set(gcf, 'Color', 'white')
axis off
axis equal
hold on

if ~isempty(obj.morton{t})
    for i = 1:L_max
        for q = obj.morton{t}{i}
            scl = 2^(L_max_forest-i);
            [x, y] = quadforest.deinterleave(q, L_max);
            plot_square([int64(scl*x)+2^L_max_forest*shift(1) int64(scl*y)+2^L_max_forest*shift(2)], 2^(L_max_forest-i))
        end
    end
end

% red bounding box
plot_square([2^L_max_forest*shift(1) 2^L_max_forest*shift(2)], 2^(L_max_forest), 'r--')
% tree label
text(2^L_max_forest*shift(1) + 0.1, 2^L_max_forest*shift(2) + 0.3, int2str(t))

end

function plot_square(top_left, side_length, varargin)

if isempty(varargin)
    opts = 'k-';
else
    opts = varargin{1};
end

x = top_left(1);
y = top_left(2);
s = side_length;
plot([x x+s], [y y], opts)
plot([x x+s], [y+s y+s], opts)
plot([x x], [y y+s], opts)
plot([x+s x+s], [y y+s], opts)

end