function [x, y] = deinterleave(xy,n)
x = uint64(0);
y = uint64(0);
for i = 1:n
    x = x + uint64(bitshift(bitand(xy,2^(2*i-2)),1-i));
    y = y + uint64(bitshift(bitand(xy,2^(2*i-1)),-i));
end
end