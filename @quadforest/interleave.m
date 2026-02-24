function result = interleave(x,y,n)
result = uint64(0);
x = uint64(x);
y = uint64(y);
for i = 1:n
    result = result + uint64(bitshift(bitand(x,2^(i-1)),i-1));
    result = result + uint64(bitshift(bitand(y,2^(i-1)),i));
end
end