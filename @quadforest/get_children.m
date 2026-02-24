function children = get_children(parents)
i = 1;
children = zeros(1,4*length(parents));
for p = parents
    shift_p = bitshift(p, 2);
    children(i) = shift_p;
    children(i+1) = shift_p + 1;
    children(i+2) = shift_p + 2;
    children(i+3) = shift_p + 3;
    i = i + 4;
end
end