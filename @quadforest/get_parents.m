function parents = get_parents(children)
num_chil = length(children);
parents = zeros(num_chil,1);
for i = 1:num_chil
    parents(i) = bitshift(children(i),-2);
end
parents = unique(parents.');
end