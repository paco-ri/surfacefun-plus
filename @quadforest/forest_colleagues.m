function c = forest_colleagues(obj, xy, forest_id, n)
[x, y] = quadforest.deinterleave(xy,n);

c = [0 0;
     0 0;
     0 0;
     0 0];

if x > 0
    c(1,1) = quadforest.interleave(x-1,y,n);
    c(1,2) = forest_id;
else % search to the left
    c(1,2) = obj.C(forest_id, 1);
    % look at where forest_id is from POV of left neighbor
    % find position of forest_id in C(c(1,2), :)
    neighbor_idx = find(obj.C(c(1,2), :) == forest_id);
    c(1,1) = quadforest.interleave(2^n-1,y,n);
    if neighbor_idx == 1 % forest_id and neighbor are both on left of each other
        c(1,1) = obj.rotate_node(n+1,c(1,1),0); % 180 degrees
    elseif neighbor_idx == 3 % from neighbor POV, forest_id is below
        c(1,1) = obj.rotate_node(n+1,c(1,1),1);
    elseif neighbor_idx == 4 % from neighbor POV, forest_id is above
        c(1,1) = obj.rotate_node(n+1,c(1,1),-1);
    end     
    % rotate
end

if x < 2^n-1 
    c(2,1) = quadforest.interleave(x+1,y,n);
    c(2,2) = forest_id;
else % search to the right
    c(2,2) = obj.C(forest_id, 2);
    % look at where forest_id is from POV of right neighbor
    neighbor_idx = find(obj.C(c(2,2), :) == forest_id);
    c(2,1) = quadforest.interleave(0,y,n);
    if neighbor_idx == 2 % forest_id and neighbor are both on right of each other
        c(2,1) = obj.rotate_node(n+1,c(2,1),0); % 180 degrees
    elseif neighbor_idx == 3 % from neighbor POV, forest_id is below
        c(2,1) = obj.rotate_node(n+1,c(2,1),-1);
    elseif neighbor_idx == 4 % from neighbor POV, forest_id is above
        c(2,1) = obj.rotate_node(n+1,c(2,1),1);
    end
end

if y < 2^n-1
    c(3,1) = quadforest.interleave(x,y+1,n);
    c(3,2) = forest_id;
else % search down
    c(3,2) = obj.C(forest_id, 3);
    % look at where forest_id is from POV of bottom neighbor
    neighbor_idx = find(obj.C(c(3,2), :) == forest_id);
    c(3,1) = quadforest.interleave(x,0,n);
    if neighbor_idx == 3 % forest_id and neighbor are both below each other
        c(3,1) = obj.rotate_node(n+1,c(3,1),0); % 180 degrees
    elseif neighbor_idx == 1 % from neighbor POV, forest_id is to the left
        c(3,1) = obj.rotate_node(n+1,c(3,1),-1);
    elseif neighbor_idx == 2 % from neighbor POV, forest_id is to the right
        c(3,1) = obj.rotate_node(n+1,c(3,1),1);
    end
end

if y > 0
    c(4,1) = quadforest.interleave(x,y-1,n);
    c(4,2) = forest_id;
else % search up
    c(4,2) = obj.C(forest_id, 4);
    % look at where forest_id is from POV of top neighbor
    neighbor_idx = find(obj.C(c(4,2), :) == forest_id);
    c(4,1) = quadforest.interleave(x,2^n-1,n);
    if neighbor_idx == 4 % forest_id and neighbor are both above each other
        c(4,1) = obj.rotate_node(n+1,c(4,1),0); % 180 degrees
    elseif neighbor_idx == 1 % from neighbor POV, forest_id is to the left
        c(4,1) = obj.rotate_node(n+1,c(4,1),1);
    elseif neighbor_idx == 2 % from neighbor POV, forest_id is to the right
        c(4,1) = obj.rotate_node(n+1,c(4,1),-1);
    end
end

end