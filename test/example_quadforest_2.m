function [quadforest,C,L_max] = example_quadforest_2()
%UNTITLED15 Summary of this function goes here
%   Detailed explanation goes here

C = [10 2 4 3;
     11 3 5 1;
     12 1 6 2;
     1 5 7 6;
     2 6 8 4;
     3 4 9 5;
     4 8 10 9;
     5 9 11 7;
     6 7 12 8;
     7 11 1 12;
     8 12 2 10; 
     9 10 3 11];

n_pat = size(C, 1);
quadforest = cell(n_pat, 1);
% L_max = 5;
% L_max = 4;
L_max = 6;
for i = 1:12
    quadforest{i} = cell(L_max,1);
    if i == 5
        % quadforest{5}{1} = uint64([0b00 0b01]);
        % quadforest{5}{2} = uint64([0b1000, 0b1001, 0b1010, 0b1100, 0b1101, 0b1110]);
        % quadforest{5}{3} = uint64([0b101100, 0b101101, 0b101110, 0b101111, 0b111101, 0b111110, 0b111111]);
        % quadforest{5}{4} = uint64([0b11110000, 0b11110001, 0b11110010, 0b11110011]);

        quadforest{5}{1} = uint64([0b01 0b10 0b11]);
        quadforest{5}{2} = uint64([0b0000 0b0001 0b0010]);
        quadforest{5}{3} = uint64([0b001100 0b001110]);
        quadforest{5}{4} = uint64([0b00110100 0b00110110 0b00110111 0b00111100 0b00111101 0b00111110]);
        quadforest{5}{5} = uint64([0b0011010100 0b0011010101 0b0011010110 0b0011010111 0b0011111100 0b0011111101 0b0011111110]);
        quadforest{5}{6} = uint64([0b001111111100 0b001111111101 0b001111111110 0b001111111111]);
    elseif i == 8
        quadforest{i}{1} = uint64([0b00 0b01]);
        quadforest{i}{2} = uint64([0b1000 0b1010 0b1011 0b1100 0b1101 0b1110 0b1111]);
        quadforest{i}{3} = uint64([0b100100 0b100101 0b100110]);
        quadforest{i}{4} = uint64([0b10011100 0b10011101 0b10011110]);
        quadforest{i}{5} = uint64([0b1001111100 0b1001111101 0b1001111110 0b1001111111]);
        quadforest{i}{6} = uint64([]);
    else
        quadforest{i}{1} = uint64([0 1 2 3]);
        for j = 2:L_max
            % ensure each quadtree has the same number of levels, even if some are empty
            quadforest{i}{j} = uint64([]);
        end
    end
end

end