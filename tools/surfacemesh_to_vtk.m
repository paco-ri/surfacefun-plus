function surfacemesh_to_vtk(S, fname, varargin)
%SURFACEMESH_TO_VTK simple surfacefun to vtk file converter for plotting
arguments (Input)
    S
    fname
end
arguments (Repeating)
    varargin
end
if isempty(varargin) > 0
    title = "surfacemesh vtk file";
else 
    title = varargin{1};
end
fid = fopen(fname, 'w');

fprintf(fid, '# vtk DataFile Version 3.0\n');
fprintf(fid, title);
fprintf(fid, '\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS ');
fprintf(fid, '%d float\n', 4*length(S));

% Extract and write the coordinates of points to the file
for i = 1:length(S)
    fprintf(fid, '%f %f %f\n', S.x{i}(1,1), S.y{i}(1,1), S.z{i}(1,1));
    fprintf(fid, '%f %f %f\n', S.x{i}(1,end), S.y{i}(1,end), S.z{i}(1,end));
    fprintf(fid, '%f %f %f\n', S.x{i}(end,end), S.y{i}(end,end), S.z{i}(end,end));
    fprintf(fid, '%f %f %f\n', S.x{i}(end,1), S.y{i}(end,1), S.z{i}(end,1));
end

% Write the cell labels to the file
fprintf(fid, 'CELLS %d %d\n', length(S), 5*length(S));
for i = 1:length(S)
    fprintf(fid, '4 %d %d %d %d\n', 4*i-4, 4*i-3, 4*i-2, 4*i-1);
end

% Write the cell types to the file
fprintf(fid, 'CELL_TYPES %d\n', length(S));
for i = 1:length(S)
    fprintf(fid, '9\n'); % Assuming all cells are tetrahedra (type 9)
end

fclose(fid);

end