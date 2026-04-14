function surfacemesh_to_vtk(S, fname, f, opts)
%SURFACEMESH_TO_VTK simple surfacefun to vtk file converter for plotting
%   SURFACEMESH_TO_VTK(S, FNAME) writes the surfacemesh S to a VTK file.
%   SURFACEMESH_TO_VTK(S, FNAME, F) also writes the surfacefun F as scalar
%   point data, or the surfacefunv F as vector point data.
%   SURFACEMESH_TO_VTK(..., 'Points', PTS) adds PTS (N x 3) as standalone
%   3D points in the same VTK file.
%   SURFACEMESH_TO_VTK(..., 'Title', T) sets the VTK file title.
arguments (Input)
    S
    fname
    f = []
    opts.Title (1,1) string = "surfacemesh vtk file"
    opts.Points (:,3) double = zeros(0,3)
end
title = opts.Title;
extraPts = opts.Points;

nsurf = length(S);
nsurfpts = 9*nsurf;
nextra = size(extraPts, 1);
npts = nsurfpts + nextra;

fid = fopen(fname, 'w');

fprintf(fid, '# vtk DataFile Version 3.0\n');
fprintf(fid, title);
fprintf(fid, '\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS ');
fprintf(fid, '%d float\n', npts);

% Get the number of points in each direction on each surface patch
n = size(S.x{1}, 1); % Assuming all patches have the same number of points
n2 = idivide(n, int32(2)); % n/2, used for edge points

% Extract and write the coordinates of points to the file
for i = 1:nsurf
    % Four corners of the surface patch
    fprintf(fid, '%f %f %f\n', S.x{i}(1,1), S.y{i}(1,1), S.z{i}(1,1));
    fprintf(fid, '%f %f %f\n', S.x{i}(1,n), S.y{i}(1,n), S.z{i}(1,n));
    fprintf(fid, '%f %f %f\n', S.x{i}(n,n), S.y{i}(n,n), S.z{i}(n,n));
    fprintf(fid, '%f %f %f\n', S.x{i}(n,1), S.y{i}(n,1), S.z{i}(n,1));
    % Four edge points on the surface patch
    fprintf(fid, '%f %f %f\n', S.x{i}(1,n2), S.y{i}(1,n2), S.z{i}(1,n2));
    fprintf(fid, '%f %f %f\n', S.x{i}(n2,n), S.y{i}(n2,n), S.z{i}(n2,n));
    fprintf(fid, '%f %f %f\n', S.x{i}(n,n2), S.y{i}(n,n2), S.z{i}(n,n2));
    fprintf(fid, '%f %f %f\n', S.x{i}(n2,1), S.y{i}(n2,1), S.z{i}(n2,1));
    % Center point of the surface patch
    fprintf(fid, '%f %f %f\n', S.x{i}(n2,n2), S.y{i}(n2,n2), S.z{i}(n2,n2));
end

% Optional standalone points
for i = 1:nextra
    fprintf(fid, '%f %f %f\n', extraPts(i,1), extraPts(i,2), extraPts(i,3));
end

% Write the cell labels to the file
fprintf(fid, 'CELLS %d %d\n', nsurf + nextra, 10*nsurf + 2*nextra);
for i = 1:nsurf
    fprintf(fid, '9 %d %d %d %d %d %d %d %d %d\n', 9*i-9, 9*i-8, 9*i-7, 9*i-6, 9*i-5, 9*i-4, 9*i-3, 9*i-2, 9*i-1);
end
for i = 1:nextra
    fprintf(fid, '1 %d\n', nsurfpts + i - 1);
end

% Write the cell types to the file
fprintf(fid, 'CELL_TYPES %d\n', nsurf + nextra);
for i = 1:nsurf
    fprintf(fid, '28\n'); % Each cell is a quadratic quadrilateral with a center point (type 28)
end
for i = 1:nextra
    fprintf(fid, '1\n'); % VTK_VERTEX for standalone points
end

% Optionally write point data (scalar or vector)
if ~isempty(f)
    fprintf(fid, 'POINT_DATA %d\n', npts);

    if isa(f, 'surfacefun')
        % Scalar field
        fprintf(fid, 'SCALARS scalar_field float 1\n');
        fprintf(fid, 'LOOKUP_TABLE default\n');
        for i = 1:nsurf
            v = f.vals{i};
            % Same sampling order as geometry: 4 corners, 4 edge mids, 1 center
            fprintf(fid, '%f\n', v(1,1));
            fprintf(fid, '%f\n', v(1,n));
            fprintf(fid, '%f\n', v(n,n));
            fprintf(fid, '%f\n', v(n,1));
            fprintf(fid, '%f\n', v(1,n2));
            fprintf(fid, '%f\n', v(n2,n));
            fprintf(fid, '%f\n', v(n,n2));
            fprintf(fid, '%f\n', v(n2,1));
            fprintf(fid, '%f\n', v(n2,n2));
        end
        % No scalar field is defined on optional standalone points.
        for i = 1:nextra
            fprintf(fid, 'nan\n');
        end

    elseif isa(f, 'surfacefunv')
        % Vector field
        fprintf(fid, 'VECTORS vector_field float\n');
        for i = 1:nsurf
            vx = f.components{1}.vals{i};
            vy = f.components{2}.vals{i};
            vz = f.components{3}.vals{i};
            % Same sampling order as geometry: 4 corners, 4 edge mids, 1 center
            idx = {[1,1],[1,n],[n,n],[n,1],[1,n2],[n2,n],[n,n2],[n2,1],[n2,n2]};
            for j = 1:9
                r = idx{j}(1); c = idx{j}(2);
                fprintf(fid, '%f %f %f\n', vx(r,c), vy(r,c), vz(r,c));
            end
        end
        % No vector field is defined on optional standalone points.
        for i = 1:nextra
            fprintf(fid, 'nan nan nan\n');
        end

    else
        error('SURFACEMESH_TO_VTK:invalidField', ...
              'f must be a surfacefun or surfacefunv.');
    end
end

fclose(fid);

end