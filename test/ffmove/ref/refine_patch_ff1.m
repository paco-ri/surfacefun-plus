function [x_fin,y_fin,z_fin,E_fin,F_fin,G_fin,E_err,F_err,G_err] ...
    = refine_patch_ff1(r,n,x_patch,y_patch,z_patch,E_patch,F_patch,G_patch)
%REFINE_PATCH_FF1 refine a patch based on first fundamental form
%   The four children are returned along with their first fundamental forms
%   and, in E_ERR/F_ERR/G_ERR, the discrepancy between interpolating the
%   parent's form onto each child and recomputing it there.
%
%   R is the refinement level: derivatives are taken with respect to the
%   parameters of a level-(R-1) patch, so the forms stay in one fixed unit as
%   the mesh deepens. Pass R = 1 for a base patch.
%
%   Shared by SURFACEMESH/ADAP_REF and SURFACEMESH/FF_INDICATOR.

% tools for refining a patch
x  = chebpts(n, [-1 1]);
xL = chebpts(n, [-1 0]);
xR = chebpts(n, [ 0 1]);
BL = barymat(xL, x);
BR = barymat(xR, x);

% get nodes of fine patches
x_fin = cell(4,1);
y_fin = cell(4,1);
z_fin = cell(4,1);
x_fin{1} = BL * x_patch * BL.';
y_fin{1} = BL * y_patch * BL.';
z_fin{1} = BL * z_patch * BL.';
x_fin{2} = BL * x_patch * BR.';
y_fin{2} = BL * y_patch * BR.';
z_fin{2} = BL * z_patch * BR.';
x_fin{3} = BR * x_patch * BL.';
y_fin{3} = BR * y_patch * BL.';
z_fin{3} = BR * z_patch * BL.';
x_fin{4} = BR * x_patch * BR.';
y_fin{4} = BR * y_patch * BR.';
z_fin{4} = BR * z_patch * BR.';

% get derivatives and forms on fine patch
D = diffmat(n, 1, [0 1]);
D = 2^(r-1) .* D;
xu_fin = cell(4,1); xv_fin = cell(4,1);
yu_fin = cell(4,1); yv_fin = cell(4,1);
zu_fin = cell(4,1); zv_fin = cell(4,1);
E_fin = cell(4,1); F_fin = cell(4,1); G_fin = cell(4,1);
E_c2f = cell(4,1); F_c2f = cell(4,1); G_c2f = cell(4,1);

for i = 1:4
    xu_fin{i} = x_fin{i} * D.'; xv_fin{i} = D * x_fin{i};
    yu_fin{i} = y_fin{i} * D.'; yv_fin{i} = D * y_fin{i};
    zu_fin{i} = z_fin{i} * D.'; zv_fin{i} = D * z_fin{i};

    E_fin{i} = xu_fin{i}.*xu_fin{i} + yu_fin{i}.*yu_fin{i} + zu_fin{i}.*zu_fin{i};
    F_fin{i} = xu_fin{i}.*xv_fin{i} + yu_fin{i}.*yv_fin{i} + zu_fin{i}.*zv_fin{i};
    G_fin{i} = xv_fin{i}.*xv_fin{i} + yv_fin{i}.*yv_fin{i} + zv_fin{i}.*zv_fin{i};
end

% eval coarse-grid E at fine-grid nodes
E_c2f{1} = BL * E_patch * BL.';
F_c2f{1} = BL * F_patch * BL.';
G_c2f{1} = BL * G_patch * BL.';
E_c2f{2} = BL * E_patch * BR.';
F_c2f{2} = BL * F_patch * BR.';
G_c2f{2} = BL * G_patch * BR.';
E_c2f{3} = BR * E_patch * BL.';
F_c2f{3} = BR * F_patch * BL.';
G_c2f{3} = BR * G_patch * BL.';
E_c2f{4} = BR * E_patch * BR.';
F_c2f{4} = BR * F_patch * BR.';
G_c2f{4} = BR * G_patch * BR.';

E_err = cell(4,1); F_err = cell(4,1); G_err = cell(4,1);
for i = 1:4
    E_err{i} = abs(E_c2f{i} - E_fin{i});
    F_err{i} = abs(F_c2f{i} - F_fin{i});
    G_err{i} = abs(G_c2f{i} - G_fin{i});
end

end
