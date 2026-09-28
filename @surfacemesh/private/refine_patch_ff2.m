function [x_fin,y_fin,z_fin,...
    E_fin,F_fin,G_fin,...
    L_fin,M_fin,N_fin,...
    E_err,F_err,G_err,...
    L_err,M_err,N_err] ...
    = refine_patch_ff2(...
    r,n,x_patch,y_patch,z_patch,...
    E_patch,F_patch,G_patch,...
    L_patch,M_patch,N_patch)
%REFINE_PATCH_FF2 refine a patch based on second fundamental form
%   As REFINE_PATCH_FF1, but the first and second fundamental forms are both
%   returned, along with the discrepancy between interpolating the parent's
%   forms onto each child and recomputing them there.
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
D2 = diffmat(n, 2, [0 1]);
D2 = 4^(r-1) .* D2;
xu_fin = cell(4,1); xv_fin = cell(4,1);
yu_fin = cell(4,1); yv_fin = cell(4,1);
zu_fin = cell(4,1); zv_fin = cell(4,1);
xuu_fin = cell(4,1); xuv_fin = cell(4,1); xvv_fin = cell(4,1);
yuu_fin = cell(4,1); yuv_fin = cell(4,1); yvv_fin = cell(4,1);
zuu_fin = cell(4,1); zuv_fin = cell(4,1); zvv_fin = cell(4,1);
E_fin = cell(4,1); F_fin = cell(4,1); G_fin = cell(4,1);
L_fin = cell(4,1); M_fin = cell(4,1); N_fin = cell(4,1);
E_c2f = cell(4,1); F_c2f = cell(4,1); G_c2f = cell(4,1);
L_c2f = cell(4,1); M_c2f = cell(4,1); N_c2f = cell(4,1);

for i = 1:4
    xu_fin{i} = x_fin{i} * D.'; xv_fin{i} = D * x_fin{i};
    yu_fin{i} = y_fin{i} * D.'; yv_fin{i} = D * y_fin{i};
    zu_fin{i} = z_fin{i} * D.'; zv_fin{i} = D * z_fin{i};

    xuu_fin{i} = x_fin{i} * D2.';
    xuv_fin{i} = D * x_fin{i} * D.';
    xvv_fin{i} = D2 * x_fin{i};
    yuu_fin{i} = y_fin{i} * D2.';
    yuv_fin{i} = D * y_fin{i} * D.';
    yvv_fin{i} = D2 * y_fin{i};
    zuu_fin{i} = z_fin{i} * D2.';
    zuv_fin{i} = D * z_fin{i} * D.';
    zvv_fin{i} = D2 * z_fin{i};

    E_fin{i} = xu_fin{i}.*xu_fin{i} + yu_fin{i}.*yu_fin{i} + zu_fin{i}.*zu_fin{i};
    F_fin{i} = xu_fin{i}.*xv_fin{i} + yu_fin{i}.*yv_fin{i} + zu_fin{i}.*zv_fin{i};
    G_fin{i} = xv_fin{i}.*xv_fin{i} + yv_fin{i}.*yv_fin{i} + zv_fin{i}.*zv_fin{i};

    x_nml_fin = yu_fin{i}.*zv_fin{i} - zu_fin{i}.*yv_fin{i};
    y_nml_fin = zu_fin{i}.*xv_fin{i} - xu_fin{i}.*zv_fin{i};
    z_nml_fin = xu_fin{i}.*yv_fin{i} - yu_fin{i}.*xv_fin{i};
    scl = sqrt(x_nml_fin.^2 + y_nml_fin.^2 + z_nml_fin.^2);
    x_nml_fin = x_nml_fin./scl;
    y_nml_fin = y_nml_fin./scl;
    z_nml_fin = z_nml_fin./scl;

    L_fin{i} = xuu_fin{i}.*x_nml_fin + yuu_fin{i}.*y_nml_fin + zuu_fin{i}.*z_nml_fin;
    M_fin{i} = xuv_fin{i}.*x_nml_fin + yuv_fin{i}.*y_nml_fin + zuv_fin{i}.*z_nml_fin;
    N_fin{i} = xvv_fin{i}.*x_nml_fin + yvv_fin{i}.*y_nml_fin + zvv_fin{i}.*z_nml_fin;
end

% eval coarse-grid E at fine-grid nodes
E_c2f{1} = BL * E_patch * BL.';
F_c2f{1} = BL * F_patch * BL.';
G_c2f{1} = BL * G_patch * BL.';
L_c2f{1} = BL * L_patch * BL.';
M_c2f{1} = BL * M_patch * BL.';
N_c2f{1} = BL * N_patch * BL.';
E_c2f{2} = BL * E_patch * BR.';
F_c2f{2} = BL * F_patch * BR.';
G_c2f{2} = BL * G_patch * BR.';
L_c2f{2} = BL * L_patch * BR.';
M_c2f{2} = BL * M_patch * BR.';
N_c2f{2} = BL * N_patch * BR.';
E_c2f{3} = BR * E_patch * BL.';
F_c2f{3} = BR * F_patch * BL.';
G_c2f{3} = BR * G_patch * BL.';
L_c2f{3} = BR * L_patch * BL.';
M_c2f{3} = BR * M_patch * BL.';
N_c2f{3} = BR * N_patch * BL.';
E_c2f{4} = BR * E_patch * BR.';
F_c2f{4} = BR * F_patch * BR.';
G_c2f{4} = BR * G_patch * BR.';
L_c2f{4} = BR * L_patch * BR.';
M_c2f{4} = BR * M_patch * BR.';
N_c2f{4} = BR * N_patch * BR.';

E_err = cell(4,1); F_err = cell(4,1); G_err = cell(4,1);
L_err = cell(4,1); M_err = cell(4,1); N_err = cell(4,1);
for i = 1:4
    E_err{i} = abs(E_c2f{i} - E_fin{i});
    F_err{i} = abs(F_c2f{i} - F_fin{i});
    G_err{i} = abs(G_c2f{i} - G_fin{i});
    L_err{i} = abs(L_c2f{i} - L_fin{i});
    M_err{i} = abs(M_c2f{i} - M_fin{i});
    N_err{i} = abs(N_c2f{i} - N_fin{i});
end

end
