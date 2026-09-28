function [L, M, N] = patch_ff2(n, x, y, z)
%PATCH_FF2 second fundamental form of a single patch
%   [L, M, N] = PATCH_FF2(N, X, Y, Z) returns the second fundamental form of
%   the patch with nodes X, Y, Z, differentiated with respect to the patch's
%   own [-1,1] parameter square.
%
%   Shared by SURFACEMESH/ADAP_REF and SURFACEMESH/FF_INDICATOR.

D  = diffmat(n, 1);
D2 = diffmat(n, 2);

xu = x * D.';  xv = D * x;
yu = y * D.';  yv = D * y;
zu = z * D.';  zv = D * z;

xuu = x * D2.';  xuv = D * x * D.';  xvv = D2 * x;
yuu = y * D2.';  yuv = D * y * D.';  yvv = D2 * y;
zuu = z * D2.';  zuv = D * z * D.';  zvv = D2 * z;

x_nml = yu.*zv - zu.*yv;
y_nml = zu.*xv - xu.*zv;
z_nml = xu.*yv - yu.*xv;
scl = sqrt(x_nml.^2 + y_nml.^2 + z_nml.^2);
x_nml = x_nml./scl;
y_nml = y_nml./scl;
z_nml = z_nml./scl;

L = xuu.*x_nml + yuu.*y_nml + zuu.*z_nml;
M = xuv.*x_nml + yuv.*y_nml + zuv.*z_nml;
N = xvv.*x_nml + yvv.*y_nml + zvv.*z_nml;

end
