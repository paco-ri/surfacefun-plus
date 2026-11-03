function dom = taylor_stellarator(n, nu, nv)
%TAYLOR_STELLARATOR The one-surface stellarator of fmm-taylor's PREPARE_STELLARATOR.
%   DOM = TAYLOR_STELLARATOR(N, NU, NV) has NU patches toroidally (u) and NV
%   poloidally (v), each on an N x N Chebyshev grid. It is not
%   SURFACEMESH.STELLARATOR: that one makes v the toroidal angle and uses
%   4 instead of 2 for the major-radius coefficient d(3,2). The hanging-edge
%   tests were first measured on this surface.

ubreaks = linspace(0, 2*pi, nu+1);
vbreaks = linspace(0, 2*pi, nv+1);
x = cell(nu*nv, 1); y = x; z = x;
k = 1;
for ku = 1:nu
    for kv = 1:nv
        [uu, vv] = chebpts2(n, n, [ubreaks(ku:ku+1) vbreaks(kv:kv+1)]);
        Q = 3;
        d = [0.15  0.09  0.00  0.00  0.00;
             0.00  1.00  0.03 -0.01  0.00;
             0.08  2.00 -0.01 -0.02  0.00;
             0.01 -0.28 -0.28  0.03  0.02;
             0.00  0.09 -0.03  0.06  0.00;
             0.01 -0.02  0.02  0.00 -0.02];
        rz = zeros(size(vv));
        for j = -1:4
            for kk = -1:3
                rz = rz + d(j+2, kk+2)*exp(-1i*j*vv + 1i*kk*Q*uu);
            end
        end
        rz = exp(1i*vv).*rz;
        [x{k}, y{k}] = pol2cart(uu, real(rz));
        z{k} = imag(rz);
        k = k + 1;
    end
end
if ( bitand(nv, nv-1) == 0 && bitand(nu, nu-1) == 0 )
    ordering = morton(nv, nu);
    ordering = ordering(:);
    x(ordering) = x;
    y(ordering) = y;
    z(ordering) = z;
end
dom = surfacemesh(x, y, z);
end
