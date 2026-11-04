function nd = binormal_d(X, TDu, TDv, side)
%BINORMAL_D   Outward binormal derivative along one side of a leaf.
%   ND = SURFACEOP.LEAF.BINORMAL_D(X, TDU, TDV, SIDE) returns the operator that
%   maps the n^2 grid values of a leaf to their derivative along the outward
%   binormal (tangent to the surface, normal to the side) at points on side
%   SIDE = 1, 2, 3, 4 (left, right, down, up). X holds the n^2 x 3 grid
%   coordinates; TDU and TDV evaluate the u- and v-derivatives of grid values
%   at the points. The binormal and metric are evaluated at the points, not
%   interpolated to them: neither is a polynomial, and either may be poorly
%   resolved along a side.

xu = TDu * X;
xv = TDv * X;
E = sum(xu.^2, 2);
F = sum(xu.*xv, 2);
G = sum(xv.^2, 2);
J = E.*G - F.^2;

if ( side <= 2 )
    b = (2*side-3) * xu;
    t = xv;
else
    b = (2*side-7) * xv;
    t = xu;
end
t = t ./ vecnorm(t, 2, 2);
b = b - t .* sum(b.*t, 2);
b = b ./ vecnorm(b, 2, 2);

% b . grad u, with grad u = ((G xu - F xv) u_u + (E xv - F xu) u_v) / J
bu = sum(b.*xu, 2);
bv = sum(b.*xv, 2);
nd = ((G.*bu - F.*bv)./J) .* TDu + ((E.*bv - F.*bu)./J) .* TDv;

end
