function I = patchL2norm(fvals, J)
%PATCHL2NORM   L^2 norm of a function over a single surface patch.
%   I = PATCHL2NORM(FVALS, J) returns the L^2 norm of the function whose
%   values at the tensor-product second-kind Chebyshev nodes of a patch are
%   FVALS.
%
%   J is the Jacobian determinant; sqrt(J) is the area element.
%
%   FVALS and J are n x n matrices at the same nodes.

[nv, nu] = size(fvals);
wu = chebtech2.quadwts(nu); wu = wu(:);
wv = chebtech2.quadwts(nv); wv = wv(:);
I = sqrt(sum(sum(abs(fvals).^2 .* wv .* wu.' .* sqrt(J))));

end
