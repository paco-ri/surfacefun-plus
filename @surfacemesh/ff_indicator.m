function eta = ff_indicator(dom, p2q, mode)
%FF_INDICATOR   Per-patch fundamental-form refinement indicator.
%   ETA = FF_INDICATOR(DOM, P2Q) returns an NPAT x 1 vector whose entry i is
%   the following error indicator of patch i: interpolate the patch's first
%   fundamental form onto its four children, recompute the form directly on
%   each child, and take the largest patch L^2 norm of the difference over
%   the children and the three coefficients E, F, G. A
%   patch whose ETA is below a tolerance already resolves the surface.
%
%   ETA = FF_INDICATOR(DOM, P2Q, MODE) with MODE = 2 includes the second
%   fundamental form L, M, N as well. The default is MODE = 1.
%
%   Forms are taken with respect to the parameters of the patch's base
%   (level-0) patch, so ETA stays in one unit as the mesh deepens.
%
%   P2Q is the leaf list of DOM, as returned by REFINE_LEAVES or ADAP_REF:
%   row i is [tree_root, level, morton] for patch i.
%
%   See also SURFACEMESH/ADAP_REF, SURFACEMESH/REFINE_LEAVES.

arguments (Input)
    dom
    p2q  (:,3) double
    mode (1,1) double = 1
end

arguments (Output)
    eta
end

npat = length(dom.x);
if ( size(p2q, 1) ~= npat )
    error('SURFACEMESH:ff_indicator:p2q', ...
        'P2Q has %d rows but DOM has %d patches.', size(p2q, 1), npat);
end

n = size(dom.x{1}, 1);
D  = diffmat(n, 1);
D2 = diffmat(n, 2);
t = chebpts(n, [-1 1]);
B = {barymat(chebpts(n, [-1 0]), t), barymat(chebpts(n, [0 1]), t)};

eta = zeros(npat, 1);
for i = 1:npat
    r = p2q(i, 2) + 1;
    x = dom.x{i}; y = dom.y{i}; z = dom.z{i};
    ff = patch_forms(x, y, z, mode, D, D2);
    for jv = 1:2
        for ju = 1:2
            c = @(f) B{jv} * f * B{ju}.';
            % A child's own parameters are half its parent's.
            ffc = patch_forms(c(x), c(y), c(z), mode, D, D2);
            J = 16^r * (ffc{1}.*ffc{3} - ffc{2}.^2);
            for k = 1:numel(ff)
                err = abs(4^(r-1) * c(ff{k}) - 4^r * ffc{k});
                eta(i) = max(eta(i), patchL2norm(err, J));
            end
        end
    end
end

end

function ff = patch_forms(x, y, z, mode, D, D2)
%PATCH_FORMS   {E, F, G} of a patch, plus {L, M, N} if MODE == 2, with
%   respect to the patch's own [-1,1]^2 parameters.

xu = x * D.';  xv = D * x;
yu = y * D.';  yv = D * y;
zu = z * D.';  zv = D * z;
ff = {xu.*xu + yu.*yu + zu.*zu, ...
      xu.*xv + yu.*yv + zu.*zv, ...
      xv.*xv + yv.*yv + zv.*zv};
if ( mode ~= 2 )
    return
end

nx = yu.*zv - zu.*yv;
ny = zu.*xv - xu.*zv;
nz = xu.*yv - yu.*xv;
s = sqrt(nx.^2 + ny.^2 + nz.^2);
nx = nx./s; ny = ny./s; nz = nz./s;
ff(4:6) = {(x * D2.').*nx     + (y * D2.').*ny     + (z * D2.').*nz, ...
           (D * x * D.').*nx  + (D * y * D.').*ny  + (D * z * D.').*nz, ...
           (D2 * x).*nx       + (D2 * y).*ny       + (D2 * z).*nz};

end
