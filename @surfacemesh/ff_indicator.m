function eta = ff_indicator(dom, p2q, mode)
%FF_INDICATOR   Per-patch fundamental-form refinement indicator.
%   ETA = FF_INDICATOR(DOM, P2Q) returns an NPAT x 1 vector whose entry i is
%   the following error indicator of patch i: interpolate the patch's first
%   fundamental form onto its four children, recompute the form directly on
%   each child, and take the largest patch L^2 norm of the difference over
%   the children and the three coefficients E, F, G. A
%   patch whose ETA is below a tolerance already resolves the surface.
%
%   ETA = FF_INDICATOR(DOM, P2Q, MODE) with MODE = 2 involved the second
%   fundamental form into the indicator as well. The default is MODE = 1.
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
eta = zeros(npat, 1);

for i = 1:npat
    r = p2q(i, 2) + 1;
    scl = 4^(r-1); % scaling so that errors on children are comparable to those on the parent patch
    E = scl*dom.E{i};
    F = scl*dom.F{i};
    G = scl*dom.G{i};

    if ( mode == 2 )
        [L, M, N] = patch_ff2(n, dom.x{i}, dom.y{i}, dom.z{i});
        [~,~,~, E_fin, F_fin, G_fin, ~,~,~, ...
            E_err, F_err, G_err, L_err, M_err, N_err] = ...
            refine_patch_ff2(r, n, dom.x{i}, dom.y{i}, dom.z{i}, ...
            E, F, G, scl*L, scl*M, scl*N);
        err = @(j) {E_err{j} F_err{j} G_err{j} L_err{j} M_err{j} N_err{j}};
    else
        [~,~,~, E_fin, F_fin, G_fin, E_err, F_err, G_err] = ...
            refine_patch_ff1(r, n, dom.x{i}, dom.y{i}, dom.z{i}, E, F, G);
        err = @(j) {E_err{j} F_err{j} G_err{j}};
    end

    for j = 1:4
        J = E_fin{j}.*G_fin{j} - F_fin{j}.^2;
        errj = err(j);
        for k = 1:numel(errj)
            eta(i) = max(eta(i), patchL2norm(errj{k}, J));
        end
    end
end

end
