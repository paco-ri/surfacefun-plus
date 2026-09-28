% New adap_ref vs the committed one (adap_ref_old), threshold marking.
n = 5; nv = 3; nu = 3*nv;
st = prepare_stellarator(n, nu, nv, 16, 40);
sh = prepare_torus(n, nu, nv, n, nu, nv, 1.0, 0.6, 16, 40, 40);
doms = {st{1}, sh{1}};
names = {'stellarator', 'torus outer'};
for d = 1:2
    dom0 = doms{d}; np0 = length(dom0.x);
    eta0 = surfacemesh.ff_indicator(dom0, [(1:np0).', zeros(np0,2)], 1);
    for mode = 1:2
        eta0 = surfacemesh.ff_indicator(dom0, [(1:np0).', zeros(np0,2)], mode);
        for q = [0.5 0.8]
            tol = quantile(eta0, q);
            for rmax = 1:2
                evalc('[~, ~, pold] = adap_ref_old(dom0, tol, rmax, mode, 1:np0);');
                [~, ~, pnew] = surfacemesh.adap_ref(dom0, tol, rmax, mode);
                onlyold = size(setdiff(pold, pnew, 'rows'), 1);
                onlynew = size(setdiff(pnew, pold, 'rows'), 1);
                fprintf('%-11s mode %d q %.1f rmax %d: old %3d leaves, new %3d; only-old %d only-new %d\n', ...
                    names{d}, mode, q, rmax, size(pold,1), size(pnew,1), onlyold, onlynew);
            end
        end
    end
    % depth 3: old breaks, new must complete and solve
    for mk = {5, [5 14], 1:np0}
        for rmax = 3
            [domR, ~, p2q, hist] = surfacemesh.adap_ref(dom0, 0, rmax, 1, mk{1});
            ok = 'solve OK';
            try
                vn = normal(domR);
                phihat = surfacefunv(@(x,y,z) -y./sqrt(x.^2+y.^2), ...
                    @(x,y,z) x./sqrt(x.^2+y.^2), @(x,y,z) 0.*z, domR);
                [~, ~, vH] = hodge(cross(vn, phihat)); %#ok<ASGLU>
            catch ME
                ok = ['SOLVE FAIL: ' ME.message];
            end
            outside = sum(~ismember(p2q(p2q(:,2) > 0, 1), mk{1}));
            fprintf('%-11s rmax 3 tol 0 marked %-10s: %4d leaves, depth %d, %d passes, %d refined leaves outside MARKED roots (balance), %s\n', ...
                names{d}, mat2str(mk{1}(1:min(end,3))), size(p2q,1), max(p2q(:,2)), numel(hist)-1, outside, ok);
        end
    end
end
disp('DONE_SENTINEL')
