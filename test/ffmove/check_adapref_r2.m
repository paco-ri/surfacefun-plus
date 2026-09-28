% rmax = 2 with thresholds low enough that level-1 leaves get marked.
n = 5; nv = 3; nu = 3*nv;
st = prepare_stellarator(n, nu, nv, 16, 40);
sh = prepare_torus(n, nu, nv, n, nu, nv, 1.0, 0.6, 16, 40, 40);
doms = {st{1}, sh{1}};
names = {'stellarator', 'torus outer'};
for d = 1:2
    dom0 = doms{d}; np0 = length(dom0.x);
    for mode = 1:2
        [dom1, ~, p1] = surfacemesh.refine_leaves(dom0, [(1:np0).', zeros(np0,2)], 1:np0, 2);
        eta1 = surfacemesh.ff_indicator(dom1, p1, mode);
        for q = [0.5 0.8 0.95]
            tol = quantile(eta1, q);
            evalc('[~, ~, pold] = adap_ref_old(dom0, tol, 2, mode, 1:np0);');
            [~, ~, pnew, h] = surfacemesh.adap_ref(dom0, tol, 2, mode);
            onlyold = setdiff(pold, pnew, 'rows');
            onlynew = setdiff(pnew, pold, 'rows');
            fprintf('%-11s mode %d tol %.2e: old %3d leaves (depth-2 %3d), new %3d (depth-2 %3d); only-old %d only-new %d\n', ...
                names{d}, mode, tol, size(pold,1), sum(pold(:,2)==2), size(pnew,1), sum(pnew(:,2)==2), ...
                size(onlyold,1), size(onlynew,1));
            if ~isempty(onlyold) || ~isempty(onlynew)
                % Explain: is every depth-2 leaf only in new a child of a
                % level-1 leaf with eta > tol that old never re-examined?
                h1 = h(2); key = @(P) P(:,1)*1e6 + P(:,2)*1e3 + P(:,3);
                par = [onlynew(:,1), onlynew(:,2)-1, floor(onlynew(:,3)/4)];
                par = unique(par(par(:,2)==1,:), 'rows');
                [tf, loc] = ismember(key(par), key(h1.p2q));
                fprintf('    only-new depth-2 parents at level 1: %d, found %d, with eta > tol %d\n', ...
                    size(par,1), sum(tf), sum(h1.eta(loc(tf)) > tol));
            end
        end
    end
end
disp('DONE_SENTINEL')
