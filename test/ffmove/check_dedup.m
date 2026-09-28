n = 5; nv = 3; nu = 3*nv;
st = prepare_stellarator(n, nu, nv, 16, 40);
sh = prepare_torus(n, nu, nv, n, nu, nv, 1.0, 0.6, 16, 40, 40);
doms = {st{1}, sh{1}, sh{2}};
names = {'stellarator', 'torus outer', 'torus inner'};
worst = 0;
for d = 1:3
    dom0 = doms{d};
    npat0 = length(dom0.x);
    dom = dom0; p2q = [(1:npat0).', zeros(npat0, 2)];
    for lev = 0:2
        if lev > 0
            [dom, ~, p2q] = surfacemesh.refine_leaves(dom0, p2q, 1:size(p2q,1), 3);
        end
        for mode = 1:2
            a = ff_indicator_old(dom, p2q, mode);
            b = surfacemesh.ff_indicator(dom, p2q, mode);
            rel = max(abs(a - b) ./ abs(a));
            worst = max(worst, rel);
            fprintf('%-12s lev %d mode %d: %4d patches, max eta %.4e, max rel diff %.1e\n', ...
                names{d}, lev, mode, numel(a), max(a), rel);
        end
    end
end
fprintf('WORST rel diff %.2e\n', worst);
disp('DONE_SENTINEL')
