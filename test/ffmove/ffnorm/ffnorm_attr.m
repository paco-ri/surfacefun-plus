% Split the 606d1e0 indicator change into its two parts on the level-0 meshes.
n = 5; nv = 3; nu = 3*nv;
st = prepare_stellarator(n, nu, nv, 16, 40);
sh = prepare_torus(n, nu, nv, n, nu, nv, 1.0, 0.6, 16, 40, 40);
doms = {st{1}, sh{1}, sh{2}};
names = {'stellarator', 'shell outer', 'shell inner'};
for d = 1:3
    dom = doms{d};
    npat = length(dom.x);
    eta = zeros(npat, 6);   % [old  E-only  E F G  eachE eachF eachG]
    for i = 1:npat
        [~,~,~, Ef, Ff, Gf, Ee, Fe, Ge] = refine_patch_ff1(1, n, ...
            dom.x{i}, dom.y{i}, dom.z{i}, dom.E{i}, dom.F{i}, dom.G{i});
        for j = 1:4
            J = Ef{j}.*Gf{j} - Ff{j}.^2;
            for c = 1:n
                eta(i,1) = max(eta(i,1), patchL2norm(Ee{j}(:,c), J));
            end
            e = [patchL2norm(Ee{j}, J), patchL2norm(Fe{j}, J), patchL2norm(Ge{j}, J)];
            eta(i,2) = max(eta(i,2), e(1));
            eta(i,3) = max(eta(i,3), max(e));
            eta(i,4:6) = max(eta(i,4:6), e);
        end
    end
    cur = surfacemesh.ff_indicator(dom, [(1:npat).', zeros(npat,2)], 1);
    [~, argmax] = max(eta(:,4:6), [], 2);
    fprintf('\n== %s (%d patches)\n', names{d}, npat);
    fprintf('max eta: old %.4e | fixed norm, E only %.4e | fixed norm, E F G %.4e\n', max(eta(:,1:3)));
    fprintf('ff_indicator now: %.4e (diff vs E F G col %.1e)\n', max(cur), max(abs(cur - eta(:,3))));
    fprintf('max per coef: E %.4e  F %.4e  G %.4e\n', max(eta(:,4:6)));
    fprintf('patches whose largest term is E/F/G: %d / %d / %d\n', sum(argmax==1), sum(argmax==2), sum(argmax==3));
    fprintf('ratio fixed(EFG)/old: min %.2f median %.2f max %.2f\n', ...
        min(eta(:,3)./eta(:,1)), median(eta(:,3)./eta(:,1)), max(eta(:,3)./eta(:,1)));
    fprintf('ratio fixed(E)/old:   min %.2f median %.2f max %.2f\n', ...
        min(eta(:,2)./eta(:,1)), median(eta(:,2)./eta(:,1)), max(eta(:,2)./eta(:,1)));
    % Dorfler sets: does the change alter which patches get marked?
    for th = [0.5 0.6]
        mk = cell(1,3);
        for v = 1:3
            [es, ix] = sort(eta(:,v).^2, 'descend');
            k = find(cumsum(es) >= th*sum(es), 1);
            mk{v} = sort(ix(1:k));
        end
        fprintf('Dorfler theta=%.1f marks: old %d, E-only %d, EFG %d; old==EFG set: %d\n', ...
            th, numel(mk{1}), numel(mk{2}), numel(mk{3}), isequal(mk{1}, mk{3}));
    end
end
disp('DONE_SENTINEL')
