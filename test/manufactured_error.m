function [e, perr] = manufactured_error(dom, ufun)
%MANUFACTURED_ERROR Error of a Laplace-Beltrami solve against a known solution.
%   E = MANUFACTURED_ERROR(DOM, UFUN) forms f = lap_Gamma(u), u = UFUN(x,y,z),
%   on DOM, solves lap_Gamma(uh) = f with SURFACEOP, and returns
%   max|uh - u| / max|u - mean(u)| over the nodes, after removing the
%   constant the closed-surface problem leaves free.
%
%   [E, PERR] = MANUFACTURED_ERROR(...) also returns PERR(k), the same
%   quantity over the nodes of patch k only, so E = max(PERR).
%
%   f comes from local per-patch differentiation, which does not use the
%   hanging-edge (split) flags; the SURFACEOP solve is global and does. A
%   wrong split flag therefore corrupts uh but not f, and shows up in E. A
%   2:1 balance violation makes the solve throw.

ue = surfacefun(ufun, dom);
pdo = []; pdo.lap = 1;
L = surfaceop(dom, pdo, lap(ue));
L.rankdef = true;
L.build();
uh = solve(L);
% surfacefun overloads subsref, so index one patch at a time
npat = length(dom.x);
dv = cell(npat, 1); ev = cell(npat, 1);
for k = 1:npat
    dv{k} = uh.vals{k}(:) - ue.vals{k}(:);
    ev{k} = ue.vals{k}(:);
end
mu = mean(vertcat(dv{:}));
sc = max(abs(vertcat(ev{:}) - mean(vertcat(ev{:}))));
perr = cellfun(@(v) max(abs(v - mu)), dv) / sc;
e = max(perr);
end
