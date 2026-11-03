function steps = random_refinement_steps(dom0, mname, seed, rmax)
%RANDOM_REFINEMENT_STEPS Meshes of one seeded random or clustered refinement.
%   STEPS = RANDOM_REFINEMENT_STEPS(DOM0, MNAME, SEED, RMAX) refines every
%   patch of DOM0 once, then at each further level marks 10% of the leaves
%   shallower than RMAX:
%     'random'     chosen uniformly at random, scattered over the surface
%     'clustered'  the ones nearest a randomly chosen leaf, so one region is
%                  refined repeatedly
%   STEPS(lev) holds DOM, the level-lev mesh, and P2Q, its leaf set, so
%   REFINE_LEAVES(DOM0, STEPS(lev).P2Q, [], RMAX) rebuilds it. 'random' leaf
%   sets depend only on SEED, not on the polynomial order of DOM0.

rng(seed);
np0 = length(dom0.x);
p2q = [(1:np0).', zeros(np0, 2)];
steps = struct('dom', {}, 'p2q', {});
for lev = 1:rmax
    elig = find(p2q(:, 2) < rmax);
    if isempty(elig), break, end
    k = max(1, round(0.10*numel(elig)));
    if lev == 1
        marked = (1:size(p2q, 1)).';
    elseif strcmp(mname, 'random')
        marked = elig(randperm(numel(elig), k));
    else
        dom = steps(lev-1).dom;
        cx = cellfun(@(a) mean(a(:)), dom.x);
        cy = cellfun(@(a) mean(a(:)), dom.y);
        cz = cellfun(@(a) mean(a(:)), dom.z);
        s = elig(randi(numel(elig)));
        dd = (cx(elig)-cx(s)).^2 + (cy(elig)-cy(s)).^2 + (cz(elig)-cz(s)).^2;
        [~, ord] = sort(dd);
        marked = elig(ord(1:k));
    end
    [steps(lev).dom, ~, p2q] = surfacemesh.refine_leaves(dom0, p2q, marked, rmax);
    steps(lev).p2q = p2q;
end
end
