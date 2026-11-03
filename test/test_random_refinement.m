% refine_leaves on random and clustered markings, checked by a manufactured
% Laplace-Beltrami solve on every refined mesh (MANUFACTURED_ERROR).
%
% Each trial is RANDOM_REFINEMENT_STEPS: every patch refined once, then 10%
% of the leaves below rmax marked at each further level. A trial fails if a
% solve throws (a 2:1 balance violation) or if a mesh's error exceeds ERRFAC
% times the error on uniform level 1 (a wrong hanging-edge flag). Every
% adaptive mesh contains uniform level 1, so a correct one cannot be much
% less accurate than it.

fprintf('=== test_random_refinement start ===\n');

rmax = 4; % max quadtree depth
n = 8;    % at n = 5 a wrong split flag hides below the discretization error
nv = 3; nu = 3*nv;
errfac = 10;
ufun = @(x,y,z) sin(0.7*x).*cos(0.5*y) + 0.3*z.^2 - 0.2*x.*z;
dom0 = taylor_stellarator(n, nu, nv); % base mesh
np0 = length(dom0.x); % number of base patches

dom1 = surfacemesh.refine_leaves(dom0, [(1:np0).', zeros(np0,2)], 1:np0, rmax);
err1 = manufactured_error(dom1, ufun);
fprintf('stellarator, %d base patches; uniform level 1: %d patches, error %.3e\n', ...
    np0, length(dom1.x), err1);

nfail = 0; ntrial = 0; worst = 0;
for mname = {'random', 'clustered'}
    for trial = 1:10
        ntrial = ntrial + 1;
        steps = random_refinement_steps(dom0, mname{1}, trial, rmax);
        ratio = 0;
        for lev = 1:numel(steps)
            domR = steps(lev).dom;
            try
                ratio = max(ratio, manufactured_error(domR, ufun) / err1);
            catch ME
                fprintf('FAIL mode=%s trial=%d level=%d npat=%d : %s\n', ...
                    mname{1}, trial, lev, length(domR.x), ME.message);
                nfail = nfail + 1;
                break
            end
            if ratio > errfac
                fprintf(['FAIL mode=%s trial=%d level=%d npat=%d : error ' ...
                    '%.1f x uniform level 1\n'], mname{1}, trial, lev, ...
                    length(domR.x), ratio);
                nfail = nfail + 1;
                break
            end
        end
        worst = max(worst, ratio);
        fprintf('mode=%s trial %2d done (%d patches, max error / uniform level 1 = %.2f)\n', ...
            mname{1}, trial, length(steps(end).dom.x), ratio);
    end
end

fprintf('\nfailed trials: %d of %d; worst error ratio %.2f (limit %g)\n', ...
    nfail, ntrial, worst, errfac);
if nfail == 0
    fprintf('VERDICT: PASS\n');
else
    fprintf('VERDICT: FAIL\n');
end
fprintf('DONE_SENTINEL\n');
