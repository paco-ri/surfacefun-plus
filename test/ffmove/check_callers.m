% Plan item 4: the four adap_ref callers. For each, compare the new leaf set
% with the committed adap_ref_old at the test's own parameters, then run the
% test script unchanged.
ffmove = fileparts(mfilename('fullpath'));
addpath(ffmove, fullfile(ffmove, 'ref'));
sfp = fullfile(ffmove, '..');
fmt = fullfile(ffmove, '..', '..', '..', 'fmm-taylor', 'test');
addpath(fullfile(fmt, '..', 'examples'));

cases = {
    'test_adap_sphere',        sfp, @() surfacemesh.sphere(8, 1),          1e-8,  2, 1
    'test_refine_stellarator', sfp, @() surfacemesh.stellarator(8, 5, 15), 1e-8,  5, 1:10
    'testintacycadap',         fmt, @() circulartorus(7, 30, 10, 2.0, 4.0), 1e-10, 4, 1:8
    'testintbcycadap',         fmt, @() circulartorus(8, 30, 10, 2.0, 4.0), 1e-10, 4, 1:10
    };

for k = 1:size(cases, 1)
    [name, folder, mk, tol, rmax, marked] = cases{k, :};
    dom0 = mk();
    t0 = tic;
    [~, ~, pnew, hist] = surfacemesh.adap_ref(dom0, tol, rmax, 1, marked);
    tnew = toc(t0);
    fprintf('%s: new %d leaves, depth %d, %d passes (%.1f s)\n', name, ...
        size(pnew, 1), max(pnew(:, 2)), numel(hist) - 1, tnew);
    try
        evalc('[~, ~, pold] = adap_ref_old(dom0, tol, rmax, 1, marked);');
        fprintf('  old %d leaves, depth %d; only-old %d only-new %d\n', ...
            size(pold, 1), max(pold(:, 2)), ...
            size(setdiff(pold, pnew, 'rows'), 1), ...
            size(setdiff(pnew, pold, 'rows'), 1));
    catch ME
        fprintf('  old FAILS: %s\n', ME.message);
    end
    fprintf('  running %s\n', name);
    try
        run_test(folder, name);
        fprintf('  %s completed\n', name);
    catch ME
        fprintf('  %s ERROR: %s\n', name, ME.message);
    end
    close all
end
disp('DONE_SENTINEL')

function run_test(folder, name)
% A function workspace, so the test's CLEAR doesn't wipe the loop. RUN cds
% into FOLDER for the test's relative paths.
run(fullfile(folder, [name '.m']));
end
