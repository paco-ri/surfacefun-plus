fprintf('=== regress start\n');
tdir = '/home/paco/fmm-taylor/test/sigma_amr';
for name = {'test_refine_leaves', 'test_prolong', 'test_random_refinement'}
    t0 = tic;
    try
        out = evalc(['run(fullfile(tdir, ''' name{1} '.m''))']);
        lines = strsplit(out, newline);
        keep = lines(~cellfun(@isempty, regexp(lines, 'PASS|FAIL|VERDICT|rror|failed|max|err', 'once')));
        keep = keep(cellfun(@isempty, strfind(keep, 'DONE_SENTINEL')));
        fprintf('--- %s (%.0f s)\n', name{1}, toc(t0));
        fprintf('   %s\n', keep{:});
    catch ME
        fprintf('--- %s ERROR: %s\n', name{1}, ME.message);
    end
end
fprintf('DONE_SENTINEL\n');
