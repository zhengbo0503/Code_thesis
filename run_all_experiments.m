function summary = run_all_experiments(varargin)
%RUN_ALL_EXPERIMENTS Run the thesis experiments and collect local output files.
%
%   run_all_experiments
%   run_all_experiments('DryRun', true)
%   run_all_experiments('Only', {'svals_timing', 'svals_ssd'})
%
%   Collects fresh CSV, MAT, PNG, PDF, FIG and log files in Code_thesis/output/.
%   Native experiment files are retained for compatibility with later steps.
%   No thesis directory or figure-data manifest is required or modified.
%   Files with matching names from earlier runs are replaced; unrelated older
%   files remain. summary.steps(k).outputFiles identifies this run's outputs.
%   The run log has a timestamped name; status.csv and summary.mat describe
%   the latest invocation. This function does not exit MATLAB.
%
%   Options:
%       DryRun             - Print the plan without running or writing files.
%       IncludeComplex     - Include complex EVD/SVD experiments (default true).
%       Only               - Names of experiments to run.
%       ContinueOnError    - Attempt independent later steps (default true).
%       InteractivePauses  - Allow interactive pause() checks (default false).
%       AdvanpixPath       - Directory containing the Advanpix mp class.
%       JuliaExecutable    - Path to Julia.
%       OutputDir          - Collection directory (default Code_thesis/output).
%
%   Each script has its own function workspace. Existing solver checks and
%   completion markers are retained. Numeric CSV outputs must be nonempty,
%   finite and consistent with their own header; row counts and column names
%   are not prescribed by a separate schema. Unchanged input/stale files are
%   not collected. Different outputs with the same basename are rejected.
%   Failed steps do not publish scientific results to the collection folder.
%
%   Individual experiments can still write their native data/output folders.
%   Optional MAT saves remain optional: this runner does not change savedata
%   flags or export figures that an experiment only displays.
%   A completed step is not proof of numerical correctness; inspect its
%   numerical diagnostics, including retained plain-Jacobi sweep-limit flags.

root = fileparts(mfilename('fullpath'));
parser = inputParser;
addParameter(parser, 'DryRun', false, @(x) islogical(x) && isscalar(x));
addParameter(parser, 'IncludeComplex', true, @(x) islogical(x) && isscalar(x));
addParameter(parser, 'Only', {}, @(x) iscellstr(x) || isstring(x) || ischar(x));
addParameter(parser, 'ContinueOnError', true, @(x) islogical(x) && isscalar(x));
addParameter(parser, 'InteractivePauses', false, @(x) islogical(x) && isscalar(x));
addParameter(parser, 'AdvanpixPath', '/Users/cyae/Dropbox/MATLAB/AdvanpixMCT_Mac');
addParameter(parser, 'JuliaExecutable', fullfile(getenv('HOME'), '.juliaup', 'bin', 'julia'));
addParameter(parser, 'OutputDir', fullfile(root, 'output'));
parse(parser, varargin{:});
options = parser.Results;
for optionName = {'AdvanpixPath','JuliaExecutable','OutputDir'}
    field = optionName{1};
    value = char(options.(field));
    assert(~isempty(value), '%s must be a nonempty path.', field);
    if ~startsWith(value, filesep), value = fullfile(root, value); end
    options.(field) = value;
end

jobs = experiment_plan(options.IncludeComplex);
if ~isempty(options.Only)
    requested = cellstr(options.Only);
    unknown = setdiff(requested, {jobs.name});
    assert(isempty(unknown), 'Unknown experiment(s): %s', strjoin(unknown, ', '));
    jobs = jobs(ismember({jobs.name}, requested));
end
assert(~isempty(jobs), 'No experiments selected.');
for k = 1:numel(jobs)
    assert(isfile(fullfile(root, jobs(k).relative)), ...
        'Experiment file is missing: %s', jobs(k).relative);
end

fprintf('Experiment plan (%d steps):\n', numel(jobs));
for k = 1:numel(jobs)
    fprintf('  %2d. %-28s %s\n', k, jobs(k).name, jobs(k).relative);
end
fprintf('Local output folder: %s\n', options.OutputDir);
if options.DryRun
    summary = struct('dryRun', true, 'plan', jobs, 'options', options);
    fprintf('Dry run only: no experiments run and no files written.\n');
    return
end

originalDir = pwd;
originalPath = path;
originalFigsEnv = getenv('THESIS_FIGS_DIR');
restoreState = onCleanup(@() restore_environment(originalDir, originalPath, originalFigsEnv));
if any(~strcmp({jobs.kind}, 'julia') & ~strcmp({jobs.name}, 'evals_timing'))
    assert(isfolder(options.AdvanpixPath), 'Advanpix directory is missing: %s', options.AdvanpixPath);
    addpath(options.AdvanpixPath);
end
if ~options.InteractivePauses
    addpath(fullfile(root, 'tools', 'batch_pause_guard'), '-begin');
end
runnerPath = path;

runDir = options.OutputDir;
if ~isfolder(runDir), mkdir(runDir); end
stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss_SSS'));
logFile = fullfile(runDir, ['run_' stamp '.log']);
assert(~isfile(logFile), 'Run log already exists: %s', logFile);
workDir = tempname;
mkdir(workDir);
removeWork = onCleanup(@() remove_work_directory(workDir));
owners = containers.Map('KeyType', 'char', 'ValueType', 'char');
for reserved = {'status.csv', 'summary.mat', ['run_' stamp '.log']}
    owners(lower(reserved{1})) = 'runner';
end
emptyStep = struct('name', '', 'state', 'pending', 'elapsedSeconds', 0, ...
    'message', '', 'detailLog', '', 'csvFiles', {{}}, 'outputFiles', {{}});
summary = struct('runDirectory', runDir, 'startedAt', timestamp_now(), ...
    'finishedAt', '', 'allCompleted', false, 'options', options, ...
    'environment', environment_snapshot(), ...
    'steps', repmat(emptyStep, 1, numel(jobs)));
for k = 1:numel(jobs), summary.steps(k).name = jobs(k).name; end
checkpoint(runDir, summary);
diary(logFile);
fprintf('\nRun started: %s\nLog: %s\n', summary.startedAt, logFile);
fprintf('MATLAB: %s\nAnymatrix path: %s\n', ...
    summary.environment.matlab, summary.environment.anymatrixPath);
fprintf('Anymatrix Add-On version: %s\n', summary.environment.anymatrixAddonVersion);
fprintf('Completion requires a return, checked outputs, and no recorded failure diagnostics.\n');
if ~options.InteractivePauses
    fprintf('Interactive pause() checks will fail that step; independent steps can continue.\n');
end

for k = 1:numel(jobs)
    job = jobs(k);
    dependency = find(strcmp({jobs.name}, job.dependency), 1);
    if ~isempty(dependency) && ~strcmp(summary.steps(dependency).state, 'completed')
        summary.steps(k).state = 'skipped';
        summary.steps(k).message = ['Dependency did not complete: ' job.dependency];
        fprintf('[SKIP] %s: %s\n', job.name, summary.steps(k).message);
        checkpoint(runDir, summary);
        continue
    end

    stageDir = fullfile(workDir, job.name);
    mkdir(stageDir);
    % Capture the unchanged exporters' secondary CSV copies here, rather
    % than letting their legacy default write into a thesis directory.
    setenv('THESIS_FIGS_DIR', stageDir);
    path(runnerPath);
    experimentFile = fullfile(root, job.relative);
    cd(fileparts(experimentFile));
    folders = output_folders(job, experimentFile);
    before = snapshot_outputs(folders);
    diary('off');
    offset = file_size(logFile);
    diary(logFile);
    summary.steps(k).state = 'running';
    checkpoint(runDir, summary);
    fprintf('\n[RUN %d/%d] %s at %s\n', k, numel(jobs), job.name, timestamp_now());
    started = now;
    timer = tic;
    failure = [];
    try
        switch job.kind
            case 'matlab'
                run_script_isolated(experimentFile);
            case 'julia'
                detailLog = fullfile(runDir, ['julia_timings_' stamp '.log']);
                summary.steps(k).detailLog = detailLog;
                run_julia(options.JuliaExecutable, fileparts(experimentFile), detailLog);
                validate_julia_outputs(fileparts(experimentFile), started);
            case 'complex'
                result = run_complex_isolated(job.name);
                assert(isstruct(result) && isfield(result, 'passed') && result.passed, ...
                    'Complex driver did not report completion.');
                summary.steps(k).detailLog = result.logFile;
        end
    catch exception
        failure = exception;
    end
    diary('off');
    stageText = read_log_since(logFile, offset);
    diary(logFile);
    nativeFiles = changed_outputs(folders, before);

    try
        if ~isempty(failure), rethrow(failure); end
        if ~isempty(regexp(stageText, ...
                'Error exceed predefined bound|Error:[^\r\n]*does not converge', 'once'))
            error('run_all_experiments:NumericalDiagnostic', ...
                'The script printed a bound/convergence failure; inspect run.log.');
        end
        if strcmp(job.kind, 'matlab')
            assert(contains(stageText, ['[END] ' job.name]), ...
                'The expected script END marker is missing.');
        end
        files = collect_outputs(stageDir, nativeFiles, runDir, owners, job.name);
        summary.steps(k).outputFiles = files;
        summary.steps(k).csvFiles = files(endsWith(lower(string(files)), '.csv'));
        if strcmp(job.kind, 'complex') && ~isempty(summary.steps(k).detailLog)
            [~, base, ext] = fileparts(summary.steps(k).detailLog);
            summary.steps(k).detailLog = fullfile(runDir, [base ext]);
        end
        summary.steps(k).state = 'completed';
        summary.steps(k).message = 'Execution and output checks completed.';
    catch exception
        summary.steps(k).state = 'failed';
        summary.steps(k).message = sprintf('%s: %s', exception.identifier, exception.message);
        % Keep fresh nested logs even when scientific output is incomplete.
        logs = nativeFiles(endsWith(lower(string(nativeFiles)), '.log'));
        for j = 1:numel(logs)
            [~, base, ext] = fileparts(logs{j});
            name = ['failed_' job.name '_' base ext];
            try
                publish_file(logs{j}, fullfile(runDir, name));
                summary.steps(k).detailLog = fullfile(runDir, name);
            catch logError
                fprintf(2, 'Could not collect failure log: %s\n', logError.message);
            end
        end
        fprintf(2, '[FAILURE] %s\n%s\n', job.name, getReport(exception, 'extended', 'hyperlinks', 'off'));
    end
    summary.steps(k).elapsedSeconds = toc(timer);
    fprintf('[%s] %s (%.1f seconds)\n', upper(summary.steps(k).state), ...
        job.name, summary.steps(k).elapsedSeconds);
    if ~isempty(summary.steps(k).detailLog)
        fprintf('Additional log: %s\n', summary.steps(k).detailLog);
    end
    checkpoint(runDir, summary);
    if strcmp(summary.steps(k).state, 'failed') && ~options.ContinueOnError
        for j = k+1:numel(jobs)
            summary.steps(j).state = 'skipped';
            summary.steps(j).message = 'Stopped after an earlier failure.';
        end
        break
    end
end

summary.finishedAt = timestamp_now();
summary.allCompleted = all(strcmp({summary.steps.state}, 'completed'));
summary.outputFiles = [summary.steps.outputFiles];
checkpoint(runDir, summary);
fprintf('\n================ RUN SUMMARY ================\n');
for k = 1:numel(summary.steps)
    fprintf('%-11s %-28s %9.1f s\n', summary.steps(k).state, ...
        summary.steps(k).name, summary.steps(k).elapsedSeconds);
end
fprintf('Finished: %s\nSummary: %s\nLog: %s\n', summary.finishedAt, ...
    fullfile(runDir, 'status.csv'), logFile);
diary('off');
if ~summary.allCompleted
    error('run_all_experiments:IncompleteRun', ...
        'Some steps failed or were skipped. Results are saved in %s.', runDir);
end
clear restoreState
end

function jobs = experiment_plan(includeComplex)
rows = {
    'test_off_quantity', 'ch3_preconditioning/experiments/test_off_quantity.m', 'matlab', '';
    'test_scond_reduction', 'ch3_preconditioning/experiments/test_scond_reduction.m', 'matlab', '';
    'evals_varying_cond', 'ch4_eigenvalue_eigenvector/experiments/evals_varying_cond.m', 'matlab', '';
    'evals_varying_dim', 'ch4_eigenvalue_eigenvector/experiments/evals_varying_dim.m', 'matlab', '';
    'evals_special_matrix', 'ch4_eigenvalue_eigenvector/experiments/evals_special_matrix.m', 'matlab', '';
    'julia_timings', 'ch4_eigenvalue_eigenvector/julia/test.jl', 'julia', '';
    'evals_timing', 'ch4_eigenvalue_eigenvector/experiments/evals_timing.m', 'matlab', 'julia_timings';
    'evecs_varying_cond', 'ch4_eigenvalue_eigenvector/experiments/evecs_varying_cond.m', 'matlab', '';
    'evecs_varying_dim', 'ch4_eigenvalue_eigenvector/experiments/evecs_varying_dim.m', 'matlab', '';
    'evecs_special_matrix', 'ch4_eigenvalue_eigenvector/experiments/evecs_special_matrix.m', 'matlab', '';
    'svals_varying_cond', 'ch5_singular_value_singular_vector/experiments/svals_varying_cond.m', 'matlab', '';
    'svals_varying_cols', 'ch5_singular_value_singular_vector/experiments/svals_varying_cols.m', 'matlab', '';
    'svals_special_matrix', 'ch5_singular_value_singular_vector/experiments/svals_special_matrix.m', 'matlab', '';
    'svals_timing', 'ch5_singular_value_singular_vector/experiments/svals_timing.m', 'matlab', '';
    'svals_ssd', 'ch5_singular_value_singular_vector/experiments/svals_ssd.m', 'matlab', '';
    'svecs_varying_cond', 'ch5_singular_value_singular_vector/experiments/svecs_varying_cond.m', 'matlab', '';
    'svecs_varying_cols', 'ch5_singular_value_singular_vector/experiments/svecs_varying_cols.m', 'matlab', '';
    'svecs_special_matrix', 'ch5_singular_value_singular_vector/experiments/svecs_special_matrix.m', 'matlab', '';
    'svecs_left_varying_cond', 'ch5_singular_value_singular_vector/experiments/svecs_left_varying_cond.m', 'matlab', '';
};
if includeComplex
    rows(end+1,:) = {'test_complex_evals_varying_kappa', ...
        'ch4_eigenvalue_eigenvector/complex/test_complex_evals_varying_kappa.m', 'complex', ''};
    rows(end+1,:) = {'test_complex_evals_varying_n', ...
        'ch4_eigenvalue_eigenvector/complex/test_complex_evals_varying_n.m', 'complex', ''};
    rows(end+1,:) = {'run_all_complex_full', ...
        'ch5_singular_value_singular_vector/complex/run_all_complex_full.m', 'complex', ''};
end
jobs = cell2struct(rows, {'name','relative','kind','dependency'}, 2);
jobs = jobs.';
end

function run_script_isolated(scriptFile)
% The experiment may clear this workspace, but cannot clear its caller.
run(scriptFile);
end

function result = run_complex_isolated(name)
result = feval(name);
end

function run_julia(executable, folder, logFile)
assert(isfile(executable), 'Julia executable is missing: %s', executable);
if ~isfolder(fullfile(folder, 'result')), mkdir(fullfile(folder, 'result')); end
code = ['using InteractiveUtils, LinearAlgebra, Quadmath; versioninfo(); ' ...
    'println("Quadmath version: ", Base.pkgversion(Quadmath)); ' ...
    'println("BLAS configuration: ", BLAS.get_config()); ' ...
    'println("BLAS threads: ", BLAS.get_num_threads()); ' ...
    'include("test.jl"); println("[RUNNER_END] julia_timings")'];
inner = [shell_quote(executable) ' --startup-file=no --project=' ...
    shell_quote(fullfile(folder, 'JacobiEigen')) ' -e ' shell_quote(code) ...
    ' 2>&1 | tee ' shell_quote(logFile)];
command = ['/bin/bash -o pipefail -c ' shell_quote(inner)];
[status, output] = system(command, '-echo');
assert(status == 0, 'Julia process failed with exit status %d; see %s.', status, logFile);
assert(contains(output, '[RUNNER_END] julia_timings'), 'Julia completion marker is missing.');
end

function validate_julia_outputs(folder, started)
names = {'timing_order.csv','timing_cnd.csv'};
for k = 1:numel(names)
    file = fullfile(folder, 'result', names{k});
    info = dir(file);
    assert(~isempty(info) && info.datenum >= started - 2/86400, ...
        'Julia did not refresh %s.', file);
    validate_numeric_csv(file);
end
end

function folders = output_folders(job, experimentFile)
if strcmp(job.kind, 'julia')
    folders = {fullfile(fileparts(experimentFile), 'result')};
else
    base = fileparts(experimentFile);
    if strcmp(job.kind, 'matlab'), base = fileparts(base); end
    folders = {fullfile(base, 'data'), fullfile(base, 'output'), fullfile(base, 'plots')};
end
end

function snapshot = snapshot_outputs(folders)
snapshot = containers.Map('KeyType', 'char', 'ValueType', 'any');
extensions = {'.csv','.mat','.png','.pdf','.fig','.log'};
for k = 1:numel(folders)
    if ~isfolder(folders{k}), continue; end
    entries = dir(fullfile(folders{k}, '**', '*'));
    for j = 1:numel(entries)
        entry = entries(j);
        if entry.isdir, continue; end
        [~,~,ext] = fileparts(entry.name);
        if ~ismember(lower(ext), extensions), continue; end
        file = fullfile(entry.folder, entry.name);
        snapshot(file) = [entry.bytes, entry.datenum];
    end
end
end

function files = changed_outputs(folders, before)
after = snapshot_outputs(folders);
files = {};
names = sort(keys(after));
for k = 1:numel(names)
    name = names{k};
    if ~isKey(before, name) || ~isequal(before(name), after(name))
        files{end+1} = name; %#ok<AGROW>
    end
end
end

function names = collect_outputs(stageDir, nativeFiles, destination, owners, jobName)
% Staging contains only fresh exporter CSVs, even if a file was rewritten
% with identical bytes during a filesystem timestamp's resolution window.
staged = dir(fullfile(stageDir, '*.csv'));
sources = [arrayfun(@(f) fullfile(f.folder,f.name), staged, ...
    'UniformOutput', false).', nativeFiles];
assert(~isempty(sources), 'run_all_experiments:MissingOutputs', ...
    'The experiment produced no new saved output files.');
names = {};
uniqueSources = {};
seen = containers.Map('KeyType', 'char', 'ValueType', 'double');
for k = 1:numel(sources)
    source = sources{k};
    [~,base,ext] = fileparts(source);
    name = [base ext];
    key = lower(name);
    if isKey(seen, key)
        previous = uniqueSources{seen(key)};
        assert(files_equal(previous, source), ...
            'run_all_experiments:OutputNameCollision', ...
            'Different output files share the name %s: %s and %s.', ...
            name, previous, source);
        continue
    end
    assert(~isKey(owners, key), 'run_all_experiments:OutputNameCollision', ...
        'Output name %s was already collected by another step in this run.', name);
    info = dir(source);
    assert(~isempty(info) && info.bytes > 0, 'Empty output file: %s.', source);
    if strcmpi(ext, '.csv')
        validate_numeric_csv(source);
    elseif strcmpi(ext, '.mat')
        contents = whos('-file', source);
        assert(~isempty(contents), 'No variables in MAT output: %s.', source);
    end
    names{end+1} = name; %#ok<AGROW>
    uniqueSources{end+1} = source; %#ok<AGROW>
    seen(key) = numel(names);
end
% Validate all files and names before replacing any collected scientific file.
for k = 1:numel(names)
    publish_file(uniqueSources{k}, fullfile(destination, names{k}));
    owners(lower(names{k})) = jobName;
    fprintf('[COLLECTED] %s\n', names{k});
end
end

function validate_numeric_csv(file)
fid = fopen(file, 'r');
assert(fid >= 0, 'Cannot read %s.', file);
closer = onCleanup(@() fclose(fid));
header = fgetl(fid);
assert(ischar(header) && ~isempty(strtrim(header)), 'Empty CSV: %s.', file);
columns = strsplit(strtrim(header), ',');
clear closer
values = readmatrix(file, 'NumHeaderLines', 1);
assert(~isempty(values) && size(values,2) == numel(columns), ...
    'CSV data dimensions do not match its header: %s.', file);
assert(all(isfinite(values), 'all'), 'CSV contains non-finite results: %s.', file);
end

function equal = files_equal(first, second)
a = dir(first); b = dir(second);
equal = false;
if isempty(a) || isempty(b) || a.bytes ~= b.bytes, return; end
fa = fopen(first, 'rb'); fb = fopen(second, 'rb');
if fa < 0 || fb < 0
    if fa >= 0, fclose(fa); end
    if fb >= 0, fclose(fb); end
    return
end
closeA = onCleanup(@() fclose(fa));
closeB = onCleanup(@() fclose(fb));
while ~feof(fa)
    if ~isequal(fread(fa,1048576,'*uint8'),fread(fb,1048576,'*uint8')), return; end
end
equal = true;
end

function publish_file(source, destination)
folder = fileparts(destination);
tmp = tempname(folder);
cleanup = onCleanup(@() delete_if_present(tmp));
[ok, message] = copyfile(source, tmp);
assert(ok, 'Could not stage output file: %s', message);
[ok, message] = movefile(tmp, destination, 'f');
assert(ok, 'Could not collect output file: %s', message);
clear cleanup
end

function remove_work_directory(folder)
if isfolder(folder), rmdir(folder, 's'); end
end

function checkpoint(folder, summary)
tmp = [tempname(folder) '.mat'];
save(tmp, 'summary');
[ok, message] = movefile(tmp, fullfile(folder, 'summary.mat'), 'f');
assert(ok, 'Could not save run checkpoint: %s', message);
steps = summary.steps;
status = table({steps.name}.', {steps.state}.', [steps.elapsedSeconds].', ...
    {steps.message}.', {steps.detailLog}.', ...
    'VariableNames', {'Experiment','State','ElapsedSeconds','Message','DetailLog'});
tmp = [tempname(folder) '.csv'];
writetable(status, tmp);
[ok, message] = movefile(tmp, fullfile(folder, 'status.csv'), 'f');
assert(ok, 'Could not save status CSV: %s', message);
end

function info = environment_snapshot()
info = struct('matlab', version, 'release', version('-release'), ...
    'computer', computer, 'advanpixPath', which('mp'), ...
    'anymatrixPath', which('anymatrix'), 'anymatrixAddonVersion', 'unknown');
if ~isempty(info.anymatrixPath)
    metadata = fullfile(fileparts(info.anymatrixPath), 'resources', 'addons_core.xml');
    if isfile(metadata)
        match = regexp(fileread(metadata), '<version>([^<]+)</version>', 'tokens', 'once');
        if ~isempty(match), info.anymatrixAddonVersion = match{1}; end
    end
end
end

function text = read_log_since(file, offset)
fid = fopen(file, 'rb');
if fid < 0, text = ''; return; end
cleanup = onCleanup(@() fclose(fid));
fseek(fid, offset, 'bof');
text = char(fread(fid, Inf, '*uint8').');
end

function n = file_size(file)
info = dir(file);
if isempty(info), n = 0; else, n = info.bytes; end
end

function text = timestamp_now()
text = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
end

function quoted = shell_quote(text)
single = char(39);
quoted = [single strrep(char(text), single, [single '"' single '"' single]) single];
end

function delete_if_present(file)
if isfile(file), delete(file); end
end

function restore_environment(folder, oldPath, figsEnv)
diary('off');
path(oldPath);
cd(folder);
setenv('THESIS_FIGS_DIR', figsEnv);
end
