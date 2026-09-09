function summary = run_all_complex_full()
%RUN_ALL_COMPLEX_FULL - Validate and run the complex full profiles unattended
%
%   Usage:
%       summary = run_all_complex_full()
%
%   Purpose:
%       RUN_ALL_COMPLEX_FULL prepares the MEX wrappers, validates the complex
%       stack, and runs the varying-kappa, varying-n, and SSD full profiles in
%       sequence. During a full run it uses macOS caffeinate to prevent system
%       sleep, writes a timestamped diary, and updates an atomic status file
%       after each stage.
%
%   Arguments:
%       None.
%
%   Outputs:
%       (1) summary - Structure containing the run timestamp, stage
%           status and durations, completion and error fields, output paths,
%           pass flag, and MATLAB release.
%

root = setup_complex_paths();
assert(ismember(version('-release'),{'2025b','2026a'}), ...
    'The complex experiments target MATLAB R2025b or R2026a (found R%s).', version('-release'));
dataDir = fullfile(root,'data');
if ~isfolder(dataDir), mkdir(dataDir); end

runStamp = char(datetime('now','Format','yyyyMMdd_HHmmss'));
logFile = fullfile(dataDir,['complex_full_run_' runStamp '.log']);
statusFile = fullfile(dataDir,['complex_full_run_' runStamp '_status.mat']);

diary(logFile);
stopDiary = onCleanup(@() diary('off'));
caffeinatePid = start_caffeinate();
stopCaffeinate = onCleanup(@() stop_caffeinate(caffeinatePid));

summary = struct('runStamp',runStamp,'startedAt',timestamp_now(),'finishedAt','', ...
    'currentStage','initialising','completedStages',{{}}, ...
    'stageDurationsSeconds',struct(), ...
    'passed',false,'errorIdentifier','','errorMessage','', ...
    'logFile',logFile,'statusFile',statusFile, ...
    'matlabRelease',version('-release'));
save_summary_atomic(statusFile,summary);

fprintf('Complex full-run driver started at %s.\n',summary.startedAt);
fprintf('Log: %s\nStatus: %s\n',logFile,statusFile);

try
    ensure_mex_current(root);

    summary.currentStage = 'validation';
    save_summary_atomic(statusFile,summary);
    stageTimer = tic;
    validate_complex_stack();
    summary.stageDurationsSeconds.validation = toc(stageTimer);
    summary.completedStages{end+1} = 'validation';
    save_summary_atomic(statusFile,summary);

    stageNames = {'varying_kappa','varying_n','ssd_timing'};
    stageFunctions = {@test_complex_varying_kappa, ...
        @test_complex_varying_n,@test_complex_ssd_timing};
    for k = 1:numel(stageNames)
        stage = stageNames{k};
        summary.currentStage = stage;
        save_summary_atomic(statusFile,summary);
        fprintf('Starting %s at %s.\n',stage,timestamp_now());
        stageTimer = tic;
        stageFunctions{k}();
        summary.stageDurationsSeconds.(stage) = toc(stageTimer);
        summary.completedStages{end+1} = stage;
        save_summary_atomic(statusFile,summary);
        fprintf('Completed %s at %s.\n',stage,timestamp_now());
    end

    summary.currentStage = 'complete';
    summary.finishedAt = timestamp_now();
    summary.passed = true;
    save_summary_atomic(statusFile,summary);
    fprintf('All complex full experiments completed.\n');
    fprintf('Finished at %s.\n',summary.finishedAt);
catch exception
    summary.finishedAt = timestamp_now();
    summary.passed = false;
    summary.errorIdentifier = exception.identifier;
    summary.errorMessage = exception.message;
    save_summary_atomic(statusFile,summary);
    fprintf(2,'COMPLEX_FULL_RUN_FAILED during %s: %s\n', ...
        summary.currentStage,exception.message);
    rethrow(exception)
end

clear stopCaffeinate stopDiary
end

function ensure_mex_current(root)
%ENSURE_MEX_CURRENT - Rebuild all complex MEX files when one is stale.
names = {'zgesvj','zgejsv','cgesvj','cgejsv'};
needsBuild = false;
for k = 1:numel(names)
    source = fullfile(root,['get_' names{k}],[names{k} '_mex.c']);
    binary = fullfile(root,[names{k} '_mex.' mexext]);
    sourceInfo = dir(source);
    binaryInfo = dir(binary);
    if isempty(binaryInfo) || binaryInfo.datenum < sourceInfo.datenum
        needsBuild = true;
        break
    end
end
if needsBuild
    fprintf('A complex MEX binary is missing or stale; rebuilding all four.\n');
    build_all_complex_mex();
end
end

function pid = start_caffeinate()
%START_CAFFEINATE - Start a sleep-prevention process tied to MATLAB.
matlabPid = feature('getpid');
command = sprintf( ...
    'caffeinate -ims -w %d >/dev/null 2>&1 & echo $!',matlabPid);
[status,output] = system(command);
pid = str2double(strtrim(output));
assert(status == 0 && isfinite(pid) && pid > 0, ...
    'run_all_complex_full:CaffeinateFailed', ...
    'Could not start macOS caffeinate; no experiment was started.');
fprintf('macOS sleep prevention active (caffeinate PID %d).\n',pid);
end

function stop_caffeinate(pid)
%STOP_CAFFEINATE - Stop the sleep-prevention process when it is running.
if isfinite(pid) && pid > 0
    system(sprintf('kill %d >/dev/null 2>&1',floor(pid)));
end
end

function save_summary_atomic(statusFile,summary)
%SAVE_SUMMARY_ATOMIC - Save and atomically replace the run-status file.
temporaryFile = [tempname(fileparts(statusFile)) '.mat'];
removeTemporaryFile = onCleanup(@() delete_if_present(temporaryFile));
save(temporaryFile,'summary','-v7.3');
[moved,message] = movefile(temporaryFile,statusFile,'f');
assert(moved,'run_all_complex_full:StatusMoveFailed', ...
    'Could not replace run-status file: %s',message);
clear removeTemporaryFile
end

function delete_if_present(file)
%DELETE_IF_PRESENT - Delete a file when it exists.
if isfile(file)
    delete(file);
end
end

function value = timestamp_now()
%TIMESTAMP_NOW - Return a local timestamp with its UTC offset.
value = char(datetime('now','TimeZone','local', ...
    'Format','yyyy-MM-dd HH:mm:ss Z'));
end
