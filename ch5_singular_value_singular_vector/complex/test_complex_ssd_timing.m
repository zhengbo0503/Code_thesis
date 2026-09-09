function results = test_complex_ssd_timing()
%TEST_COMPLEX_SSD_TIMING - Test complex SSD accuracy and runtime
%
%   Usage:
%       results = test_complex_ssd_timing()
%
%   Purpose:
%       TEST_COMPLEX_SSD_TIMING compares MP3JacobiSVD-SSD, CGESVJ, and
%       CGEJSV on complex single-precision matrices. At each matrix size it
%       discards the first execution of each method and reports the mean of
%       the second and third executions. It also evaluates Figure 5-style
%       forward errors using 71-digit reference singular values and the
%       practical scond(At) reference outside the timed runs. The function
%       saves the result structure and exports PNG/PDF plots.
%
%   Arguments:
%       None.
%
%   Outputs:
%       (1) results - Structure containing the profile parameters, methods,
%           seeds, raw and measured timings, reported mean times, accuracy
%           metrics, practical references, branch indicators, LAPACK INFO
%           values, precision settings, and complete timing protocols.
%

profile = 'full';
root = setup_complex_paths();
m = 1000;
nvalues = unique(round(logspace(2,3,10)));
discardedRuns = 1;
measuredRuns = 2;
totalRuns = discardedRuns+measuredRuns;
protocol = ['For each matrix size, run methods in the fixed order ' ...
    'MP3JacobiSVD-SSD, CGESVJ, CGEJSV three times; discard run 1 ' ...
    'for every method and report the arithmetic mean of runs 2 and 3. ' ...
    'Reference singular values and scond(At) are evaluated outside all ' ...
    'timed runs.'];
accuracyProtocol = ['Forward errors use one 71-digit reference per matrix; ' ...
    'the practical reference is sqrt(m*n)*u*scond(At), with single ' ...
    'unit roundoff and At evaluated before demotion to single.'];
kappa = 1e6;
mode = 3;
workingUnitRoundoff = eps('single')/2;
referenceMpDigits = 71;
methods = {'MP3JacobiSVD-SSD','CGESVJ','CGEJSV'};
nm = numel(methods);
nn = numel(nvalues);
rawTimes = nan(nn,totalRuns,nm);
seeds = 300000+(1:nn);
forward = nan(nn,nm);
residual = nan(nn,nm);
scalecondAt = nan(nn,1);
boundAt = nan(nn,1);
info = nan(nn,2);
doqr = m >= (11*nvalues)/6;

% One untimed call removes MEX loading and JIT setup from measured runs.
Awarm = complex_randsvd(12,6,1e3,3,'single',300001);
mposj_ssd_complex(Awarm,false);
cgesvj_mex(Awarm,'G');
cgejsv_mex(Awarm);

for in = 1:nn
    n = nvalues(in);
    A = complex_randsvd(m,n,kappa,mode,'single',seeds(in));

    % Discard run 1 for every method and average runs 2 and 3.
    for r = 1:totalRuns
        t = tic;
        [U1,S1,V1] = mposj_ssd_complex(A,false);
        rawTimes(in,r,1) = toc(t);

        t = tic;
        [U2,S2,V2,~,~,info2] = cgesvj_mex(A,'G');
        rawTimes(in,r,2) = toc(t);

        t = tic;
        [U3,S3,V3,~,~,~,info3] = cgejsv_mex(A);
        rawTimes(in,r,3) = toc(t);

    end

    info(in,:) = [info2,info3];
    assert(info2 == 0 && info3 == 0,'A LAPACK routine did not converge.');

    % The Figure 5 accuracy diagnostics are deliberately outside the
    % measured runs. The singular-value reference uses 71-digit arithmetic,
    % and scalecondAt is evaluated on the double-precision preconditioned
    % matrix before its demotion to single precision.
    sref = reference_singular_values(A);
    [forward(in,1),residual(in,1)] = compute_error_svals( ...
        double(A),double(U1),double(S1),double(V1),sref);
    [forward(in,2),residual(in,2)] = compute_error_svals( ...
        double(A),double(U2),double(S2),double(V2),sref);
    [forward(in,3),residual(in,3)] = compute_error_svals( ...
        double(A),double(U3),double(S3),double(V3),sref);
    [~,~,~,~,scalecondAt(in)] = mposj_ssd_complex(A,true);
    boundAt(in) = sqrt(m*n)*workingUnitRoundoff*scalecondAt(in);

    fprintf('SSD timing: n=%d, %d/%d\n',n,in,nn);
end

measuredTimes = rawTimes(:,discardedRuns+1:end,:);
meanTimes = reshape(mean(measuredTimes,2),nn,nm);
results = struct('profile',profile,'m',m,'nvalues',nvalues,'kappa',kappa, ...
    'mode',mode,'repeats',measuredRuns,'discardedRuns',discardedRuns, ...
    'measuredRuns',measuredRuns,'totalRuns',totalRuns,'protocol',protocol, ...
    'accuracyProtocol',accuracyProtocol, ...
    'methods',{methods}, ...
    'rawTimes',rawTimes,'measuredTimes',measuredTimes, ...
    'meanTimes',meanTimes,'forward',forward,'residual',residual, ...
    'scalecondAt',scalecondAt,'boundAt',boundAt, ...
    'workingUnitRoundoff',workingUnitRoundoff, ...
    'referenceMpDigits',referenceMpDigits, ...
    'info',info,'doqr',doqr,'branch_threshold',6*m/11, ...
    'seeds',seeds,'matlab_release',version('-release'));

dataDir = fullfile(root,'data'); plotDir = fullfile(root,'plots');
if ~isfolder(dataDir), mkdir(dataDir); end
if ~isfolder(plotDir), mkdir(plotDir); end
save(fullfile(dataDir,['complex_ssd_timing_' profile '.mat']),'results','-v7.3');
plot_complex_results(results,'ssd_timing',plotDir);
end
