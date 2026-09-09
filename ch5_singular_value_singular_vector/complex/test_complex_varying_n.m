function results = test_complex_varying_n()
%TEST_COMPLEX_VARYING_N - Test complex SVD accuracy as n varies
%
%   Usage:
%       results = test_complex_varying_n()
%
%   Purpose:
%       TEST_COMPLEX_VARYING_N runs the complex counterpart of Figure 2 for
%       MP3JacobiSVD, ZGESVJ, ZGEJSV, and MATLAB svd. It measures forward
%       error, residual, and loss of unitarity using 71-digit reference
%       singular values, evaluates the practical scond(At) reference, and
%       records the QR/direct branch selected at each matrix size. The
%       function saves the result structure and exports PNG/PDF plots.
%
%   Arguments:
%       None.
%
%   Outputs:
%       (1) results - Structure containing the profile parameters, method
%           names, seeds, actual condition numbers, accuracy metrics,
%           practical references, branch indicators, LAPACK INFO values,
%           sweep counts, completion flags, precision settings, and MATLAB
%           release.
%

profile = 'full';
root = setup_complex_paths();
m = 1000; nvalues = unique(round(logspace(1,3,15))); modes = 1:5;
kappa = 1e8;
methods = {'MP3JacobiSVD','ZGESVJ','ZGEJSV','MATLAB svd'};
nm = numel(methods); nn = numel(nvalues); nmode = numel(modes);
unitRoundoff = eps('double')/2;
algorithmMpDigits = 34;
referenceMpDigits = 71;
seeds = 200000+1000*reshape(modes,1,[])+reshape(1:nn,[],1);
actualKappas = nan(nn,nmode);
scalecondAt = nan(nn,nmode);
boundAt = nan(nn,nmode);
completed = false(nn,nmode);
forward = nan(nn,nmode,nm);
residual = nan(nn,nmode,nm);
unitarityU = nan(nn,nmode,nm);
unitarityV = nan(nn,nmode,nm);
info = nan(nn,nmode,2);
sweeps = nan(nn,nmode,2);
doqr = m >= (11*nvalues)/6;
[dataDir,plotDir] = output_dirs(root);

for jm = 1:nmode
    mode = modes(jm);
    for in = 1:nn
        n = nvalues(in);
        seed = seeds(in,jm);
        [A,~,matrixMetadata] = complex_randsvd( ...
            m,n,kappa,mode,'double',seed);
        actualKappas(in,jm) = matrixMetadata.actual_kappa;
        sref = reference_singular_values(A);

        [U,S,V,nos,scalecondAt(in,jm)] = mposj_complex(A,3,true);
        [forward(in,jm,1),residual(in,jm,1),unitarityU(in,jm,1),unitarityV(in,jm,1)] = ...
            compute_error_svals(A,U,S,V,sref);
        sweeps(in,jm,1) = nos;
        boundAt(in,jm) = scalecondAt(in,jm)*sqrt(m*n)*unitRoundoff;

        [U,S,V,~,rw,infoz] = zgesvj_mex(A,'G');
        [forward(in,jm,2),residual(in,jm,2),unitarityU(in,jm,2),unitarityV(in,jm,2)] = ...
            compute_error_svals(A,U,S,V,sref);
        info(in,jm,1) = infoz;
        sweeps(in,jm,2) = double(rw(4));

        [U,S,V,~,~,~,infoe] = zgejsv_mex(A);
        [forward(in,jm,3),residual(in,jm,3),unitarityU(in,jm,3),unitarityV(in,jm,3)] = ...
            compute_error_svals(A,U,S,V,sref);
        info(in,jm,2) = infoe;

        [U,S,V] = svd(A,'econ');
        [forward(in,jm,4),residual(in,jm,4),unitarityU(in,jm,4),unitarityV(in,jm,4)] = ...
            compute_error_svals(A,U,S,V,sref);

        assert(infoz == 0 && infoe == 0,'A LAPACK routine did not converge.');
        completed(in,jm) = true;
        fprintf('varying n: MODE %d, %d/%d\n',mode,in,nn);
    end
end

results = struct('profile',profile,'m',m,'nvalues',nvalues,'kappa',kappa, ...
    'modes',modes,'methods',{methods},'forward',forward,'residual',residual, ...
    'unitarityU',unitarityU,'unitarityV',unitarityV,'info',info, ...
    'sweeps',sweeps,'doqr',doqr,'branch_threshold',6*m/11, ...
    'seeds',seeds,'actualKappas',actualKappas, ...
    'completed',completed,'unitRoundoff',unitRoundoff, ...
    'algorithmMpDigits',algorithmMpDigits,'referenceMpDigits',referenceMpDigits, ...
    'scalecondAt',scalecondAt,'boundAt',boundAt, ...
    'matlab_release',version('-release'));

save(fullfile(dataDir,['complex_varying_n_' profile '.mat']),'results','-v7.3');
plot_complex_results(results,'varying_n',plotDir);
end

function [dataDir,plotDir] = output_dirs(root)
%OUTPUT_DIRS - Create and return the experiment data and plot directories.
dataDir = fullfile(root,'data'); plotDir = fullfile(root,'plots');
if ~isfolder(dataDir), mkdir(dataDir); end
if ~isfolder(plotDir), mkdir(plotDir); end
end
