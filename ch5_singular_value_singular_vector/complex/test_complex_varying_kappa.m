function results = test_complex_varying_kappa()
%TEST_COMPLEX_VARYING_KAPPA - Test complex SVD accuracy as kappa varies
%
%   Usage:
%       results = test_complex_varying_kappa()
%
%   Purpose:
%       TEST_COMPLEX_VARYING_KAPPA runs the complex counterpart of Figure 1
%       for MP3JacobiSVD, ZGESVJ, ZGEJSV, and MATLAB svd. It measures forward
%       error, residual, and loss of unitarity using 71-digit reference
%       singular values, and evaluates the practical scond(At) reference.
%       The function saves the result structure and exports PNG/PDF plots.
%
%   Arguments:
%       None.
%
%   Outputs:
%       (1) results - Structure containing the profile parameters, method
%           names, seeds, actual condition numbers, accuracy metrics,
%           practical references, LAPACK INFO values, sweep counts,
%           completion flags, precision settings, and MATLAB release.
%

profile = 'full';
root = setup_complex_paths();
m = 1000; n = 800; kappas = logspace(3,15,20); modes = 1:5;

methods = {'MP3JacobiSVD','ZGESVJ','ZGEJSV','MATLAB svd'};
nm = numel(methods); nk = numel(kappas); nmode = numel(modes);
unitRoundoff = eps('double')/2;
algorithmMpDigits = 34;
referenceMpDigits = 71;
seeds = 100000+1000*reshape(modes,1,[])+reshape(1:nk,[],1);
actualKappas = nan(nk,nmode);
scalecondAt = nan(nk,nmode);
boundAt = nan(nk,nmode);
completed = false(nk,nmode);
forward = nan(nk,nmode,nm);
residual = nan(nk,nmode,nm);
unitarityU = nan(nk,nmode,nm);
unitarityV = nan(nk,nmode,nm);
info = nan(nk,nmode,2);
sweeps = nan(nk,nmode,2);
[dataDir,plotDir] = output_dirs(root);

for jm = 1:nmode
    mode = modes(jm);
    for ik = 1:nk
        seed = seeds(ik,jm);
        [A,~,matrixMetadata] = complex_randsvd( ...
            m,n,kappas(ik),mode,'double',seed);
        actualKappas(ik,jm) = matrixMetadata.actual_kappa;
        sref = reference_singular_values(A);

        [U,S,V,nos,scalecondAt(ik,jm)] = mposj_complex(A,3,true);
        [forward(ik,jm,1),residual(ik,jm,1),unitarityU(ik,jm,1),unitarityV(ik,jm,1)] = ...
            compute_error_svals(A,U,S,V,sref);
        sweeps(ik,jm,1) = nos;
        boundAt(ik,jm) = scalecondAt(ik,jm)*sqrt(m*n)*unitRoundoff;

        [U,S,V,~,rw,infoz] = zgesvj_mex(A,'G');
        [forward(ik,jm,2),residual(ik,jm,2),unitarityU(ik,jm,2),unitarityV(ik,jm,2)] = ...
            compute_error_svals(A,U,S,V,sref);
        info(ik,jm,1) = infoz;
        sweeps(ik,jm,2) = double(rw(4));

        [U,S,V,~,~,~,infoe] = zgejsv_mex(A);
        [forward(ik,jm,3),residual(ik,jm,3),unitarityU(ik,jm,3),unitarityV(ik,jm,3)] = ...
            compute_error_svals(A,U,S,V,sref);
        info(ik,jm,2) = infoe;

        [U,S,V] = svd(A,'econ');
        [forward(ik,jm,4),residual(ik,jm,4),unitarityU(ik,jm,4),unitarityV(ik,jm,4)] = ...
            compute_error_svals(A,U,S,V,sref);

        assert(infoz == 0 && infoe == 0,'A LAPACK routine did not converge.');
        completed(ik,jm) = true;
        fprintf('varying kappa: MODE %d, %d/%d\n',mode,ik,nk);
    end
end

results = struct('profile',profile,'m',m,'n',n,'kappas',kappas, ...
    'modes',modes,'methods',{methods},'forward',forward,'residual',residual, ...
    'unitarityU',unitarityU,'unitarityV',unitarityV,'info',info, ...
    'sweeps',sweeps,'seeds',seeds,'actualKappas',actualKappas, ...
    'completed',completed,'unitRoundoff',unitRoundoff, ...
    'algorithmMpDigits',algorithmMpDigits,'referenceMpDigits',referenceMpDigits, ...
    'scalecondAt',scalecondAt,'boundAt',boundAt, ...
    'matlab_release',version('-release'));

save(fullfile(dataDir,['complex_varying_kappa_' profile '.mat']),'results','-v7.3');
plot_complex_results(results,'varying_kappa',plotDir);
end

function [dataDir,plotDir] = output_dirs(root)
%OUTPUT_DIRS - Create and return the experiment data and plot directories.
dataDir = fullfile(root,'data'); plotDir = fullfile(root,'plots');
if ~isfolder(dataDir), mkdir(dataDir); end
if ~isfolder(plotDir), mkdir(plotDir); end
end
