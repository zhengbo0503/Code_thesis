function results = test_complex_evals_varying_n(profile,outputRoot)
%TEST_COMPLEX_EVALS_VARYING_N - Complex two-sided eigenvalue accuracy sweep
%
%   Usage:
%       results = test_complex_evals_varying_n()
%       results = test_complex_evals_varying_n('smoke',outputRoot)
%   Purpose:
%       Adapt test_complex_varying_n from the one-sided complex folder.
%       Keep its modes, full grid, seeds, metrics and tiled plot structure;
%       use square Hermitian matrices and two-sided EVD routines.
%       Compare MP3Jacobi, Jacobi, MP2Jacobi and MATLAB eig to 71-digit
%       eigenvalues of the stored input. No runtime measurements are made.
%   Arguments:
%       profile defaults to 'full'; 'smoke' reduces the grid for validation.
%       outputRoot optionally relocates the data/ and plots/ directories.
%   Outputs:
%       results contains parameters, seeds, errors, unitarity, scaled
%       conditions, empirical references, INFO, rotations and sweep counts.
%       Saves MAT/CSV data and PNG/PDF plots after all cases complete.
%       Plain Jacobi may reach its sweep limit; retain its INFO and final
%       iterate. MP2/MP3 nonconvergence still fails the experiment.
%
if nargin < 1, profile = 'full'; end
profile = validatestring(profile,{'full','smoke'});
root = setup_complex_evd_paths();
if nargin >= 2 && ~isempty(outputRoot), root = outputRoot; end
fprintf('[START] test_complex_evals_varying_n (%s)\n',profile);
nvalues = unique(round(logspace(1,3,15))); modes = 1:5;
kappa = 1e8;
if strcmp(profile,'smoke'), nvalues = [6 10]; end
methods = {'MP3Jacobi','Jacobi','MP2Jacobi','MATLAB eig'};
nm = numel(methods); nn = numel(nvalues); nmode = numel(modes);
unitRoundoff = eps('double')/2;
algorithmMpDigits = 34;
referenceMpDigits = 71;
seeds = 200000+1000*reshape(modes,1,[])+reshape(1:nn,[],1);
actualKappas = nan(nn,nmode);
prescribedKappas = nan(nn,nmode);
scalecondAt = nan(nn,nmode);
boundAt = nan(nn,nmode);
completed = false(nn,nmode);
forward = nan(nn,nmode,nm);
residual = nan(nn,nmode,nm);
unitarityV = nan(nn,nmode,nm);
info = nan(nn,nmode,3);
sweeps = nan(nn,nmode,3);
rotations = nan(nn,nmode,3);
[dataDir,plotDir] = output_dirs(root);

for jm = 1:nmode
    mode = modes(jm);
    for in = 1:nn
        n = nvalues(in);
        seed = seeds(in,jm);
        fprintf('Starting MODE %d, case %d/%d: n=%d, kappa=%.3e\n', ...
            mode,in,nn,n,kappa);
        [A,~,matrixMetadata] = complex_randhpd( ...
            n,kappa,mode,'double',seed);
        prescribedKappas(in,jm) = matrixMetadata.prescribed_kappa;
        dref = reference_eigenvalues_complex(A);
        actualKappas(in,jm) = double(max(dref)/min(dref));

        [V,D,nrot,nos,~,scalecondAt(in,jm),info3] = mp_pjacobi_complex(A,"mp3");
        [forward(in,jm,1),residual(in,jm,1),unitarityV(in,jm,1)] = ...
            compute_error_evals_complex(A,V,D,dref);
        info(in,jm,1) = info3;
        sweeps(in,jm,1) = nos; rotations(in,jm,1) = nrot;
        % The real two-sided leading term is an empirical complex reference.
        boundAt(in,jm) = 7*n*unitRoundoff*scalecondAt(in,jm);

        [V,D,nrot,nos,infoj] = cjacobi_complex(A);
        [forward(in,jm,2),residual(in,jm,2),unitarityV(in,jm,2)] = ...
            compute_error_evals_complex(A,V,D,dref);
        info(in,jm,2) = infoj;
        sweeps(in,jm,2) = nos; rotations(in,jm,2) = nrot;

        [V,D,nrot,nos,~,~,info2] = mp_pjacobi_complex(A,"mp2");
        [forward(in,jm,3),residual(in,jm,3),unitarityV(in,jm,3)] = ...
            compute_error_evals_complex(A,V,D,dref);
        info(in,jm,3) = info2;
        sweeps(in,jm,3) = nos; rotations(in,jm,3) = nrot;

        [V,D] = eig(A,'vector');
        [forward(in,jm,4),residual(in,jm,4),unitarityV(in,jm,4)] = ...
            compute_error_evals_complex(A,V,D,dref);

        assert(info3 == 0 && info2 == 0, ...
            'A preconditioned two-sided Jacobi routine did not converge.');
        if infoj ~= 0
            warning('complex_evd:PlainJacobiNotConverged', ...
                ['Plain Jacobi reached its sweep limit: n=%d, kappa=%.17g, ' ...
                 'mode=%d, seed=%d, INFO=%d. Retaining the final iterate.'], ...
                n,kappa,mode,seed,infoj);
        end
        completed(in,jm) = true;
        fprintf('varying n: MODE %d, %d/%d\n',mode,in,nn);
    end
end

results = struct('profile',profile,'nvalues',nvalues,'kappa',kappa, ...
    'modes',modes,'methods',{methods},'forward',forward,'residual',residual, ...
    'unitarityV',unitarityV,'info',info, ...
    'sweeps',sweeps,'rotations',rotations, ...
    'seeds',seeds,'actualKappas',actualKappas, ...
    'prescribedKappas',prescribedKappas,'completed',completed,'unitRoundoff',unitRoundoff, ...
    'algorithmMpDigits',algorithmMpDigits,'referenceMpDigits',referenceMpDigits, ...
    'scalecondAt',scalecondAt,'boundAt',boundAt, ...
    'referenceMeaning','Real-case leading term; empirical in complex arithmetic', ...
    'passed',all(completed(:)),'allJacobiConverged',all(info(:) == 0),'logFile','','matlab_release',version('-release'));

save(fullfile(dataDir,['complex_evals_varying_n_' profile '.mat']),'results','-v7.3');
write_complex_evd_csv(results,'varying_n',dataDir);
plot_complex_evd_results(results,'varying_n',plotDir);
fprintf('[END] test_complex_evals_varying_n (%s)\n',profile);
end

function [dataDir,plotDir] = output_dirs(root)
%OUTPUT_DIRS - Create and return the experiment data and plot directories.
dataDir = fullfile(root,'data'); plotDir = fullfile(root,'plots');
if ~isfolder(dataDir), mkdir(dataDir); end
if ~isfolder(plotDir), mkdir(plotDir); end
end
