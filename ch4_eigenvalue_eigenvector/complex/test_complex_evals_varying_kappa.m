function results = test_complex_evals_varying_kappa(profile,outputRoot)
%TEST_COMPLEX_EVALS_VARYING_KAPPA - Complex two-sided eigenvalue accuracy sweep
%
%   Usage:
%       results = test_complex_evals_varying_kappa()
%       results = test_complex_evals_varying_kappa('smoke',outputRoot)
%   Purpose:
%       Adapt test_complex_varying_kappa from the one-sided complex folder.
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
fprintf('[START] test_complex_evals_varying_kappa (%s)\n',profile);
n = 800; kappas = logspace(3,15,20); modes = 1:5;

if strcmp(profile,'smoke'), n = 8; kappas = [1e3 1e8]; end
methods = {'MP3Jacobi','Jacobi','MP2Jacobi','MATLAB eig'};
nm = numel(methods); nk = numel(kappas); nmode = numel(modes);
unitRoundoff = eps('double')/2;
algorithmMpDigits = 34;
referenceMpDigits = 71;
seeds = 100000+1000*reshape(modes,1,[])+reshape(1:nk,[],1);
actualKappas = nan(nk,nmode);
prescribedKappas = nan(nk,nmode);
scalecondAt = nan(nk,nmode);
boundAt = nan(nk,nmode);
completed = false(nk,nmode);
forward = nan(nk,nmode,nm);
residual = nan(nk,nmode,nm);
unitarityV = nan(nk,nmode,nm);
info = nan(nk,nmode,3);
sweeps = nan(nk,nmode,3);
rotations = nan(nk,nmode,3);
[dataDir,plotDir] = output_dirs(root);

for jm = 1:nmode
    mode = modes(jm);
    for ik = 1:nk
        seed = seeds(ik,jm);
        fprintf('Starting MODE %d, case %d/%d: n=%d, kappa=%.3e\n', ...
            mode,ik,nk,n,kappas(ik));
        [A,~,matrixMetadata] = complex_randhpd( ...
            n,kappas(ik),mode,'double',seed);
        prescribedKappas(ik,jm) = matrixMetadata.prescribed_kappa;
        dref = reference_eigenvalues_complex(A);
        actualKappas(ik,jm) = double(max(dref)/min(dref));

        [V,D,nrot,nos,~,scalecondAt(ik,jm),info3] = mp_pjacobi_complex(A,"mp3");
        [forward(ik,jm,1),residual(ik,jm,1),unitarityV(ik,jm,1)] = ...
            compute_error_evals_complex(A,V,D,dref);
        info(ik,jm,1) = info3;
        sweeps(ik,jm,1) = nos; rotations(ik,jm,1) = nrot;
        % The real two-sided leading term is an empirical complex reference.
        boundAt(ik,jm) = 7*n*unitRoundoff*scalecondAt(ik,jm);

        [V,D,nrot,nos,infoj] = cjacobi_complex(A);
        [forward(ik,jm,2),residual(ik,jm,2),unitarityV(ik,jm,2)] = ...
            compute_error_evals_complex(A,V,D,dref);
        info(ik,jm,2) = infoj;
        sweeps(ik,jm,2) = nos; rotations(ik,jm,2) = nrot;

        [V,D,nrot,nos,~,~,info2] = mp_pjacobi_complex(A,"mp2");
        [forward(ik,jm,3),residual(ik,jm,3),unitarityV(ik,jm,3)] = ...
            compute_error_evals_complex(A,V,D,dref);
        info(ik,jm,3) = info2;
        sweeps(ik,jm,3) = nos; rotations(ik,jm,3) = nrot;

        [V,D] = eig(A,'vector');
        [forward(ik,jm,4),residual(ik,jm,4),unitarityV(ik,jm,4)] = ...
            compute_error_evals_complex(A,V,D,dref);

        assert(info3 == 0 && info2 == 0, ...
            'A preconditioned two-sided Jacobi routine did not converge.');
        if infoj ~= 0
            warning('complex_evd:PlainJacobiNotConverged', ...
                ['Plain Jacobi reached its sweep limit: n=%d, kappa=%.17g, ' ...
                 'mode=%d, seed=%d, INFO=%d. Retaining the final iterate.'], ...
                n,kappas(ik),mode,seed,infoj);
        end
        completed(ik,jm) = true;
        fprintf('varying kappa: MODE %d, %d/%d\n',mode,ik,nk);
    end
end

results = struct('profile',profile,'n',n,'kappas',kappas, ...
    'modes',modes,'methods',{methods},'forward',forward,'residual',residual, ...
    'unitarityV',unitarityV,'info',info, ...
    'sweeps',sweeps,'rotations',rotations,'seeds',seeds,'actualKappas',actualKappas, ...
    'prescribedKappas',prescribedKappas,'completed',completed,'unitRoundoff',unitRoundoff, ...
    'algorithmMpDigits',algorithmMpDigits,'referenceMpDigits',referenceMpDigits, ...
    'scalecondAt',scalecondAt,'boundAt',boundAt, ...
    'referenceMeaning','Real-case leading term; empirical in complex arithmetic', ...
    'passed',all(completed(:)),'allJacobiConverged',all(info(:) == 0),'logFile','','matlab_release',version('-release'));

save(fullfile(dataDir,['complex_evals_varying_kappa_' profile '.mat']),'results','-v7.3');
write_complex_evd_csv(results,'varying_kappa',dataDir);
plot_complex_evd_results(results,'varying_kappa',plotDir);
fprintf('[END] test_complex_evals_varying_kappa (%s)\n',profile);
end

function [dataDir,plotDir] = output_dirs(root)
%OUTPUT_DIRS - Create and return the experiment data and plot directories.
dataDir = fullfile(root,'data'); plotDir = fullfile(root,'plots');
if ~isfolder(dataDir), mkdir(dataDir); end
if ~isfolder(plotDir), mkdir(plotDir); end
end
