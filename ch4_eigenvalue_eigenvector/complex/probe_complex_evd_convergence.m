function report = probe_complex_evd_convergence(outputRoot)
%PROBE_COMPLEX_EVD_CONVERGENCE - Check plain Jacobi on one full-grid endpoint
%   Usage: report = probe_complex_evd_convergence(outputRoot)
%   Purpose: Exercise the n=800 condition-sweep regime that smoke tests miss.
%       Run one mode-3 endpoint with its exact full-grid seed, recording INFO
%       rather than requiring the unpreconditioned baseline to converge.
%       This is a convergence/backward-error probe, not a forward-error test.
%   Arguments: outputRoot defaults to a fresh temporary directory.
%   Outputs: report contains parameters, INFO, counts, residual, unitarity,
%       and diagnostic-check status; also saved as plain_jacobi_large_probe.mat.
setup_complex_evd_paths();
if nargin < 1 || isempty(outputRoot), outputRoot = tempname; end
if ~isfolder(outputRoot), mkdir(outputRoot); end
n = 800; kappas = logspace(3,15,20); mode = 3; caseIndex = 20;
kappa = kappas(caseIndex); seed = 100000+1000*mode+caseIndex;
fprintf('[START] probe_complex_evd_convergence n=%d kappa=%.17g mode=%d seed=%d\n', ...
    n,kappa,mode,seed);
A = complex_randhpd(n,kappa,mode,'double',seed);
fprintf('Large probe: matrix generated; applying plain Jacobi.\n');
[V,D,rotations,sweeps,info] = cjacobi_complex(A);
residual = norm(A-V*diag(D)*V','fro')/norm(A,'fro');
unitarity = norm(V'*V-eye(n),inf);
tolerance = 100*n*(eps('double')/2);
passed = ismember(info,[0 -1]) && all(isfinite([D(:);residual;unitarity])) ...
    && residual < tolerance && unitarity < tolerance ...
    && (info == 0 || sweeps == 30);
report = struct('n',n,'kappa',kappa,'mode',mode,'caseIndex',caseIndex, ...
    'seed',seed,'info',info,'sweeps',sweeps,'rotations',rotations, ...
    'residual',residual,'unitarity',unitarity,'tolerance',tolerance, ...
    'converged',info == 0,'passed',passed,'outputRoot',outputRoot);
save(fullfile(outputRoot,'plain_jacobi_large_probe.mat'),'report');
assert(passed,'Large-probe residual, unitarity or status check failed.');
fprintf(['[END] probe_complex_evd_convergence INFO=%d sweeps=%d ' ...
    'rotations=%d residual=%.17g unitarity=%.17g\n'], ...
    info,sweeps,rotations,residual,unitarity);
end
