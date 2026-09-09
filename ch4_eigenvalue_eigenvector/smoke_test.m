function nfail = smoke_test()
%SMOKE_TEST - Chapter 4 smoke test
%
%   Usage:
%       nfail = smoke_test()
%
%   Purpose:
%       Exercises every function in ch4_eigenvalue_eigenvector/shared on
%       small inputs and checks the defining property of each: cjacobi and
%       mp_pjacobi reconstruct the matrix, produce orthogonal eigenvectors
%       and descending eigenvalues, and beat eig on relative accuracy for
%       an ill-conditioned SPD matrix.
%
%       Checks needing Advanpix (mp) are skipped with a notice if it is not
%       on the path. Add it with
%           addpath('/Users/cyae/Dropbox/MATLAB/AdvanpixMCT_Mac')
%
%   Author:
%       Smoke test for the thesis code reorganisation.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'shared'));
R = {};
hasmp = ~isempty(which('mp'));

fprintf('=== Chapter 4: eigenvalues and eigenvectors ===\n');
fprintf('  Advanpix mp: %d\n', hasmp);
rng(0);
n = 30;
A = gallery('randsvd', n, -1e8, 3);       % SPD, kappa 1e8
A = (A + A')/2;

% --- float_params / off / jacobi_pair ----------------------------------
try
    u = float_params('d');
    R = chk(R,'float_params(d)', abs(u - 2^-53) < eps, sprintf('%.3e', u));
    R = chk(R,'off(diagonal) == 0', off(diag(rand(n,1))) < 1e-14, '');
    [c,s] = jacobi_pair(A,1,2);
    R = chk(R,'jacobi_pair is a rotation', abs(c^2 + s^2 - 1) < 1e-14, sprintf('c^2+s^2-1=%.2e', c^2+s^2-1));
    J = [c s; -s c]; B = A([1 2],[1 2]); Bt = J'*B*J;
    R = chk(R,'jacobi_pair zeroes A(p,q)', abs(Bt(1,2)) < 1e-12*norm(B,'fro'), sprintf('%.2e', abs(Bt(1,2))));
catch e
    R = bad(R,'float_params/off/jacobi_pair', e);
end

% --- eigsort ------------------------------------------------------------
try
    [Vs,Ds] = eigsort(A);
    d = diag(Ds);
    R = chk(R,'eigsort descending', all(diff(d) <= 1e-12*max(abs(d))), '');
    R = chk(R,'eigsort reconstruction', relerr(A, Vs*Ds*Vs') < 1e-12, sprintf('%.3e', relerr(A,Vs*Ds*Vs')));
catch e
    R = bad(R,'eigsort', e);
end

% --- cjacobi (the merged version: faithful Demmel-Veselic criterion) ----
try
    [V,D,NROT,NSWEEP,INFO] = cjacobi(A);
    d = asvec(D);
    R = chk(R,'cjacobi INFO == 0', INFO == 0, sprintf('INFO=%d NROT=%d NSWEEP=%d', INFO, NROT, NSWEEP));
    R = chk(R,'cjacobi converged in < 30 sweeps', NSWEEP < 30, sprintf('%d sweeps', NSWEEP));
    R = chk(R,'cjacobi V orthogonal', norm(V'*V - eye(n),'fro') < 1e-12, sprintf('%.3e', norm(V'*V-eye(n),'fro')));
    R = chk(R,'cjacobi reconstruction', relerr(A, V*diag(d)*V') < 1e-12, sprintf('%.3e', relerr(A,V*diag(d)*V')));
    R = chk(R,'cjacobi eigenvalues descending', all(diff(d) <= 1e-10*max(abs(d))), '');
    R = chk(R,'cjacobi eigenvalues positive (SPD)', all(d > 0), sprintf('min=%.3e', min(d)));
catch e
    R = bad(R,'cjacobi', e);
end

% --- mp_pjacobi (needs mp for "mp3"; "mp2" also uses mp internally) -----
if hasmp
    for meth = ["mp2","mp3"]
        try
            % BOUND/SCOND are documented as mp3-only, so ask for them only there
            if strcmp(meth,"mp3")
                [V,D,NROT,NSWEEP,BOUND,SCOND] = mp_pjacobi(A, meth);
            else
                [V,D,NROT,NSWEEP] = mp_pjacobi(A, meth);
            end
            d = asvec(D); V = double(V); d = double(d);
            R = chk(R,sprintf('mp_pjacobi %s V orthogonal', meth), ...
                norm(V'*V - eye(n),'fro') < 1e-12, sprintf('%.3e', norm(V'*V-eye(n),'fro')));
            R = chk(R,sprintf('mp_pjacobi %s reconstruction', meth), ...
                relerr(A, V*diag(d)*V') < 1e-11, sprintf('%.3e', relerr(A,V*diag(d)*V')));
            R = chk(R,sprintf('mp_pjacobi %s eigenvalues descending', meth), ...
                all(diff(d) <= 1e-10*max(abs(d))), '');
            if strcmp(meth,"mp3")
                R = chk(R,'mp_pjacobi mp3 returns BOUND/SCOND', ...
                    isfinite(double(BOUND)) && isfinite(double(SCOND)) && double(SCOND) > 0, ...
                    sprintf('BOUND=%.3e SCOND=%.3e', double(BOUND), double(SCOND)));
            end
        catch e
            R = bad(R,sprintf('mp_pjacobi %s', meth), e);
        end
    end
else
    fprintf('  SKIP mp_pjacobi (Advanpix mp not on path)\n');
end

% --- the accuracy claim: mp3 beats eig on relative eigenvalue error -----
if hasmp
    try
        mp.Digits(34);
        dref = sort(double(eig(mp(A))), 'descend');
        [~,D3] = mp_pjacobi(A,"mp3");   d3 = sort(double(asvec(D3)),'descend');
        dml   = sort(eig(A), 'descend');
        e3 = max(abs(d3 - dref)./abs(dref));
        em = max(abs(dml - dref)./abs(dref));
        R = chk(R,'mp3 relative error < eig relative error', e3 < em, ...
            sprintf('mp3=%.3e eig=%.3e', e3, em));
    catch e
        R = bad(R,'mp3 vs eig accuracy', e);
    end
end

% --- get_reference / compute_error_evecs / compute_error_evals ---------
if hasmp
    try
        [Vref, S] = get_reference(A, 'E');
        R = chk(R,'get_reference E sizes', size(Vref,1)==n && numel(S)==n, ...
            sprintf('%dx%d, %d vals', size(Vref,1), size(Vref,2), numel(S)));
        [Vc,~] = eigsort(A);
        [err, relgap] = compute_error_evecs(A, Vc, 'E', Vref, S);
        R = chk(R,'compute_error_evecs shape', numel(err)==n && numel(relgap)==n, sprintf('%d,%d', numel(err), numel(relgap)));
        R = chk(R,'compute_error_evecs in [0,1]', all(err >= 0) && all(err <= 1+1e-10), ...
            sprintf('max=%.3e', max(err)));
        % exact vectors against themselves must give ~0 error
        [err0,~] = compute_error_evecs(A, double(Vref), 'E', Vref, S);
        R = chk(R,'compute_error_evecs(ref,ref) ~ 0', max(err0) < 1e-12, sprintf('%.3e', max(err0)));
    catch e
        R = bad(R,'get_reference/compute_error_evecs', e);
    end
end
try
    dref = sort(eig(A),'descend');
    e0 = compute_error_evals(dref, dref);
    R = chk(R,'compute_error_evals(x,x) ~ 0', max(double(e0)) < 1e-14, sprintf('%.3e', max(double(e0))));
    e1 = compute_error_evals(dref, dref*(1+1e-6));
    R = chk(R,'compute_error_evals detects 1e-6 perturbation', ...
        max(double(e1)) > 1e-7 && max(double(e1)) < 1e-5, sprintf('%.3e', max(double(e1))));
catch e
    R = bad(R,'compute_error_evals', e);
end

% --- get_testmatrix (KMS / Lehmer) -------------------------------------
try
    K = get_testmatrix(1); L = get_testmatrix(2);
    R = chk(R,'get_testmatrix(1) KMS symmetric', issymmetric(K), sprintf('%dx%d', size(K,1), size(K,2)));
    R = chk(R,'get_testmatrix(2) Lehmer symmetric', issymmetric(L), sprintf('%dx%d', size(L,1), size(L,2)));
    R = chk(R,'both test matrices SPD', min(eig(K)) > 0 && min(eig(L)) > 0, ...
        sprintf('minK=%.3e minL=%.3e', min(eig(K)), min(eig(L))));
catch e
    R = bad(R,'get_testmatrix', e);
end

% --- scond --------------------------------------------------------------
try
    R = chk(R,'scond(A,D) positive finite', isfinite(scond(A,'D')) && scond(A,'D') > 0, sprintf('%.3e', scond(A,'D')));
catch e
    R = bad(R,'scond', e);
end

nfail = summary(R, 'Chapter 4');
end

% ---------------------------------------------------------------- helpers
function R = chk(R, name, cond, detail)
R{end+1} = struct('name',char(name),'pass',logical(cond),'detail',char(detail));
fprintf('  [%s] %-46s %s\n', tf(cond), char(name), char(detail));
end
function R = bad(R, name, e)
R{end+1} = struct('name',char(name),'pass',false,'detail',e.message);
fprintf('  [FAIL] %-46s ERROR: %s\n', char(name), e.message);
end
function s = tf(c), if c, s='PASS'; else, s='FAIL'; end, end
function r = relerr(A,B), r = norm(A-B,'fro')/max(norm(A,'fro'),realmin); end
function v = asvec(D), if isvector(D), v = D(:); else, v = diag(D); end, end
function nfail = summary(R, label)
p = sum(cellfun(@(x)x.pass, R)); t = numel(R); nfail = t - p;
fprintf('  ---- %s: %d/%d passed, %d failed ----\n', label, p, t, nfail);
end
