function nfail = smoke_test()
%SMOKE_TEST - Chapter 5 smoke test (algorithms and utilities)
%
%   Usage:
%       nfail = smoke_test()
%
%   Purpose:
%       Exercises every function in ch5_singular_value_singular_vector/shared
%       except the raw MEX wrappers, which have their own isolated test:
%           smoke_test_mex
%       Run that FIRST. Everything here depends on the MEX wrappers, so if
%       they are broken this test cannot mean anything.
%
%       Covers mposj (Algorithm: mp precond one-sided Jacobi), mposj_svec
%       (the transposed variant used for singular vectors), mposj_ssd, the
%       two compute_error variants, both reference routines and both
%       test-matrix generators.
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
hasam = ~isempty(which('anymatrix'));

fprintf('=== Chapter 5: singular values and singular vectors ===\n');
fprintf('  Advanpix mp: %d   Anymatrix: %d\n', hasmp, hasam);
rng(0);
m = 60; n = 40;
A = gallery('randsvd', [m,n], 1e8, 3);

% --- scond --------------------------------------------------------------
try
    sC = scond(A,'C');
    R = chk(R,'scond(A,C) positive finite', isfinite(sC) && sC > 0, sprintf('%.3e', sC));
catch e, R = bad(R,'scond', e); end

% --- reference_singular_values (71 digits) ------------------------------
if hasmp
    try
        sref = reference_singular_values(A);
        R = chk(R,'reference_singular_values length', numel(sref) == n, sprintf('%d', numel(sref)));
        R = chk(R,'reference_singular_values descending', all(diff(double(sref)) <= 0), '');
        R = chk(R,'reference kappa ~ 1e8', ...
            abs(log10(double(sref(1)/sref(end))) - 8) < 0.5, ...
            sprintf('kappa=%.3e', double(sref(1)/sref(end))));
    catch e, R = bad(R,'reference_singular_values', e); end
else
    fprintf('  SKIP reference_singular_values (needs mp)\n');
end

% --- mposj : the singular-value algorithm -------------------------------
if hasmp
    for nop = [2 3]
        try
            [U,S,V,nos,scalecond,timing] = mposj(A, nop, true);
            s = asvec(S);
            R = chk(R,sprintf('mposj nop=%d backward error', nop), ...
                relerr(A, U*diagm(S)*V') < 1e-13, sprintf('%.3e', relerr(A,U*diagm(S)*V')));
            R = chk(R,sprintf('mposj nop=%d U orthogonal', nop), ...
                norm(U'*U - eye(n),'fro') < 1e-12, sprintf('%.3e', norm(U'*U-eye(n),'fro')));
            R = chk(R,sprintf('mposj nop=%d V orthogonal', nop), ...
                norm(V'*V - eye(n),'fro') < 1e-12, sprintf('%.3e', norm(V'*V-eye(n),'fro')));
            R = chk(R,sprintf('mposj nop=%d sv descending', nop), all(diff(s) <= 1e-10*max(s)), '');
            R = chk(R,sprintf('mposj nop=%d scalecond/timing', nop), ...
                isfinite(scalecond) && scalecond > 0 && numel(timing) == 4 && all(timing >= 0), ...
                sprintf('scond=%.3e nos=%d t=%.3fs', scalecond, nos, sum(timing)));
        catch e, R = bad(R,sprintf('mposj nop=%d', nop), e); end
    end
    % wantscond default must be off, and scalecond then empty
    try
        [~,~,~,~,sc0] = mposj(A);
        R = chk(R,'mposj scalecond opt-in (empty by default)', isempty(sc0), sprintf('%s', mat2str(size(sc0))));
    catch e, R = bad(R,'mposj default wantscond', e); end
    % the accuracy claim: mposj beats MATLAB svd in relative sv error
    try
        sref = double(reference_singular_values(A));
        [~,S3,~] = mposj(A,3,false);  s3 = sort(asvec(S3),'descend');
        sm = sort(svd(A),'descend');
        e3 = max(abs(s3 - sref)./sref);  em = max(abs(sm - sref)./sref);
        R = chk(R,'mposj beats MATLAB svd (rel sv err)', e3 < em, sprintf('mposj=%.3e svd=%.3e', e3, em));
    catch e, R = bad(R,'mposj accuracy claim', e); end
else
    fprintf('  SKIP mposj (needs mp)\n');
end

% --- mposj_svec : the transposed algorithm for singular vectors ---------
if hasmp
    try
        [U,S,V,nos,scalecond] = mposj_svec(A);
        s = asvec(S);
        R = chk(R,'mposj_svec backward error', relerr(A, U*diagm(S)*V') < 1e-13, ...
            sprintf('%.3e', relerr(A,U*diagm(S)*V')));
        R = chk(R,'mposj_svec U orthogonal', norm(U'*U - eye(n),'fro') < 1e-12, ...
            sprintf('%.3e', norm(U'*U-eye(n),'fro')));
        R = chk(R,'mposj_svec V orthogonal', norm(V'*V - eye(n),'fro') < 1e-12, ...
            sprintf('%.3e', norm(V'*V-eye(n),'fro')));
        R = chk(R,'mposj_svec sv descending', all(diff(s) <= 1e-10*max(s)), '');
        R = chk(R,'mposj_svec returns scalecond', isfinite(scalecond) && scalecond > 0, ...
            sprintf('%.3e (nos=%d)', scalecond, nos));
    catch e, R = bad(R,'mposj_svec', e); end
    % mposj and mposj_svec are different algorithms but must agree on sv
    try
        [~,Sa,~] = mposj(A,3,false);  [~,Sb,~] = mposj_svec(A);
        sa = sort(asvec(Sa),'descend'); sb = sort(asvec(Sb),'descend');
        R = chk(R,'mposj vs mposj_svec agree on sv', max(abs(sa-sb)./sa) < 1e-12, ...
            sprintf('%.3e', max(abs(sa-sb)./sa)));
    catch e, R = bad(R,'mposj vs mposj_svec', e); end
else
    fprintf('  SKIP mposj_svec (needs mp)\n');
end

% --- mposj_ssd : single+double, no mp needed ---------------------------
try
    [U,S,V,nos,scalecond] = mposj_ssd(A, true);
    s = asvec(S);
    R = chk(R,'mposj_ssd backward error', relerr(A, U*diagm(S)*V') < 1e-5, ...
        sprintf('%.3e', relerr(A,U*diagm(S)*V')));
    R = chk(R,'mposj_ssd V orthogonal', norm(double(V)'*double(V) - eye(n),'fro') < 1e-4, ...
        sprintf('%.3e', norm(double(V)'*double(V)-eye(n),'fro')));
    R = chk(R,'mposj_ssd sv descending', all(diff(double(s)) <= 1e-6*max(double(s))), '');
    R = chk(R,'mposj_ssd scalecond opt-in', isfinite(double(scalecond)) && double(scalecond) > 0, ...
        sprintf('%.3e (nos=%d)', double(scalecond), nos));
    [~,~,~,~,sc0] = mposj_ssd(A);
    R = chk(R,'mposj_ssd scalecond empty by default', isempty(sc0), '');
catch e, R = bad(R,'mposj_ssd', e); end

% --- compute_error_svals ------------------------------------------------
if hasmp
    try
        sref = reference_singular_values(A);
        [Ue,Se,Ve] = svd(A,'econ');
        [f, r, oU, oV] = compute_error_svals(A, Ue, Se, Ve, sref);
        R = chk(R,'compute_error_svals outputs finite', ...
            isfinite(f) && isfinite(r) && isfinite(oU) && isfinite(oV), ...
            sprintf('f=%.3e r=%.3e oU=%.3e oV=%.3e', f, r, oU, oV));
        R = chk(R,'compute_error_svals backward/orth small for svd', ...
            r < 1e-14 && oU < 1e-13 && oV < 1e-13, sprintf('r=%.3e', r));
    catch e, R = bad(R,'compute_error_svals', e); end
end

% --- get_reference / compute_error_svecs, all three Vtypes -------------
if hasmp
    for vt = ['S','L']
        try
            [Vref, S] = get_reference(A, vt);
            [Uc,~,Vc] = svd(A,'econ');
            if vt == 'S', Vin = Vc; else, Vin = Uc; end
            [err, relgap] = compute_error_svecs(A, Vin, vt, Vref, S);
            R = chk(R,sprintf('compute_error_svecs %c shape', vt), ...
                numel(err) == n && numel(relgap) == n, sprintf('%d', numel(err)));
            R = chk(R,sprintf('compute_error_svecs %c in [0,1]', vt), ...
                all(err >= 0) && all(err <= 1+1e-10), sprintf('max=%.3e', max(err)));
            [err0,~] = compute_error_svecs(A, double(Vref), vt, Vref, S);
            R = chk(R,sprintf('compute_error_svecs %c (ref,ref) ~ 0', vt), ...
                max(err0) < 1e-12, sprintf('%.3e', max(err0)));
        catch e, R = bad(R,sprintf('compute_error_svecs %c', vt), e); end
    end
end

% --- test matrix generators --------------------------------------------
try
    K = get_testmatrix_kms_lehmer(1); L = get_testmatrix_kms_lehmer(2);
    R = chk(R,'get_testmatrix_kms_lehmer(1) KMS', issymmetric(K), sprintf('%dx%d', size(K,1), size(K,2)));
    R = chk(R,'get_testmatrix_kms_lehmer(2) Lehmer', issymmetric(L), sprintf('%dx%d', size(L,1), size(L,2)));
catch e, R = bad(R,'get_testmatrix_kms_lehmer', e); end

if hasam
    for id = 1:4
        try
            M = get_testmatrix_anymatrix(id);
            R = chk(R,sprintf('get_testmatrix_anymatrix(%d)', id), ...
                ~isempty(M) && all(isfinite(M(:))), sprintf('%dx%d cond=%.2e', size(M,1), size(M,2), cond(M)));
        catch e, R = bad(R,sprintf('get_testmatrix_anymatrix(%d)', id), e); end
    end
else
    fprintf('  SKIP get_testmatrix_anymatrix (Anymatrix not on path)\n');
end

nfail = summary(R, 'Chapter 5');
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
function r = relerr(A,B), r = norm(double(A)-double(B),'fro')/max(norm(double(A),'fro'),realmin); end
function v = asvec(S), if isvector(S), v = double(S(:)); else, v = double(diag(S)); end, end
function D = diagm(S), if isvector(S), D = diag(double(S)); else, D = double(S); end, end
function nfail = summary(R, label)
p = sum(cellfun(@(x)x.pass, R)); t = numel(R); nfail = t - p;
fprintf('  ---- %s: %d/%d passed, %d failed ----\n', label, p, t, nfail);
end
