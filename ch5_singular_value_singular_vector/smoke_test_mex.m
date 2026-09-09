function nfail = smoke_test_mex()
%SMOKE_TEST_MEX - Chapter 5 MEX wrapper smoke test
%
%   Usage:
%       nfail = smoke_test_mex()
%
%   Purpose:
%       Exercises the four LAPACK MEX wrappers in shared/ on small matrices,
%       in every JOBA mode the chapter's code actually uses:
%           dgesvj : 'G' (general), 'U' (upper triangular), 'L' (lower)
%           dgejsv : 'C'
%           sgesvj : 'G', 'U'   (single precision)
%           sgejsv : 'C'        (single precision)
%       and checks the returned singular values against MATLAB's svd.
%
%       RUN THIS IN ITS OWN MATLAB PROCESS. A mismatched OpenBLAS/libomp
%       will segfault MATLAB rather than raise an error, so isolating it
%       keeps a crash from taking down the other chapter tests.
%
%       Before calling, it reports which shared library each binary links,
%       because that is the usual cause of a crash here: OpenBLAS >= 0.3.31
%       from Homebrew pulls in libomp, which collides with MATLAB's own.
%
%   Author:
%       Smoke test for the thesis code reorganisation.

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'shared'));
R = {};
fprintf('=== Chapter 5: MEX wrappers ===\n');
fprintf('  MATLAB %s on %s\n', version('-release'), computer('arch'));

names = {'dgesvj_mex','dgejsv_mex','sgesvj_mex','sgejsv_mex'};

% --- stage 1: are they present and recognised as MEX? ------------------
for i = 1:numel(names)
    w = which(names{i});
    R = chk(R, sprintf('%s found', names{i}), exist(names{i},'file') == 3, w);
end

% --- stage 2: double precision ----------------------------------------
rng(0);
m = 50; n = 35;
A  = gallery('randsvd', [m,n], 1e6, 3);
hasmp = ~isempty(which('mp'));
if hasmp, mp.Digits(34); sref = double(svd(mp(A))); else, sref = svd(A); end
fprintf('  reference: %s\n', ternary(hasmp,'Advanpix mp (34 digits)','MATLAB svd (NOT accurate for small sv)'));

% dgesvj, JOBA = 'G'  (as used by the svals_*/svecs_* scripts)
try
    [U,S,V,sva,work,info] = dgesvj_mex(A,'G','U','V',n,eye(n),max(6,m+n));
    R = svdchk(R,'dgesvj G', A, U, S, V, sref, info, 1e-12);
catch e, R = bad(R,'dgesvj G', e); end

% dgesvj, JOBA = 'U'  (as used by mposj on the QR factor Rt)
try
    [~,Rt] = qr(A,'econ');                 % upper triangular, n-by-n
    optlwork = max(6,m+n); work = zeros(optlwork,1);
    [U,S,V,~,work,info] = dgesvj_mex(Rt,'U','U','V',n,eye(n),optlwork,work);
    R = svdchk(R,'dgesvj U (mposj path)', Rt, U, S, V, refsv(Rt,hasmp), info, 1e-12);
catch e, R = bad(R,'dgesvj U', e); end

% dgesvj, JOBA = 'L'  (as used by mposj_svec on Rt')
try
    [~,Rt] = qr(A,'econ');
    L = Rt';                                % lower triangular
    optlwork = max(6,m+n); work = zeros(optlwork,1);
    [V2,S2,U2,~,work,info] = dgesvj_mex(L,'L','U','V',n,eye(n),optlwork,work);
    R = svdchk(R,'dgesvj L (mposj_svec path)', L, V2, S2, U2, refsv(L,hasmp), info, 1e-12);
catch e, R = bad(R,'dgesvj L', e); end

% dgejsv, JOBA = 'C'
try
    [U,S,V,sva,work,iwork,info] = dgejsv_mex(A,'C','U','V','R','N','N');
    R = svdchk(R,'dgejsv C', A, U, S, V, sref, info, 1e-12);
catch e, R = bad(R,'dgejsv C', e); end

% --- stage 3: single precision ----------------------------------------
As = single(A);
if hasmp, srefs = double(svd(mp(double(As)))); else, srefs = svd(double(As)); end
try
    [U,S,V,sva,work,info] = sgesvj_mex(As,'G','U','V',n,eye(n,'single'),max(6,m+n));
    R = svdchk(R,'sgesvj G', double(As), double(U), double(S), double(V), srefs, info, 1e-5, double(svd(single(A))));
catch e, R = bad(R,'sgesvj G', e); end

try
    [~,Rs] = qr(As,'econ');
    optlwork = max(6,m+n); work = zeros(optlwork,1,'single');
    [U,S,V,~,work,info] = sgesvj_mex(Rs,'U','U','V',n,eye(n,'single'),optlwork,work);
    R = svdchk(R,'sgesvj U (mposj_ssd path)', double(Rs), double(U), double(S), double(V), refsv(double(Rs),hasmp), info, 1e-5, double(svd(Rs)));
catch e, R = bad(R,'sgesvj U', e); end

try
    [U,S,V,sva,work,iwork,info] = sgejsv_mex(As,'C','U','V','R','N','N');
    R = svdchk(R,'sgejsv C', double(As), double(U), double(S), double(V), srefs, info, 1e-5, double(svd(single(A))));
catch e, R = bad(R,'sgejsv C', e); end

% --- stage 4: the accuracy claim these wrappers exist to support -------
% On an ill-conditioned matrix DGESVJ should beat MATLAB's svd on the
% relative error of the SMALL singular values.
try
    B = gallery('randsvd', [60,40], 1e12, 3);
    if hasmp, bref = double(svd(mp(B))); else, bref = svd(B); end
    [~,Sj,~,~,~,infoj] = dgesvj_mex(B,'G','U','V',40,eye(40),max(6,100));
    sj = sort(asvec(Sj),'descend');
    sm = svd(B);
    ej = max(abs(sj - bref)./bref);
    em = max(abs(sm - bref)./bref);
    R = chk(R,'dgesvj info==0 on kappa=1e12', infoj == 0, sprintf('info=%d', infoj));
    if hasmp
        R = chk(R,'dgesvj beats svd at kappa=1e12', ej < em, sprintf('dgesvj=%.3e  matlab svd=%.3e', ej, em));
    else
        R = chk(R,'dgesvj rel err finite (no mp ref)', isfinite(ej), sprintf('%.3e', ej));
    end
catch e, R = bad(R,'dgesvj ill-conditioned', e); end

nfail = summary(R, 'Chapter 5 MEX');
end

% ---------------------------------------------------------------- helpers
function R = svdchk(R, name, A, U, S, V, sref, info, tol, smlin)
% sref must be an ACCURATE reference (Advanpix mp where available).
% MATLAB's own svd is not a valid reference here: it loses relative
% accuracy on the small singular values of an ill-conditioned matrix,
% which is exactly what these Jacobi routines exist to avoid. So the
% singular-value check is stated as "no worse than MATLAB's svd", which
% is self-calibrating and is the property the thesis actually claims.
s = sort(asvec(S),'descend');
sref = sort(sref(:),'descend');
k = min(numel(s), numel(sref));
relsv  = max(abs(s(1:k) - sref(1:k)) ./ max(sref(1:k), realmin));
if nargin < 10 || isempty(smlin), smlin = svd(A); end
sml    = sort(double(smlin),'descend');
relml  = max(abs(sml(1:k) - sref(1:k)) ./ max(sref(1:k), realmin));
bwd    = norm(A - U*diagm(S)*V','fro')/norm(A,'fro');
orth   = norm(U'*U - eye(size(U,2)),'fro');
R = chk(R, [name ' info==0'], info == 0, sprintf('info=%d', info));
R = chk(R, [name ' sv no worse than svd'], relsv <= max(relml*10, tol), ...
    sprintf('mex=%.3e  matlab=%.3e', relsv, relml));
R = chk(R, [name ' backward error'],  bwd  < tol*1e2, sprintf('%.3e', bwd));
R = chk(R, [name ' U orthogonality'], orth < tol*1e2, sprintf('%.3e', orth));
end
function D = diagm(S), if isvector(S), D = diag(S); else, D = S; end, end
function v = asvec(S), if isvector(S), v = S(:); else, v = diag(S); end, end
function R = chk(R, name, cond, detail)
R{end+1} = struct('name',char(name),'pass',logical(cond),'detail',char(detail));
fprintf('  [%s] %-38s %s\n', tf(cond), char(name), char(detail));
end
function R = bad(R, name, e)
R{end+1} = struct('name',char(name),'pass',false,'detail',e.message);
fprintf('  [FAIL] %-38s ERROR: %s\n', char(name), e.message);
end
function s = tf(c), if c, s='PASS'; else, s='FAIL'; end, end
function o = ternary(c,a,b), if c, o=a; else, o=b; end, end
function s = refsv(A,hasmp)
if hasmp, mp.Digits(34); s = double(svd(mp(A))); else, s = svd(A); end
end
function nfail = summary(R, label)
p = sum(cellfun(@(x)x.pass, R)); t = numel(R); nfail = t - p;
fprintf('  ---- %s: %d/%d passed, %d failed ----\n', label, p, t, nfail);
end
