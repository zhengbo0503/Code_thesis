function [U,S,V,nos,scalecond,timing] = mposj_complex(G,nop,wantscond)
%MPOSJ_COMPLEX - Mixed precision one-sided Jacobi SVD for complex matrices
%
%   Usage:
%       [U,S,V] = mposj_complex(G)
%       [U,S,V,nos] = mposj_complex(G)
%       [U,S,V] = mposj_complex(G, nop)
%       [U,S,V,nos,scalecond,timing] = mposj_complex(G, nop, true)
%
%   Purpose:
%       MPOSJ_COMPLEX computes the SVD of a general complex matrix using
%       ZGESVJ as the final one-sided Jacobi SVD routine. It constructs the
%       preconditioner at complex single precision, applies it at complex
%       double or 34-digit multiprecision, and performs the final Jacobi
%       computation at complex double precision.
%
%   Arguments:
%       (1) G - Complex, double matrix.
%       (2) nop - Integer, and by default 3.
%           If nop = 3, MPOSJ_COMPLEX applies the preconditioner at 34-digit
%           multiprecision; if nop = 2, it uses double precision.
%       (3) wantscond - Logical, and by default false.
%           Whether to compute the one-sided column-scaled condition number
%           of the preconditioned matrix. This diagnostic must be requested
%           explicitly.
%
%   Outputs:
%       (1) U,S,V - Economy-size SVD factors. If G is mxn and
%           r = min(m,n), then U is mxr, S is rxr, and V is nxr.
%       (2) nos - Integer giving the number of ZGESVJ sweeps.
%       (3) scalecond - Real, double scalar containing the scaled condition
%           number, or NaN when wantscond is false.
%       (4) timing - Real, double vector containing the time used for
%           constructing the preconditioner, applying the preconditioner,
%           applying the Jacobi algorithm, and everything else.
%

if nargin < 2 || isempty(nop)
    nop = 3;
end
if nargin < 3 || isempty(wantscond)
    wantscond = false;
end
if nop ~= 2 && nop ~= 3
    error('mposj_complex:InvalidPrecisionCount', ...
        'The second argument must be 2 or 3.');
end
if ~isnumeric(G) || ~isa(G,'double') || ~isfloat(G) || isreal(G)
    error('mposj_complex:InvalidInput', ...
        'G must be a complex double matrix.');
end

tTotal = tic;
scalecond = NaN;
algorithmMpDigits = 34;
[m,n] = size(G);

if m < n
    [Vt,St,Ut,nos,scalecond,timing] = mposj_complex(G',nop,wantscond);
    U = Ut;
    S = St;
    V = Vt;
    return
end

doqr = m >= (11*n)/6;

tPreconditioner = tic;
[~,~,Vs] = svd(single(G),'econ');
[Vd,~] = qr(double(Vs));
tPreconditioner = toc(tPreconditioner);

tApply = tic;
if nop == 3
    % Object-local precision prevents ambient mp.Digits() from changing
    % the algorithmic precision.
    Gtmp = mp(G,algorithmMpDigits)*mp(Vd,algorithmMpDigits);
    Gt = double(Gtmp);
else
    Gt = G*Vd;
    Gtmp = Gt;
end
tApply = toc(tApply);

if wantscond
    scalecond = double(scond_columns_complex(Gtmp));
end

tJacobi = tic;
if doqr
    [Qgt,Rt] = qr(Gt,'econ');
    [U,S,V,~,rwork,info] = zgesvj_mex(Rt,'U');
    U = Qgt*U;
else
    [U,S,V,~,rwork,info] = zgesvj_mex(Gt,'G');
end
tJacobi = toc(tJacobi);

if info < 0
    error('mposj_complex:ZGESVJInput','ZGESVJ reported an invalid input.');
elseif info > 0
    error('mposj_complex:ZGESVJNoConvergence', ...
        'ZGESVJ did not converge in 30 sweeps (INFO=%d).',info);
end

nos = double(rwork(4));
V = Vd*V;

tTotal = toc(tTotal);
tOther = tTotal-tPreconditioner-tApply-tJacobi;
timing = [tPreconditioner,tApply,tJacobi,tOther];
end
