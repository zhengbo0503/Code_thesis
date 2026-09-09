function [U,S,V,nos,scalecond,timing] = mposj_ssd_complex(G,wantscond)
%MPOSJ_SSD_COMPLEX - Complex (single, single, double) Jacobi SVD algorithm
%
%   Usage:
%       [U,S,V] = mposj_ssd_complex(G)
%       [U,S,V,nos] = mposj_ssd_complex(G)
%       [U,S,V,nos,scalecond] = mposj_ssd_complex(G, true)
%       [U,S,V,nos,scalecond,timing] = mposj_ssd_complex(G, wantscond)
%
%   Purpose:
%       MPOSJ_SSD_COMPLEX computes the SVD of a general complex matrix using
%       CGESVJ as the final one-sided Jacobi SVD routine. It constructs the
%       preconditioner at complex single precision, applies it at complex
%       double precision, demotes the preconditioned matrix to complex
%       single precision, and performs the final Jacobi computation there.
%
%   Arguments:
%       (1) G - Complex, single matrix.
%       (2) wantscond - Logical, and by default false.
%           Whether to compute the one-sided column-scaled condition number
%           of the double-precision preconditioned matrix.
%
%   Outputs:
%       (1) U,S,V - Economy-size single-precision SVD factors. If G is mxn
%           and r = min(m,n), then U is mxr, S is rxr, and V is nxr.
%       (2) nos - Integer giving the number of CGESVJ sweeps.
%       (3) scalecond - Real, double scalar containing the scaled condition
%           number, or NaN when wantscond is false.
%       (4) timing - Real, double vector containing the time used for
%           constructing the preconditioner, applying the preconditioner,
%           applying the Jacobi algorithm, and everything else. Timing is
%           measured only when this sixth output is requested.
%

if nargin < 2 || isempty(wantscond)
    wantscond = false;
end
if ~isnumeric(G) || ~isa(G,'single') || ~isfloat(G) || isreal(G)
    error('mposj_ssd_complex:InvalidInput', ...
        'G must be a complex single matrix.');
end

measureTiming = nargout >= 6;
if measureTiming
    tTotalStart = tic;
end
scalecond = NaN;
[m,n] = size(G);

if m < n
    if measureTiming
        [Vt,St,Ut,nos,scalecond,timing] = mposj_ssd_complex(G',wantscond);
    else
        [Vt,St,Ut,nos,scalecond] = mposj_ssd_complex(G',wantscond);
        timing = [];
    end
    U = Ut;
    S = St;
    V = Vt;
    return
end

doqr = m >= (11*n)/6;

if measureTiming
    tPreconditionerStart = tic;
end
[~,~,Vt] = svd(G,'econ');
if measureTiming
    tPreconditioner = toc(tPreconditionerStart);
end

if measureTiming
    tApplyStart = tic;
end
Gtmp = double(G)*double(Vt);
Gt = single(Gtmp);
if measureTiming
    tApply = toc(tApplyStart);
end

if wantscond
    scalecond = double(scond_columns_complex(Gtmp));
end

if measureTiming
    tJacobiStart = tic;
end
if doqr
    [Qgt,Rt] = qr(Gt,'econ');
    [U,S,V,~,rwork,info] = cgesvj_mex(Rt,'U');
    U = Qgt*U;
else
    [U,S,V,~,rwork,info] = cgesvj_mex(Gt,'G');
end
if measureTiming
    tJacobi = toc(tJacobiStart);
end

if info < 0
    error('mposj_ssd_complex:CGESVJInput','CGESVJ reported an invalid input.');
elseif info > 0
    error('mposj_ssd_complex:CGESVJNoConvergence', ...
        'CGESVJ did not converge in 30 sweeps (INFO=%d).',info);
end

nos = double(rwork(4));
V = Vt*V;

if measureTiming
    tTotal = toc(tTotalStart);
    tOther = tTotal-tPreconditioner-tApply-tJacobi;
    timing = [tPreconditioner,tApply,tJacobi,tOther];
else
    timing = [];
end
end
