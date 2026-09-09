function [U,S,V,nos,scalecond,timing] = mposj(G, nop, wantscond)
%MPOSJ - Mixed precision one-sided Jacobi algorithm
%		 using quadruple, double and single precisions
%
%   Usage:
%       [U,S,V] = mposj(G)
%       [U,S,V,nos] = mposj(G)
%       [U,S,V,nos,scalecond,timing] = mposj(G, nop, true)
%       [U,S,V] = mposj(G, nop)
%
%   Purpose: 
%       MPOSJ computes a SVD of a general matrix based on LAPACK routine
%       DGESVJ, the one-sided Jacobi SVD algorithm. 
%       MPOSJ first computes a preconditioner at single precision, applies
%       the preconditioner at double or quadruple precision, and finally
%       performs DGESVJ on the preconditioned matrix at double precision. 
%
%		For ill-conditioned matrix, MPOSJ can compute singular values with
%       smaller relative forward error compared to MATLAB function svd.
%
%   Arguments: 
%       (1) G - Real, double matrix.
%       (2) nop - Integer, and by default 3.
%           If nop = 3, MPOSJ will use quadruple precision to apply the
%           preconditioner; If nop = 2, MPOSJ will use double precision
%           instead. 
%       (3) wantscond - Logical, and by default false.
%           Whether to compute the scaled condition number of the
%           preconditioned matrix. This is a diagnostic that is not part
%           of the algorithm, and computing it is expensive when nop = 3
%           (it requires a quadruple precision SVD), so it must be
%           requested explicitly and is skipped in the timing tests.
%   
%   Outputs:
%       (1) U,S,V - Real, double matrix.
%           Suppose G is mxn, then U is mxn whose columns are numerically
%           orthonormal, S is a nxn diagonal matrix whose diagonal entries
%           are singular values of G, and V is a nxn numerically orthogonal
%           matrix. 
%       (2) nos - Integer.
%           Number of sweep used by DGESVJ.
%       (3) scalecond - Real, double.
%           Scaled condition number of the preconditioned matrix, using
%           the one-sided (column) scaling. Useful for some posterior
%           analysis. Empty unless wantscond is true.
%		(4) timing - Real, double vector
%			The runtime measure for four different components:
%			  (i) constructing the preconditioner,
%			 (ii) applying the preconditioner,
%			(iii) applying the Jacobi algorithm, and
%			 (iv) everything else.
%
%   Author:
%       Zhengbo Zhou, Manchester, UK, Dec 2025
%

t_total = tic;

if nargin < 2 || isempty(nop)
    nop = 3; % Use MP3SVDJacobi by default.
end
if (nop ~= 2) && (nop ~= 3)
    error("The number of precision (the second argument) should be 2 or 3.");
end
if nargin < 3 || isempty(wantscond)
    wantscond = false; % Diagnostic only, so opt in explicitly.
end

scalecond = [];

[m,n] = size(G);
if m < n   
    [Vt,St,Ut,nos,scalecond,timing] = mposj(G', nop, wantscond);
    U = Ut; S = St; V = Vt;
    return
end

idty = eye(n);
doqr = (m >= (11*n)/6); % Consistent with LAPACK choice. 

t_construct_preconditioner = tic;

% Compute the preconditioner
[~,~,Vs] = svd(single(G),'econ');
[Vd,~] = qr(double(Vs));

t_construct_preconditioner = toc(t_construct_preconditioner);

t_apply_preconditioner = tic;

% Apply the preconditioner
if nop == 3
    Gmp = mp(G); 
    Vdmp = mp(Vd); 
    Gtmp = Gmp*Vdmp; 
    Gt = double(Gtmp); 
else
    Gt = G*Vd;
    Gtmp = Gt;
end

t_apply_preconditioner = toc(t_apply_preconditioner);

% Output scaled condition number for posterior analysis. 
if wantscond
    scalecond = double(scond(Gtmp,'C')); 
end

t_Jacobi = tic; 

% Apply the one-sided Jacobi
if doqr % Apply QR factorization before doing Jacobi SVD
    [Qgt,Rt] = qr(Gt,'econ');
    optlwork = max(6,m+n);
    work = zeros(optlwork,1);
    [U,S,V,~,work,info] = dgesvj_mex(Rt,'U','U','V',n,idty,optlwork,work);
    nos = work(4);
    U = Qgt*U;
else % Plain Jacobi SVD
    optlwork = max(6,m+n);
    work = zeros(optlwork,1);
    [U,S,V,~,work,info] = dgesvj_mex(Gt,'G','U','V',n,idty,optlwork,work);
    nos = work(4);
end

t_Jacobi = toc(t_Jacobi);

if info < 0
    error("DGESVJ has invalid inputs.\n");
elseif info > 0
    warning("DGESVJ does not converged in 30 iterations.\n")
end

V = Vd*V;
t_total = toc(t_total);
t_else = t_total - t_construct_preconditioner - t_apply_preconditioner - t_Jacobi;

timing = [t_construct_preconditioner, t_apply_preconditioner, t_Jacobi, t_else];

end
