function [U,S,V,nos,scalecond] = mposj_svec(G, nop)
%MPOSJ_SVEC Mixed-precision Jacobi SVD with QR after preconditioning.
%
%   Usage:
%       [U,S,V] = mposj_svec(G)
%       [U,S,V,nos] = mposj_svec(G)
%       [U,S,V,nos,scalecond] = mposj_svec(G)
%       [U,S,V] = mposj_svec(G, nop)
%
%   Purpose:
%       Compute an economy SVD of a real matrix. A single-precision SVD
%       supplies approximate right singular vectors, which are orthogonalized
%       in double precision. Apply this preconditioner in 34-digit arithmetic
%       (nop=3) or double precision (nop=2), then form an economy QR
%       factorization in double precision. DGESVJ is applied to R' and the
%       resulting left and right singular vectors are mapped back to G.
%       Wide matrices are handled through the transpose with the same nop.
%
%   Arguments:
%       G   - Real double matrix of size m-by-n.
%       nop - Number of precisions: 3 (default) or 2.
%
%   Outputs:
%       U,S,V     - Economy SVD factors. With r=min(m,n), U is m-by-r,
%                   S is r-by-r diagonal, and V is n-by-r. The columns of
%                   U and V are numerically orthonormal.
%       nos       - Number of sweeps used by DGESVJ.
%       scalecond - Column-scaled condition number of the preconditioned
%                   matrix (of G' for wide G). For nop=3 this diagnostic is
%                   evaluated before the 34-digit product is converted to
%                   double precision; the returned scalar is double.
%
%   Author:
%       Zhengbo Zhou, Manchester, UK, Dec 2025
%


if nargin == 1
    nop = 3; % Use MP3SVDJacobi by default.
end
if (nop ~= 2) && (nop ~= 3)
    error("The number of precision (the second argument) should be 2 or 3.");
end

[m,n] = size(G);
if m < n
    [Vt,St,Ut,nos,scalecond] = mposj_svec(G', nop);
    U = Ut; S = St; V = Vt;
    return
end
idty = eye(n);

% Compute the preconditioner
[~,~,Vs] = svd(single(G),'econ');
[Vd,~] = qr(double(Vs));

% Apply the preconditioner
if nop == 3
    Gmp = mp(G, 34);
    Vdmp = mp(Vd, 34);
    Gtmp = Gmp*Vdmp;
    Gt = double(Gtmp);
else
    Gtmp = G*Vd;
    Gt = Gtmp;
end

% Use the product before demotion for the scaled-condition diagnostic.
if nargout >= 5
    scalecond = double(scond(Gtmp, 'C'));
end

% Apply the one-sided Jacobi

[Qgt,Rt] = qr(Gt,'econ');
optlwork = max(6,m+n);
work = zeros(optlwork,1);
[V,S,U,~,work,info] = dgesvj_mex(Rt','L','U','V',n,idty,optlwork,work);
nos = work(4);
U = Qgt*U;

if info < 0
    error("DGESVJ has invalid inputs.\n");
elseif info > 0
    warning("DGESVJ does not converged in 30 iterations.\n")
end

V = Vd*V;

end
