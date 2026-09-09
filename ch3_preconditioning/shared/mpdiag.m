function [Qt, Dt] = mpdiag(A)
%mpdiag - Low precision spectral decomposition with orthogonal 
%eigenvector matrix
%
%   [Qt, Dt] = MPDIAG(A) essentially computes a spectral
%   decomposition at single precision but with orthogonal eigenvectors at
%   double precision.
%   
%   The algorithm make use of the LAPACK routines, 
%       SSYTRD: Tridiagonal reduction at single precision. 
%       DORGTR: Accumulation of the Householder vectors. 
%       DSTEQR: Diagonalization of a tridiagonal matrix. 
%   These routines are called in MATLAB by using the "lapack" package. 
%
%   Input arguments:
%       A  : A symmetric matrix.
%
%   Output arguments:
%       Qt : Eigenvector matrix such that Qt'*A*Qt = Dt. Importantly, here
%           Qt is orthogonal at double precision.
%       Dt : Eigenvalues of A with accuracy at single precision.

% Test input
if ~issymmetric(A)
    A = (A + A')/2;
    warning("Input is not symmetric, is this intentional?\n" + ...
        "I will use the nearest symmetric matrix of your input\n");
end
if class(A) ~= "double"
    A = double(A);
    warning("Your input is not stored in double, is this intentional?" + ...
        "I will cast your input to double\n")
end

% Parameters
n = size(A,1);

%  Use SSYTRD to perform tridiagonal reduction in single precision
UPLO = 'U';
T_diag = zeros(n,1);
T_offd = zeros(n-1,1);
tau = zeros(n-1,1);
work = zeros(1);
lwork = -1; % Looking for optimal work size
info = 0;
C = lapack('ssytrd', UPLO, n, A, n, T_diag, T_offd, tau, work, lwork, info);
info = C{10};
work = C{8};
if info == 0
    lwork = int32(work(1));
    work = ones(lwork,1);
else
    error("SSYTRD failed\n");
end
C = lapack('ssytrd', UPLO, n, A, n, T_diag, T_offd, tau, work, lwork, info);
info = C{10};
if info ~= 0, error("SSYTRD failed\n"); end
A = double(C{3});
T_diag = double(C{5});
T_offd = double(C{6});

% Recompute nontrivial reflectors; preserve LAPACK identity reflectors.
tau_ref = double(C{7});
tau = tau_ref;
for i = 2:n-1
    if tau_ref(i) == 0
        continue
    end
    v = zeros( n,1 );
    v( i ) = 1;
    v( 1:i-1 ) = A( 1:i-1,i+1 );
    tau( i ) = 2 / ( norm(v)^2 );
end

% check if we get the correct tau
if ( norm( tau_ref - tau ) >= n * float_params('s') )
    error("TAU in SSYTRD is not correctly computed.")
end

% Generate a real orthogonal matrix which is the defined as product of
% n-1 elementary reflectors using DORGTR.
C = lapack('dorgtr', UPLO, n, A, n, tau, work, -1, info);
info = C{8};
work = C{6};
if info == 0
    lwork = int32(work(1));
    work = ones(lwork,1);
else
    error("DORGTR failed\n");
end
C = lapack('dorgtr', UPLO, n, A, n, tau, work, lwork, info);
info = C{8};
if info ~= 0, error("DORGTR failed\n"); end
Q = C{3};

% Compute the eigenvalues and corresponding eigenvectors of A using a
% tridiagonal QR algorithm using DSTEQR.
compz = 'V';
C = lapack('dsteqr', compz, n, T_diag, T_offd, Q, n, zeros(2*n - 2,1), info);
info = C{8};
if info == 0
    Dt = diag(C{3});
    Qt = C{5};
elseif info < 0
    error("the %d-th argument in DSTEQR had an illegel value.\n", -info);
elseif info > 0
    warning("DSTEQR fails to converge in %d iterations.\n" + ...
        "The i-th off diagonal entry does not converge to 0.\n" + ...
        "We return the tridiagonal matrix", 30 * n, info);
    T_diag = C{3};
    T_offd = C{4};
    Dt = diag(T_diag) + diag(T_offd, 1) + diag(T_offd, -1);
    Qt = C{5};
end

end