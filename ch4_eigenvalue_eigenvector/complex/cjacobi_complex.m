function [V,D,NROT,NSWEEP,INFO] = cjacobi_complex(A, EVEC)
%CJACOBI_COMPLEX - Cyclic-by-row Jacobi algorithm
%   
%   Usage:
%       [V,D,NROT,NSWEEP,INFO] = cjacobi_complex(A,EVEC)
%   
%   Purpose:
%       CJACOBI_COMPLEX computes the spectral decomposition (SPD) of a Hermitian, 
%       N-by-N matrix A using the cyclic-by-row Jacobi algorithm. 
%       The SPD of A is written as 
%               A = V * diag(D) * V'
%       where V is an N-by-N unitary matrix, and D is an N dimensional
%       array whose entries are eigenvalues in descending order.
%       CJACOBI_COMPLEX can sometimes compute tiny eigenvalues with high relative
%       accuracy if A is Hermitian positive definite. 
%
%   Input:
%    - A is COMPLEX matrix, dimension (N,N)
%       A is a HERMITIAN matrix. The function will terminate if the
%       input is not Hermitian. Real symmetric input is also accepted. 
%    
%    - EVEC is INTEGER
%       Default value is 1;
%        = 1: CJACOBI_COMPLEX computes both eigenvectors and eigenvalues of A.
%        = 0: CJACOBI_COMPLEX only computes eigenvalues of A. 
%
%   Output:
%    - D is REAL vector, dimension (N)
%       If INFO = 0, D(i) is the ith largest eigenvalue of A.
%
%    - V is COMPLEX matrix, dimension (N,N)
%       If INFO = 0 and EVEC = 1, V(:,i) is an normalized eigenvector 
%       corresponding to D(i). Loss of orthogonality/unitarity depends
%       on floating-point roundoff accumulated over the computation.
%       If EVEC = 0, V is a permuted identity placeholder, not eigenvectors. 
%   
%    - NROT is INTEGER 
%       The number of applied Jacobi rotation.
%
%    - NSWEEP is INTEGER 
%       The number of sweep the Jacobi rotation required to converge.
%       One sweep indicates the Jacobi algorithm has run through all A(i,j)
%       where i < j.
%
%    - INFO is INTEGER 
%       =  0: sucessfully exit
%       = -1: The Jacobi algorithm does not converge in 30 iterations.
%
%   Author:
%       Zhengbo Zhou, June 2025, Manchester, UK
%   
%   Reference 
%   Demmel & Veselic, Jacobi's method is more accurate than QR, 
%   SIAM J. Matrix Anal. Appl. 13, 1204-1245, 1992

% Minimal Hermitian adaptation of ../shared/cjacobi.m.
% Parameters
n = size(A,1);
u = float_params('d');
NROT = 0;
NSWEEP = 0;
V = eye(n);
done_rot = true;
maxiter = 30;
tol = sqrt(n) * u;

% Assign EVEC default value 
if nargin == 1
    EVEC = 1;
end

% Check Hermitian symmetry (as in the real routine, no definiteness test).
if ~ishermitian(A)
    error('cjacobi_complex:NotHermitian', 'Matrix must be HERMITIAN.');
end

while done_rot && NSWEEP < maxiter
    NSWEEP = NSWEEP + 1;
    done_rot = false; 
    
    % Main iteration 
    for p = 1:n-1
        for q = p+1:n 
            if abs(A(p,q)) > tol * sqrt(abs(A(p,p)*A(q,q)))
                done_rot = true;
                NROT = NROT + 1;

                % Form Jacobi rotation
                [c,s] = jacobi_pair_complex(A,p,q);
                J = [c,s;-conj(s),c];

                % Apply Jacobi rotation to A
                A([p,q],:) = J'*A([p,q],:);
                A(:,[p,q]) = A(:,[p,q]) * J;
                % Hermitian diagonal entries are real; discard imaginary roundoff.
                A(p,p) = real(A(p,p));
                A(q,q) = real(A(q,q));
                A(p,q) = 0;
                A(q,p) = 0;

                % Conditionally apply the Jacobi rotation to V
                if (EVEC) 
                    V(:,[p,q]) = V(:,[p,q]) * J;
                end
            end
        end
    end
end

% Extract eigenvalues of A
D = diag(A);

% Sort eigenvalues in descending order 
[~,ind] = sort(D,'descend');
D = D(ind); 
V = V(:,ind);

% INFO
if (NSWEEP >= maxiter)
    INFO = -1;
else
    INFO = 0;
end

end
