function [U,H,k] = ns(X, info)
%NS - Polar decomposition using the Newton-Schulz iteration
%   [U,H,k] = ns(X, info) returns the polar decomposition X=UH, where U is
%   unitary and H is positive definite. Optionally, the function can return
%   a positive integer k which is the no. of iteration required for the
%   Newton-Schulz iteration to converge. 
%   The function converge for any matrix as it applied on a scaled version
%   of X, which Xscale = X/norm(X,2). In such way, the Newton-Schulz
%   iteration always converge regardless of the size of sigma_1(X). 
%
%   Input argument:
%   X : A computed eigenvector matrix at single precision. 
%   info : Anything
%   
%   Output argument: 
%   U : The unitary polar factor of X. 
%   H : The Hermitian polar factor of X. 
%   k : The number of iteration required for the Newton-Schulz iteration to
%       converge. 
%   
%   Usage:
%   1. U = ns(X) computes the unitary polar factor of X.
%   2. U = ns(X, info) prints out the number of iteration required by
%       orthogonalizing X. Here info can be anything, since whether the
%       function prints out the iteration number is based on the test on
%       the number of inputs. 
%   3. [U,H] = ns(X) computes the polar decomposition X = UH. 
%   4. [U,H,k] = ns(X) computes the polar decomposition X = UH, and outputs
%       the no. of iterations k.

% Parameters 
u = float_params('d');
n = size(X,1);
cgtol = sqrt(sqrt(n)*u);
I = eye(n,n);
maxiter = 99;

if normest(X,1e-3) < sqrt(3)
    Xscale = X;
else
    Xscale = X./norm(X,'fro');
end

% Apply the Newton Schulz iteration
for i = 1:maxiter 

    Xold = Xscale;
    Xscale = 0.5*Xscale*(3*I - Xscale'*Xscale);
    
    diff_F = norm(Xscale - Xold,'fro');
    if (diff_F <= cgtol), break, end

end

U = Xscale;

if nargout > 1
    H = U'*X;
    H = (H+H')/2;
    k = i;
end

% Print convergence if required
if i == maxiter 
    fprintf("NS iteration doesn't converge after %2.0f iterations.\n", ...
        maxiter);
elseif (i ~= maxiter) && (nargin > 1)
    fprintf("NS iteration converged after %2.0f iteration.\n", i);
end

end