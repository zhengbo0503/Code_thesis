function r = scaled_cond(A)
%CONDS - Diagonally scaled 2 norm condition number
%	r = conds(A) computes the condition number of diagonally scaled A, here
%	A is scaled using a diagonal matrix in two-sided fashion. Moreover, the
%	scaled A should have diagonals all ones. 
%
%   See also cond, condest
%

% Check input
if ~(issymmetric(A)), error("Input must be symmetric."); end 

% Main
D = diag(diag(A).^(-1/2));
r = cond(D*A*D);

end

