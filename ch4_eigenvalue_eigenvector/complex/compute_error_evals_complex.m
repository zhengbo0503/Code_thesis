function [f,r,oV] = compute_error_evals_complex(A,V,D,dref)
%COMPUTE_ERROR_EVALS_COMPLEX - Eigenvalue counterpart of compute_error_svals
%   Usage: [f,r,oV] = compute_error_evals_complex(A,V,D,dref)
%   Purpose: Compare eigenvalues to the stored-input multiprecision reference.
%   Arguments: A = V*diag(D)*V' approximately; dref is the sorted reference.
%   Outputs: Maximum relative eigenvalue error, relative Frobenius residual,
%       and norm(V'*V-I,inf). Forward-error arithmetic uses dref precision.
r = norm(A - V*diag(D)*V','fro')/norm(A,'fro');
oV = norm(V'*V - eye(size(V,2)),inf);
if nargin < 4 || isempty(dref)
    dref = reference_eigenvalues_complex(A);
end
assert(numel(dref) == size(A,1),'Reference has the wrong number of eigenvalues.');
d = sort(D(:),'descend');
f = double(max(abs(dref(:) - d)./abs(dref(:))));
end
