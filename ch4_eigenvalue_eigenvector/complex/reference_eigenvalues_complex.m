function dref = reference_eigenvalues_complex(A)
%REFERENCE_EIGENVALUES_COMPLEX - Reference the stored Hermitian matrix at 71 digits
%   Usage: dref = reference_eigenvalues_complex(A)
%   Purpose: Compute eigenvalues of the actual stored input, not its
%       prescribed spectrum. Reject loss of positive definiteness.
%   Arguments: A is the stored Hermitian double matrix.
%   Outputs: dref is a descending 71-digit multiprecision column vector.
assert(ishermitian(A),'The reference matrix must be Hermitian.');
dref = eig(mp(A,71));
assert(isreal(dref),'The Hermitian reference returned nonreal eigenvalues.');
dref = sort(dref(:),'descend');
assert(all(dref > 0),'The stored matrix is not positive definite.');
end
