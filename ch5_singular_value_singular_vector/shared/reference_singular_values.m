function sref = reference_singular_values(A)
%REFERENCE_SINGULAR_VALUES Compute reference singular values at 71 digits.
%
%   sref = REFERENCE_SINGULAR_VALUES(A) converts the stored matrix A to
%   71-decimal-digit multiprecision, computes its singular values, and
%   returns them as a descending column vector.

reference_digits = 71;
sref = svd(mp(A, reference_digits), 'econ');
sref = sort(sref(:), 'descend');
end
