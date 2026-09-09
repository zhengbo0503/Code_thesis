function [f, r, oU, oV] = compute_error_svals(A, U, S, V, sref)
%COMPUTE_ERROR_SVALS - Computing error for our testing scripts 
%
%   Usage:
%       [f, r, oU, oV] = COMPUTE_ERROR(A, U, S, V)
%       [f, r, oU, oV] = COMPUTE_ERROR(A, U, S, V, sref)
%
%   Purpose:
%       Compute several errors for the computed singular value
%       decomposition. If sref is not supplied, the reference singular
%       values are computed at 71 decimal digits. The function then
%       computes the
%		following errors:
%		1. Relative forward error
%			max_k {|sigma_k - sigma_ref_k|/sigma_ref_k}
%		2. Backward error:
%			||A-U*S*V'||_F
%		3. Orthogonality errors
%			||U'*U - I||_inf and ||V'*V - I||_inf
%
%   Arguments: 
%		(1) A - Real matrix
%	       The input matrix to the SVD algorithm.
%		
%       (2,3,4) U, S, V - Real matrices
%   	    The computed singular value decomposition
%       (5) sref - Optional multiprecision column vector
%           Reference singular values for A in descending order
%
%   Output: 
%		(1) f - Real, double scalar
%			Maximum relative forward error
%		(2) r - Real, double scalar
%			Backward error
%	    (3) oU - Real, double scalar
%			Orthogonality error of the left singular vectors
%	    (4) oV - Real, double scalar
%			Orthogonality error of the right singular vectors
%		
%   Author: 
%       Zhengbo Zhou, Manchester, UK, Dec 2025
%

r = norm(A - U*S*V','fro')/norm(A,'fro');
oU = norm(U'*U - eye(size(U,2)), inf);
oV = norm(V'*V - eye(size(V,2)), inf);

% Forward error needs special treatment
if nargin < 5 || isempty(sref)
    sref = reference_singular_values(A);
else
    sref = sref(:);
end

if numel(sref) ~= min(size(A))
    error('compute_error:InvalidReferenceSize', ...
        'sref must contain min(size(A)) singular values.');
end

s = sort(diag(S), 'descend');
s = s(:);

f = double(max(abs(sref - s)./abs(sref)));
end
