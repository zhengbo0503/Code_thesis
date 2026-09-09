function err = compute_error_evals(D_ref, D_comp, varargin)
%COMPUTE_ERROR_EVALS - Computing error for our testing scripts 
%
%   Usage:
%       err = compute_error_evals(D_ref, D_comp)
%       err = compute_error_evals(D_ref, D_comp, Precision)
%
%   Purpose:
%       Compute relative forward error of two sets of eigenvalues at high
%       precision.
%
%   Input: 
%    - D_ref: (Any Precision) Vector, Dim(n,1)
%       This is the vector of reference eigenvalue that obtained either
%       from higher precision or some analytical results.
%    - D_comp: (Any Precision) Vector, Dim(n,1)
%       The computed eigenvalues.
%    - (Optional) Precision: (Any Precision) Integer, Dim(1,1) 
%       The precision that inputs into Advanpix Toolbox. Such precision is
%       used for computing the error. By default, Precision=71. 
%
%   Output: 
%    - err: (Precison predefined or user-defined) Vector, Dim(n,1)
%       The relative forward error computed via 
%           abs((D_ref - D_comp)./D_ref)
%
%   Author: 
%       Zhengbo Zhou, June 2025, Manchester, UK 
%

% Set predefined precision or user-defined precision
if nargin > 2
    prec = varargin{1};
else
    prec = 71;
end

% Sort the eigenvalues, promote to high precision
D_ref = mp(sort(D_ref), prec);
D_comp = mp(sort(D_comp), prec);

% Check the dimension to remove the possibility that one of them is a nx1
% vector and another one is 1xn vector. 
[n_ref,~] = size(D_ref);
[n_comp,~] = size(D_ref);
if (n_ref ~= n_comp)
    D_ref = D_ref';
end

% Computes the forward error 
err = abs((D_ref - D_comp)./D_ref);
end
