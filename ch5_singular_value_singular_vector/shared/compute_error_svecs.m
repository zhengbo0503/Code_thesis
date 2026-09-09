function [err,relgap] = compute_error_svecs(A, V, Vtype, Vref, S)
%COMPUTE_ERROR_SVECS compute the error of the computed singular/eigen-vector
%	Given a matrix A in mxn
%	Given singular/eigen-vectors V.
%	Given a char Vtype that is either
%		- 'S': right singular vectors
%		- 'L': left singular vectors
%		- 'E': eigenvectors
%	The function computes the error in terms  of
%		sin(angle( v(A,k), v(:,k) ))
%	where v(A,k) is the kth exact singular/eigen-vector of A.
%	Moreover, we compute the relative gap of the singular/eigen-value
%	which is defined as
%		svgap(A,k) = min{ min_{j=/=k} {|sigma(A,k)-sigma(A,j)|/sigma(A,k)},2}
%	for singular values, and for the eigenvalues, we have
%		evgap(A,k) = min_{j=/=k} {|lambda(A,k)-lambda(A,j)|/sqrt(lambda(A,k)*lambda(A,j))}
%
%   Optional arguments:
%       Vref, S — precomputed reference vectors and values from get_reference().
%                  When provided, skips the mp decomposition of A.

% Evaluate vector errors in the same 71-digit precision as the references.
V = mp(V, 71);

if nargin >= 5
    Vmp = Vref;
    n = length(S);
else
    [~,n] = size(A);
    Amp = mp(A, 71);
    switch Vtype
      case {'S', 's'}
        [~,Smp,Vmp] = svd(Amp,'econ');
        S = diag(Smp);
      case {'L', 'l'}
        [Ump,Smp,~] = svd(Amp,'econ');
        Vmp = Ump;
        S = diag(Smp);
      case {'E', 'e'}
        [Vmp,Dmp] = eigsort(Amp);
        S = diag(Dmp);
      otherwise
        error("Invalid Vtype as input!");
    end
end

err = zeros(1, n);
relgap = zeros(1, n);

switch Vtype
  case {'S', 's', 'L', 'l'}
    for i = 1:n
        x = Vmp(:,i);
        y = V(:,i);

        nx = norm(x);
        ny = norm(y);

        u = x / nx;
        v = y / ny;

        if u' * v < 0
            v = -v;
        end

        d = norm(u - v);
        err(i) = d * sqrt(max(0, 1 - d^2/4));
        gap_i = abs(S(i) - S) / S(i);
        gap_i(i) = inf;
        relgap(i) = double(min(min(gap_i), 2));
    end
  case {'E', 'e'}
    for i = 1:n
        x = Vmp(:,i);
        y = V(:,i);
        nx = norm(x);
        ny = norm(y);
        u = x / nx;
        v = y / ny;
        if u' * v < 0
            v = -v;
        end
        d = norm(u - v);
        err(i) = d * sqrt(max(0, 1 - d^2/4));
        gap_i = abs(S(i) - S) ./ sqrt(S .* S(i));
        gap_i(i) = inf;
        relgap(i) = min(gap_i);
    end
  otherwise
    error("Invalid Vtype as input!");
end
end
