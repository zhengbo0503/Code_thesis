function nrm = scond(A, type, nrmtype)
%SCOND - Scaled condition number
%	nrm = SCOND(A,type,nrmtype) computes the condition number of
%	a scaled A.
%	If A is symmetric, then it operates two-sided scaling such that
%	the diagonals of the scaled A are all one, and maintaining symmetry.
%	If A is not symmetric, then it operates column scaling, such that
%	the columns of the scaled A all have unit 2 norm.
%
%	Arguments:
%		(1) A - Double, Real matrix
%			symmetric or general matrix.
%		(2) type - {'diagonal','column'} or {'D','C'}.
%			By default, it is determined by whether A is
%			symmetric ('D') or not ('C').
%		(3) nrmtype - type of norm used to compute condition number
%			By default is 2. One can pass any type as long as
%			it is recognized by MATLAB cond.
%
%	Outputs:
%		(1) nrm - Double, Real scalar
%			the scaled condition number.
%
%	Author:
%		Zhengbo Zhou, Manchester, UK, 2025 Dec
%
%	Note (thesis merge):
%		This is the single canonical scond for the whole thesis. It
%		supersedes the earlier Code_htwz25 version whose SECOND argument
%		was the norm type rather than the scaling type. Every call site
%		has been made explicit about the scaling, so the two are
%		equivalent wherever scond was previously used.

% Check twosided or onesided
if nargin == 1	% type and nrmtype are not defined
    if (issymmetric(A))
        type = 'D';
        nrmtype = 2;
    else
        type = 'C';
        nrmtype = 2;
    end
elseif nargin == 2
    % nrmtype is not defined
    nrmtype = 2;
end

switch type
  case {'C', 'Column', 'column'}
    D1 = eye(size(A,1));
    col_norms = sqrt(sum(A.^2, 1));
    D2 = diag(1 ./ col_norms);
  case {'D', 'Diagonal', 'diagonal'}
    D1 = diag(diag(A).^(-1/2));
    D2 = D1;
  otherwise
    error("Invalid scaling type, either 'D' (diagonal scaling) or 'C' (column scaling).");
end

nrm = cond(D1*A*D2, nrmtype);

end
