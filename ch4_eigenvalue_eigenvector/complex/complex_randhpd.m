function [A,lambda,metadata] = complex_randhpd(n,kappa,mode,precision,seed)
%COMPLEX_RANDHPD - Complex Hermitian counterpart of complex_randsvd
%
%   Usage:
%       [A,lambda,metadata] = complex_randhpd(n,kappa,mode,'double',seed)
%   Purpose:
%       Use the one-sided generator's spectra and phase-normalised Gaussian
%       QR, replacing QL*diag(sigma)*QR' by Q*diag(lambda)*Q'. The rounded
%       matrix is Hermitian; its eigenvalues may differ from lambda.
%   Arguments:
%       n: order; kappa >= 1; mode: 1 large, 2 small, 3 geometric,
%       4 arithmetic, 5 random logarithmic eigenvalues; precision defaults
%       to 'double'; seed defaults to 1. Mode 5 need not attain kappa.
%   Outputs:
%       A: stored matrix; lambda: prescribed spectrum; metadata: settings
%       and prescribed_kappa. Check definiteness using the stored A.
%
if nargin < 4 || isempty(precision)
    precision = 'double';
end
if nargin < 5 || isempty(seed)
    seed = 1;
end

validateattributes(n, {'numeric'}, {'scalar','integer','positive'});
validateattributes(kappa, {'numeric'}, {'scalar','real','finite','>=',1});
validateattributes(mode, {'numeric'}, {'scalar','integer','>=',1,'<=',5});
validateattributes(seed, {'numeric'}, {'scalar','integer','nonnegative'});
precision = validatestring(precision, {'single','double'});

stream = RandStream('mt19937ar', 'Seed', double(seed));

X = complex(cast(randn(stream,n,n),precision), ...
            cast(randn(stream,n,n),precision));
[Q,R] = qr(X,0);
Q = normalise_qr_phase(Q,R);

switch mode
    case 1
        lambda = [1; repmat(1/kappa,n-1,1)];
    case 2
        lambda = [ones(n-1,1); 1/kappa];
    case 3
        if n == 1
            lambda = 1;
        else
            lambda = kappa.^(-(0:n-1)'/(n-1));
        end
    case 4
        lambda = linspace(1,1/kappa,n)';
    case 5
        if n == 1
            lambda = 1;
        else
            lambda = exp(-rand(stream,n,1)*log(kappa));
        end
end

lambda = cast(lambda,precision);
A = (Q .* reshape(lambda,1,[]))*Q';
A = (A + A')/2; % Exact Hermitian storage; reference the rounded A.

metadata = struct('n',n,'kappa',kappa,'mode',mode, ...
    'prescribed_kappa',double(max(lambda)/min(lambda)), ...
    'precision',precision,'seed',seed,'generator','complex Gaussian QR');
end

function Q = normalise_qr_phase(Q,R)
%NORMALISE_QR_PHASE - Fix the column phases implied by the diagonal of R.
d = diag(R);
phase = ones(size(d),'like',d);
nonzero = abs(d) > 0;
phase(nonzero) = d(nonzero)./abs(d(nonzero));
Q = Q .* reshape(phase,1,[]);
end
