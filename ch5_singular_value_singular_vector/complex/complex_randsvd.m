function [A, sigma, metadata] = complex_randsvd(m, n, kappa, mode, precision, seed)
%COMPLEX_RANDSVD - Generate a complex matrix with prescribed singular values
%
%   Usage:
%       A = complex_randsvd(m, n, kappa, mode)
%       A = complex_randsvd(m, n, kappa, mode, precision)
%       [A,sigma,metadata] = complex_randsvd(m, n, kappa, mode, precision, seed)
%
%   Purpose:
%       COMPLEX_RANDSVD generates an mxn complex matrix
%           A = QL*diag(sigma)*QR'
%       with prescribed singular values, where QL and QR are obtained from
%       independent complex Gaussian matrices by phase-normalised QR
%       factorizations. The singular values follow the five positive modes
%       of gallery('randsvd') used by the real-arithmetic experiments.
%
%       MATLAB's default real gallery('randsvd') uses qmult METHOD 0 to
%       generate Haar orthogonal factors rather than QR. In exact arithmetic,
%       the construction here gives the corresponding Haar unitary/Stiefel
%       distributions, but it does not reproduce MATLAB's factor-generation
%       path or the same sample for a given seed.
%
%   Arguments:
%       (1) m - Positive integer, the number of rows of A.
%       (2) n - Positive integer, the number of columns of A, with n <= m.
%       (3) kappa - Real scalar greater than or equal to 1.
%       (4) mode - Integer from 1 to 5 specifying the singular values:
%           1: one large singular value;
%           2: one small singular value;
%           3: geometrically distributed singular values;
%           4: arithmetically distributed singular values;
%           5: random singular values with uniformly distributed logarithm.
%           For mode 5, cond(A) <= kappa rather than cond(A) = kappa.
%       (5) precision - 'single' or 'double', and by default 'double'.
%       (6) seed - Nonnegative integer, and by default 1.
%
%   Outputs:
%       (1) A - Complex mxn matrix with the prescribed singular values.
%       (2) sigma - Real nx1 vector containing the singular values of A.
%       (3) metadata - Structure recording the dimensions, requested and
%           actual condition numbers, mode, precision, seed, and generator.
%

if nargin < 5 || isempty(precision)
    precision = 'double';
end
if nargin < 6 || isempty(seed)
    seed = 1;
end

validateattributes(m, {'numeric'}, {'scalar','integer','positive'});
validateattributes(n, {'numeric'}, {'scalar','integer','positive','<=',m});
validateattributes(kappa, {'numeric'}, {'scalar','real','finite','>=',1});
validateattributes(mode, {'numeric'}, {'scalar','integer','>=',1,'<=',5});
validateattributes(seed, {'numeric'}, {'scalar','integer','nonnegative'});
precision = validatestring(precision, {'single','double'});

stream = RandStream('mt19937ar', 'Seed', double(seed));

XL = complex(cast(randn(stream,m,n),precision), ...
             cast(randn(stream,m,n),precision));
XR = complex(cast(randn(stream,n,n),precision), ...
             cast(randn(stream,n,n),precision));

[QL,RL] = qr(XL,0);
[QR,RR] = qr(XR,0);
QL = normalise_qr_phase(QL,RL);
QR = normalise_qr_phase(QR,RR);

switch mode
    case 1
        sigma = [1; repmat(1/kappa,n-1,1)];
    case 2
        sigma = [ones(n-1,1); 1/kappa];
    case 3
        if n == 1
            sigma = 1;
        else
            sigma = kappa.^(-(0:n-1)'/(n-1));
        end
    case 4
        sigma = linspace(1,1/kappa,n)';
    case 5
        if n == 1
            sigma = 1;
        else
            sigma = exp(-rand(stream,n,1)*log(kappa));
        end
end

sigma = cast(sigma,precision);
A = (QL .* reshape(sigma,1,[]))*QR';

metadata = struct('m',m,'n',n,'kappa',kappa,'mode',mode, ...
    'actual_kappa',double(max(sigma)/min(sigma)), ...
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
