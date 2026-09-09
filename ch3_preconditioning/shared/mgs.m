function [Q,R] = mgs(X)
%MGS - Modified Gram Schmidt algorithm. 
%
%   Usage: 
%       [Q,R] = MGS(X)     
%    
%   Purpose: 
%       computes the QR factorization of X using a matrix-vector 
%       multiplication rich algorithm, such that 
%              || QR - X || = O( u )   --> Backward stability
%             || Q'*Q - I || = O( u )  --> Left orthogonality
%       If X is divisable by some number that is larger than 32, 
%       then we use the Block-MGS algorithm. 
%       The block size is determined by the smallest divisor of the no. of
%       col of X that is greater than 32.
%       E.g. For X with dimension 80 x 80, the block size will be 40.
%
%       "32" is chosen based on a M3 Pro Chip with 32GiB of RAM. 
%
%   Input: 
%    - X is REAL matrix, dimension (M,N) 
%       Real matrix with full column rank. If M <= N, then we will perform
%       mgs(X'). 
%   
%   Output:
%    - Qt is REAL matrix, dimension (M,N)
%       Left orthogonal, i.e. Qt'*Qt = I.
%    - R is REAL matrix, dimension (N,N) 
%       R is upper triangular. 
%
%   Author:
%       Zhengbo Zhou, Nov 2024, Manchester
%
%   Reference: 
%   [1] Barlow, Block Modified Gram--Schmidt Algorithms and Their Analysis, 
%       SIAM J. Matrix Anal. Appl. 40, 1257–1290, 2019
   

% Set parameters 
[m,n] = size( X );
rk = rank( X ); 

% Check full rank
if ( rk ~= n ) && ( rk ~= m )
    warning("Input matrix should have full column rank.");
elseif ( rk ~= n ) && ( rk == m )
    X = X'; n = m;
    warning("Input matrix have full row rank. Apply MGS to X'.");
end

% Determine block size
factors = factor( n ); 
blockSize = 1;
for i = 1:length( factors )
    blockSize = blockSize * factors( i );
    if ( blockSize >= 32 ), break; end
end

% Compute QR factorization 
if blockSize ~= n 
    [Q,R] = mgs3( X,blockSize ); 
else
    [Q,R,~] = mgs2( X );
end

end