function [Qt,R,T] = mgs2(X)
%MGS2 - Matrix-vector based modified Gram Schmidt algorithm
%
%   [Qt,R,T] = MGS2(X) computes the QR factorization of X, which satisfies
%   X = Qt*R. And T satisfies (2.10)-(2.11) in [1].
%
%   Input arguments:
%   X : Full column rank matrix. Hence, we require the column number is
%       greater than the row number. Therefore, if we see the reverse
%       happened, we will apply the MGS algorithm to X'. 
%
%   Output arguments:
%   Qt : Left orthogonal, i.e. norm( Qt'*Qt - I ) = O( u ).
%   R  : Upper triangular matrix. 
%   T  : Upper triangular matrix.
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
    X = X'; temp = m; m = n; n = temp; 
    warning("Input matrix have full row rank. Apply MGS to X'.");
end

% Allocate spaces
R = zeros( n,n ); 
Qt = zeros( m,n );
T = zeros( n,n ); 

% Initialize outputs 
R( 1,1 ) = norm( X( :,1 ) );
Qt( :,1 ) = X( :,1 ) / R( 1,1 ); 
T( 1,1 ) = 1; 

% Main 
for k = 2:n
    hk = T( 1:k-1,1:k-1 )'*Qt( :,1:k-1 )'*X( :,k ); 
    yk = X(:,k) - Qt( :,1:k-1 )*hk; 
    R( k,k ) = norm( yk ); 
    Qt( :,k ) = yk / R( k,k ); 
    gk = -1 * T( 1:k-1,1:k-1 )*Qt( :,1:k-1 )'*Qt( :,k ); 
    R( 1:k-1,k ) = hk;
    T( 1:k-1,k ) = gk;
    T( k,k ) = 1; 
end

end 