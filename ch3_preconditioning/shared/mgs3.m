function [Qt, R] = mgs3( X,p )
%MGS3 - Block version Modified Gram Schmidt Algorithm for QR factorization.

% p : Block sizes. 

[m,n] = size( X );
s = n / p; % number of blocks. 
R = zeros( n,n ); 
Qt = zeros( m,n );
T = zeros( n,n ); 

[Qt1, R11, T1] = mgs2( X( :,1:p ) );
Qt( :,1:p ) = Qt1;
R( 1:p,1:p ) = R11;
T( 1:p,1:p ) = T1; 

for k = 2:s
    Hk = T( 1:(k-1)*p,1:(k-1)*p )'*( Qt( :,1:(k-1)*p )'*X( :,(k-1)*p+1:k*p ) );
    Yk = X( :,(k-1)*p+1:k*p ) - Qt( :,1:(k-1)*p )*Hk;
    [Qk, Rk, Tk] = mgs2(Yk);
    Fk = Qt( :,1:(k-1)*p )'*Qk;
    Gk = -T( 1:(k-1)*p,1:(k-1)*p )*Fk*Tk;
    Qt( :,(k-1)*p+1:k*p ) = Qk;
    R( 1:(k-1)*p,(k-1)*p+1:k*p ) = Hk;
    R( (k-1)*p+1:k*p,(k-1)*p+1:k*p ) = Rk;
    T( 1:(k-1)*p,(k-1)*p+1:k*p ) = Gk;
    T( (k-1)*p+1:k*p,(k-1)*p+1:k*p ) = Tk;
end 

end