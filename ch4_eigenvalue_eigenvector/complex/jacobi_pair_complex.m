function [c,s] = jacobi_pair_complex(A,p,q)
%JACOBI_PAIR_COMPLEX - Unitary Jacobi rotation entries (c,s)
% Minimal Hermitian adaptation of ../shared/jacobi_pair.m.
% J = [c,s;-conj(s),c], with s carrying the phase of A(p,q).
if A(p,q) == 0
    c = 1; s = 0;
else
    rho = abs(A(p,q));
    tau = (real(A(q,q))-real(A(p,p)))/(2*rho);
    if tau >= 0
        t = 1/(tau + sqrt(1+tau*tau));
    else
        t = -1/(-tau + sqrt(1+tau*tau));
    end
    c = 1/sqrt(1+t*t);
    s = (A(p,q)/rho)*(t*c);
end
end
