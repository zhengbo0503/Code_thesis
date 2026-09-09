% check_fresh_mex - verify a freshly built dgesvj_mex in the current folder
rng(0); m = 50; n = 35;
A = gallery('randsvd',[m,n],1e6,3);
[U,S,V,sva,work,info] = dgesvj_mex(A,'G','U','V',n,eye(n),max(6,m+n));
bwd  = norm(A - U*S*V','fro')/norm(A,'fro');
orth = norm(U'*U - eye(n),'fro');
fprintf('info=%d backward=%.3e orth=%.3e\n', info, bwd, orth);
if info == 0 && bwd < 1e-13 && orth < 1e-12
    disp('FRESH_MEX_OK');
else
    disp('FRESH_MEX_BAD');
end
