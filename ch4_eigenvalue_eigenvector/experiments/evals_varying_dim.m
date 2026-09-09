%TEST2 - Test script for varying matrix dimension (FIG 1) 
%   The working precision is IEEE double, low precision is IEEE single, and
%   high precision is IEEE quadruple. 
% 
%   For each eigenvalue mode, and a fixed condition number. 
%   We create 20 matrices with varying dimensions froom 1 to 2500. 
%   For each one of them, we compute the eigenvalues and the corresponding
%   forward error. 
%

% Test initialization 
clear; close all; rng(1);
fprintf("[START] evals_varying_dim\n");
addpath("../shared/");

% Parameters 
u = float_params('d');
uh = float_params('q');
N = round(logspace(1,3,20));

% Stored data 
bound_store = zeros(length(N),1);
maxfwd_mp2 = zeros(length(N),1);
maxfwd_mp3 = zeros(length(N),1);
maxfwd_jacobi = zeros(length(N),1);
maxfwd_eig = zeros(length(N),1);
MP2 = zeros(length(N), 5);
MP3 = zeros(length(N), 5);
JACOBI = zeros(length(N), 5);
EIG = zeros(length(N), 5);
BOUNDI = zeros(length(N), 5);

% Varying variable
kappa = 1e10;

for mode = 1:5
    % Main 
    for i = 1:length(N) % loop over matrix dimension n        
        n = N(i);

	    % Generate test matrix
	    A = gallery('randsvd',n,-kappa,mode);

        % Compute the preconditioner
        [Ql,~] = eig(single(A));
        [Qt,~] = qr(double(Ql));

        % Compute reference eigenvalues at highest precision
        [~,D_mp] = eig(mp(A,71));
	    D_ref = sort(diag(D_mp));

        %  Apply eigensolvers
	    [~,D_eig] = eig(A);
        [~,D_MP3,~,~,BOUND] = mp_pjacobi(A,'mp3');
        [~,D_MP2] = mp_pjacobi(A,'mp2');
        [~,D_Jacobi] = cjacobi(A);
        
        error_eig = compute_error_evals(D_ref, diag(D_eig), 144);
        error_mp3 = compute_error_evals(D_ref, D_MP3, 144);
        error_mp2 = compute_error_evals(D_ref, D_MP2, 144);
        error_jacobi = compute_error_evals(D_ref, D_Jacobi, 144);

        maxfwd_mp2(i) = max(error_mp2);
        maxfwd_mp3(i) = max(error_mp3);
        maxfwd_jacobi(i) = max(error_jacobi);
        maxfwd_eig(i) = max(error_eig);

        % Compute the "exact" At
        At = mp(Qt,71)' * mp(A,71) * mp(Qt,71);
        At = (At + At')/2;

        bound_store(i) = 7*n*u*scond(At,'D');
        
        fprintf("Mode=%d, Finished %d of %d\n",mode, i, length(bound_store)); 
    end
    
    % Store
    MP2(:, mode) = maxfwd_mp2;
    MP3(:, mode) = maxfwd_mp3;
    JACOBI(:, mode) = maxfwd_jacobi;
    EIG(:, mode) = maxfwd_eig;
    BOUNDI(:, mode) = bound_store;
end

for mode = 1:5
    i = mode; 
    path = sprintf("../output/varydim_mode%d.csv", mode);
    Output = [N',MP2(:,i),MP3(:,i),JACOBI(:,i),EIG(:,i),BOUNDI(:,i)];
    mywritetable(Output,{'kappa','precj','uhprecj','jacobi','matlab','bound'}, path);
end

fprintf("[END] evals_varying_dim\n");
