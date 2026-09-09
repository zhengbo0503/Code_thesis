%TEST1 - Test script for varying condition number (FIG 2) 
%
%   The working precision is IEEE double, low precision is IEEE single, and
%   high precision is IEEE quadruple. 
% 
%   For each eigenvalue mode, 
%       For each condition number, we create matrices and then compute the
%       eigenvalues using four different methods. 
%       Computing its error and the corresponding bound. 
%

% Test initialization
clear; close all; rng(1);
fprintf("[START] evals_varying_cond\n");
addpath("../shared/");

% Parameters 
n = 100;
N = 1;
u = float_params('d');
ul = float_params('s');
uh = float_params('q');
gammah = n * uh/(1 - n * uh);
counter = 0;
test_kappa = logspace(1,16,20)';

% Stored data 
bound_store = zeros(length(test_kappa),1);
maxfwd_mp2 = zeros(length(test_kappa),1);
maxfwd_mp3 = zeros(length(test_kappa),1);
maxfwd_jacobi = zeros(length(test_kappa),1);
maxfwd_eig = zeros(length(test_kappa),1);
MP2 = zeros(length(test_kappa), 5);
MP3 = zeros(length(test_kappa), 5);
JACOBI = zeros(length(test_kappa), 5);
EIG = zeros(length(test_kappa), 5);
BOUNDI = zeros(length(test_kappa), 5);

for mode = 1:5
    % Main
    for i = 1:length(test_kappa) % Loop over condition numbers
        kappa = test_kappa(i);

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
        
        % Compute errors 
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

        % Compute error bound
        bound_store(i) = 7*n*u*scond(At,'D');
        
        % Print info 
        fprintf("Mode=%d, Finished %d of %d\n",mode, i, length(bound_store)); 
    end

    % Store for overleaf
    path = sprintf("../output/varykappa_mode%d.csv", mode);
    Output = [test_kappa,maxfwd_mp2,maxfwd_mp3,maxfwd_jacobi,maxfwd_eig,bound_store];
    mywritetable(Output,{'kappa','precj','uhprecj','jacobi','matlab','bound'}, path);
end

fprintf("[END] evals_varying_cond\n");
