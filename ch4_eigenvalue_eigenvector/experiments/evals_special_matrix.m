%TEST7 - Test for special matrices. (FIG 4) 

addpath("../shared/");
%% (a)
clear; close all; rng(1);
fprintf("[START] evals_special_matrix\n");
n = 100;
A = gallery('randsvd', n, -1e8, 3);

[~,D_ref] = eig(mp(A,1e3));
D_ref = diag(D_ref);

[~,D_eig] = eig(A);
[~,D_MP3,~,~,~,SCOND] = mp_pjacobi(A,'mp3');
BOUND = 7*n*(eps('double')/2)*SCOND;
[~,D_MP2] = mp_pjacobi(A,'mp2');
[~,D_Jacobi] = cjacobi(A);

error_eig = compute_error_evals(D_ref, diag(D_eig));
error_mp3 = compute_error_evals(D_ref, D_MP3);
error_mp2 = compute_error_evals(D_ref, D_MP2);
error_jacobi = compute_error_evals(D_ref, D_Jacobi);

path = sprintf("../output/random_matrix.csv");
Output = [(1:n)',flip(error_eig),flip(error_mp3),flip(error_mp2),...
    flip(error_jacobi),BOUND*ones(n,1)];
mywritetable(Output,{'k','matlab','mp3','mp2','jcb','bound'}, path);

disp("Finish (a)")

%% (b) 
clear; close all;
A = anymatrix('nessie/whiskycorr');
n = size(A,1);

[~,D_ref] = eig(mp(A,1e3));
D_ref = diag(D_ref);

[~,D_eig] = eig(A);
[~,D_MP3,~,~,~,SCOND] = mp_pjacobi(A,'mp3');
BOUND = 7*n*(eps('double')/2)*SCOND;
[~,D_MP2] = mp_pjacobi(A,'mp2');
[~,D_Jacobi] = cjacobi(A);

error_eig = compute_error_evals(D_ref, diag(D_eig));
error_mp3 = compute_error_evals(D_ref, D_MP3);
error_mp2 = compute_error_evals(D_ref, D_MP2);
error_jacobi = compute_error_evals(D_ref, D_Jacobi);

path = sprintf("../output/whiskycorr.csv");
Output = [(1:n)',flip(error_eig),flip(error_mp3),flip(error_mp2),...
    flip(error_jacobi),BOUND*ones(n,1)];
mywritetable(Output,{'k','matlab','mp3','mp2','jcb','bound'}, path);

disp("Finish (b)")

%% (c) 
clear; close all;
C = readmatrix("../data/superconductivty+data/train.csv");
A = cov(C);
n = size(A,1);
[~,D_ref] = eig(mp(A,1e3));
D_ref = diag(D_ref);

[~,D_eig] = eig(A);
[~,D_MP3,~,~,~,SCOND] = mp_pjacobi(A,'mp3');
BOUND = 7*n*(eps('double')/2)*SCOND;
[~,D_MP2] = mp_pjacobi(A,'mp2');
[~,D_Jacobi] = cjacobi(A);

error_eig = compute_error_evals(D_ref, diag(D_eig));
error_mp3 = compute_error_evals(D_ref, D_MP3);
error_mp2 = compute_error_evals(D_ref, D_MP2);
error_jacobi = compute_error_evals(D_ref, D_Jacobi);

path = sprintf("../output/covariance.csv");
Output = [(1:n)',flip(error_eig),flip(error_mp3),flip(error_mp2),flip(error_jacobi),BOUND*ones(n,1)];
mywritetable(Output,{'k','matlab','mp3','mp2','jcb','bound'}, path);

disp("Finish (c)")

%% (d) 
close all; clear;
D = gallery('lauchli',500,1e-3);
A = D'*D;
n = size(A,1); 

% Run experiments 
[~,D_ref] = eig(mp(A,1e2));
D_ref = diag(D_ref);

% D_eig = svd(D);
D_eig = eig(A);
[~,D_MP3,~,~,~,SCOND] = mp_pjacobi(A,'mp3');
BOUND = 7*n*(eps('double')/2)*SCOND;
[~,D_MP2] = mp_pjacobi(A,'mp2');
[~,D_Jacobi] = cjacobi(A);


error_eig = compute_error_evals(D_ref, D_eig);
error_mp3 = compute_error_evals(D_ref, D_MP3);
error_mp2 = compute_error_evals(D_ref, D_MP2);
error_jacobi = compute_error_evals(D_ref, D_Jacobi);

selected=1:10:n;
errors_desc = [flip(error_eig),flip(error_mp3),...
    flip(error_mp2),flip(error_jacobi)];

path = sprintf("../output/lau.csv");
Output = [selected',errors_desc(selected,:),...
    ones(length(selected),1) * BOUND];
mywritetable(Output,{'k','matlab','mp3','mp2','jcb','bound'}, path);

disp("Finish (d)")

fprintf("[END] evals_special_matrix\n");
