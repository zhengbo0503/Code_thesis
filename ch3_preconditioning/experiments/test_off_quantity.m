%The size of off quantity (FIG 1)

% Initialization
clear; close all; rng(1);
fprintf("[START] test_off_quantity\n");
addpath("../shared/");

% Parameters
NoP = 20;
n = round(logspace(1,3,NoP));
ul = 2^(-24);
u = 2^(-53);
kappa = 1e6;

off_qr = zeros(NoP,1);
off_mgs = zeros(NoP,1);
off_ns = zeros(NoP,1);
off_tridiag = zeros(NoP,1);
off_A = zeros(NoP,5);

% Main
for mode = 1:5
    for i = 1:length(n)        
        A = gallery('randsvd',n(i),-kappa,mode);
        fnrm = norm(A,'fro');
        off_A(i,mode) = off(A);

        % Orthogonalization methods
        [Ql,~] = eig(single(A));
        Qd = double(Ql);
        
        [Qt_qr,~] = qr(Qd);    % HouseQR
        [Qt_mgs,~] = mgs(Qd);  % MGS
        Qt_ns = ns(Qd);        % NS

        % tridiagonalization methods
        Qt_tridiag = mpdiag(A);

        % collect errors
        off_qr(i) = off_error(A,fnrm,Qt_qr);
        off_mgs(i) = off_error(A,fnrm,Qt_mgs);
        off_ns(i) = off_error(A,fnrm,Qt_ns);
        off_tridiag(i) = off_error(A,fnrm,Qt_tridiag);

        fprintf("Mode = %d, iteration %d of %d\n", mode, i, length(n));
    end

    path = sprintf("../output/offAt_mode%d.csv", mode);
    Output = [n',off_qr,off_mgs,off_ns,off_tridiag,(5*ul*n.^(1/2))'];
    mywritetable(Output,{'n','qr','mgs','ns','tri','bound'}, path)
end

fprintf("[END] test_off_quantity\n");

function error = off_error(A,fnrm,Qt)
A = mp(A,34);
Qt = mp(Qt,34);
error = off(Qt'* A * Qt)/fnrm;
end

