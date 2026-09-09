%Test on reduction of scaled condition number. (FIG 5)
clear; rng(1);
fprintf("[START] test_scond_reduction\n");
addpath("../shared/");

kappa_test = logspace(1,12,20)';

for mode = 1:5

    % Initialize spaces 
    kappa_original = ones(length(kappa_test),1);
    r_DAD = ones(length(kappa_test),1);
    r_DAtD = ones(length(kappa_test),1);
    r_DAtcompD = ones(length(kappa_test),1);
    
    for i = 1:length(kappa_test)

        A = gallery('randsvd', 100, -kappa_test(i), mode);

        kappa_original(i) = cond(mp(A,71));

        % Atcomp
        [Ql,~] = eig(single(A));
	    [Qt,~] = qr(double(Ql));
        At_hcomp = mp(Qt)'* mp(A) * mp(Qt);
	    At_hcomp = (At_hcomp + At_hcomp')/2;
	    Atcomp = double(At_hcomp);

        % Exact At
        At = mp(Qt,71)' * mp(A,71) * mp(Qt,71);
        At = (At + At')/2;

        r_DAD(i) = scond(mp(A,71),'D')/kappa_original(i);
        r_DAtD(i) = scond(mp(At,71),'D')/kappa_original(i);
        r_DAtcompD(i) = scond(mp(Atcomp,71),'D')/kappa_original(i);

        fprintf("MODE=%d, iteration %d of %d\n", mode, i, length(kappa_test));
        
    end

    path = sprintf("../output/reducekappa_mode%d.csv", mode);
    Output = [kappa_original, r_DAD, r_DAtD, r_DAtcompD,...
              ones(length(kappa_test),1)];
    mywritetable(Output,{'ka','DAD','DAtD','DAtcompD','one'}, path);
end

fprintf("[END] test_scond_reduction\n");
