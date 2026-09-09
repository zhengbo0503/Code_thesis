%TEST5 - mposj_ssd with (single, single, double)
addpath("../shared/");
%	compared with SGESVJ and SGEJSV

close all; clear; rng(1);
fprintf("[START] svals_ssd\n");

m = 3000*ones(10,1); 
n = round(logspace(2, log10(3000), 10));
epsln = eps('single')/2;
time_mposj = zeros(length(m),1);
time_sgesvj = zeros(length(m),1);
time_sgejsv = zeros(length(m),1);
time_matlab = zeros(length(m),1);
warmup_matrix = gallery('randsvd', [10, 5], 1e6, 3, 'single');

for i = 1:length(m)

    mm = m(i);
    nn = n(i);
    A = gallery('randsvd', [mm,nn], 1e6, 3, 'single');
    sref = reference_singular_values(A);

    % Our algorithm 
    mposj_ssd(warmup_matrix); % preload function 
    tic_mposj_tmp1 = tic;
    [~,~,~,~] = mposj_ssd(A);
    tic_mposj_tmp1 = toc(tic_mposj_tmp1);
    tic_mposj_tmp2 = tic;
    [U1,S1,V1,nos] = mposj_ssd(A);
    tic_mposj_tmp2 = toc(tic_mposj_tmp2);
    time_mposj(i) = (tic_mposj_tmp1 + tic_mposj_tmp2)/2;
    [f1(i),~,~,~] = compute_error_svals(A, U1, S1, V1, sref);

    % Scaled condition number for the reference line, computed outside
    % the timed runs so that it does not inflate time_mposj.
    [~,~,~,~,scnd] = mposj_ssd(A, true);
    
    % SGESVJ (plain Jacobi)
    sgesvj_mex(warmup_matrix, 'G', 'U', 'V', size(warmup_matrix, 2), eye(size(warmup_matrix, 2),'single'), size(warmup_matrix,1)+size(warmup_matrix,2)); 
    tic_sgesvj_tmp1 = tic;
    [~,~,~,~,~,~] = sgesvj_mex(A,'G','U','V',nn,eye(nn,'single'),max(6,mm+nn));
    tic_sgesvj_tmp1 = toc(tic_sgesvj_tmp1);
    tic_sgesvj_tmp2 = tic;
    [U2,S2,V2,sva2,work2,info2] = sgesvj_mex(A,'G','U','V',nn,eye(nn,'single'),max(6,mm+nn));
    tic_sgesvj_tmp2 = toc(tic_sgesvj_tmp2);
    if info2 ~= 0
        fprintf("Error: SGESVJ does not converge.\n");
        break;
    end
    time_sgesvj(i) = (tic_sgesvj_tmp1 + tic_sgesvj_tmp2)/2;
    [f2(i),~,~,~] = compute_error_svals(A, U2, S2, V2, sref);

    % DGEJSV (preconditioned Jacobi)
    sgejsv_mex(warmup_matrix,'C','U','V','R','N','N');
    tic_sgejsv_tmp1 = tic;
    sgejsv_mex(A,'C','U','V','R','N','N');
    tic_sgejsv_tmp1 = toc(tic_sgejsv_tmp1);
    tic_sgejsv_tmp2 = tic;
    [U3,S3,V3,sva3,work3,iwork3,info3] = sgejsv_mex(A,'C','U','V','R','N','N');
    tic_sgejsv_tmp2 = toc(tic_sgejsv_tmp2);
    if info3 ~= 0
        fprintf("Error: SGEJSV does not converge.\n");
        break;
    end
    time_sgejsv(i) = (tic_sgejsv_tmp1 + tic_sgejsv_tmp2)/2;
    [f3(i),~,~,~] = compute_error_svals(A, U3, S3, V3, sref);


    bound2(i) = scnd * sqrt(mm * nn) * epsln;

    fprintf("Finished %d of %d \n", i, length(n));

end

savedata = 1;
if savedata == 1
    save("../data/svals_ssd.mat")
end

%% Export the plotted data for the thesis figure
% mp3_* names the thesis MPJacobiSVD_SSD curve produced by mposj_ssd; this
% experiment is single/single/double and does not use quadruple precision.
csv_header = {'n','mp3_err','sgesvj_err','sgejsv_err','bound', ...
    'mp3_s','sgesvj_s','sgejsv_s'};
csv_values = [n(:), f1(:), f2(:), f3(:), bound2(:), ...
    time_mposj(:), time_sgesvj(:), time_sgejsv(:)];
write_figure_csv("../data/svals_ssd.csv", csv_values, csv_header);
%% 
close all;
C1 = "#1171BE";
C2 = "#DD5400";
C3 = "#EDB120";
C4 = "#3BAA32";

figure(1)

loglog(n,f1,'LineStyle','none','Marker','*','Color',C1);
hold on;
loglog(n,f2,'LineStyle','none','Marker','pentagram','Color',C2);
loglog(n,f3,'LineStyle','none','Marker','square','Color',C3);
loglog(n,bound2,'LineStyle',':','Marker','none','Color','k');

legend('MP3JacobiSVD', 'SGESVJ', 'SGEJSV', 'Bound');
set(findall(gcf, 'Type', 'Line'), 'LineWidth', 1);
xlabel('$n$', 'FontSize', 10);
ylabel('$\mathrm{max}_k {\varepsilon}^{(k)}_{fwd}$', 'FontSize', 10, 'Interpreter', 'latex');
ylim([1e-7, 10]);

figure(2)
loglog(n, time_mposj,'LineStyle','none','Marker','*','Color',C1); hold on;
loglog(n, time_sgesvj,'LineStyle','none','Marker','pentagram','Color',C2);
loglog(n, time_sgejsv,'LineStyle','none','Marker','square','Color',C3);
legend('MP3JacobiSVD', 'SGESVJ', 'SGEJSV')

set(findall(gcf, 'Type', 'Line'), 'LineWidth', 1);
xlabel('Number of columns', 'FontSize', 10); 
ylabel('Runtime (sec)', 'FontSize', 10); 

fprintf("[END] svals_ssd\n");
