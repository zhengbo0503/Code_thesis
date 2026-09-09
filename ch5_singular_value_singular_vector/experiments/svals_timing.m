%TEST4 -- timing test, changing number of cols
addpath("../shared/");

close all; clear; rng(1);
fprintf("[START] svals_timing\n");

m = 3000*ones(10,1); 
n = round(logspace(2, log10(3000), 10));
time_mposj = zeros(length(m),1);
time_dgesvj = zeros(length(m),1);
time_dgejsv = zeros(length(m),1);
time_matlab = zeros(length(m),1);

for i = 1:length(m)

    mm = m(i);
    nn = n(i);
    A = gallery('randsvd', [mm,nn], 1e8, 3);

    % Our algorithm 
    tic_mposj_tmp1 = tic;
    [~,~,~,~,~,timing_tmp1] = mposj(A);
    tic_mposj_tmp1 = toc(tic_mposj_tmp1);
    tic_mposj_tmp2 = tic;
    [U,S,V,nos,scalecond,timing_tmp2] = mposj(A);
    tic_mposj_tmp2 = toc(tic_mposj_tmp2);
    mposj_disect_timing(i,:) = (timing_tmp1 + timing_tmp2)/2;
    time_mposj(i) = (tic_mposj_tmp1 + tic_mposj_tmp2)/2;

    % Use the two precision version to get a standard for double precision
    % preconditioning process.
    [~,~,~,~,~,timing_twoprec_tmp1] = mposj(A,2);
    [~,~,~,~,~,timing_twoprec_tmp2] = mposj(A,2);
    time_mposj_twoprec(i,:) = (timing_twoprec_tmp1 + timing_twoprec_tmp2)/2;
    
    % DGESVJ (plain Jacobi)
    tic_dgesvj_tmp1 = tic;
    dgesvj_mex(A,'G','U','V',nn,eye(nn),max(6,mm+nn));
    tic_dgesvj_tmp1 = toc(tic_dgesvj_tmp1);
    tic_dgesvj_tmp2 = tic;
    dgesvj_mex(A,'G','U','V',nn,eye(nn),max(6,mm+nn));
    tic_dgesvj_tmp2 = toc(tic_dgesvj_tmp2);
    time_dgesvj(i) = (tic_dgesvj_tmp1 + tic_dgesvj_tmp2)/2;

    % DGEJSV (preconditioned Jacobi)
    tic_dgejsv_tmp1 = tic;
    dgejsv_mex(A,'C','U','V','R','N','N');
    tic_dgejsv_tmp1 = toc(tic_dgejsv_tmp1);
    tic_dgejsv_tmp2 = tic;
    dgejsv_mex(A,'C','U','V','R','N','N');
    tic_dgejsv_tmp2 = toc(tic_dgejsv_tmp2);
    time_dgejsv(i) = (tic_dgejsv_tmp1 + tic_dgejsv_tmp2)/2;

    % MATLAB SVD
    tic_matlab_tmp1 = tic;
    [~,~,~] = svd(A,'econ');
    tic_matlab_tmp1 = toc(tic_matlab_tmp1);
    tic_matlab_tmp2 = tic;
    [~,~,~] = svd(A,'econ');
    tic_matlab_tmp2 = toc(tic_matlab_tmp2);
    t_matlab(i) = (tic_matlab_tmp1 + tic_matlab_tmp2)/2;

    fprintf("Finished %d of %d \n", i, length(n));

end

savedata = 0;
if savedata == 1
    save('../data/svals_timing.mat');
end

%%
close all; 

C1 = "#1171BE";
C2 = "#DD5400";
C3 = "#EDB120";
C4 = "#3BAA32";

figure(1);

% a special potential time for mposj 
time_mposj_potential = time_mposj - mposj_disect_timing(:,2) + 100 * time_mposj_twoprec(:,2); 

%% Export the plotted data for the thesis figure
% Absolute seconds, as thesis panel (d) shows. The MATLAB figure below
% divides several of these by time_dgejsv; those ratios are not exported.
% The populated MATLAB timing variable is t_matlab, not time_matlab.
csv_header = {'n','mp3_s','dgesvj_s','dgejsv_s','matlab_s', ...
    'precond_ratio','mp3_jacobi_s','potential_s'};
csv_values = [n(:), time_mposj(:), time_dgesvj(:), time_dgejsv(:), t_matlab(:), ...
    mposj_disect_timing(:,2) ./ time_mposj_twoprec(:,2), ...
    mposj_disect_timing(:,3), time_mposj_potential(:)];
write_figure_csv("../data/svals_timing.csv", csv_values, csv_header);


loglog(n, time_mposj,'LineStyle','none','Marker','*','Color',C1); hold on;
loglog(n, time_dgesvj,'LineStyle','none','Marker','pentagram','Color',C2);
loglog(n, time_dgejsv,'LineStyle','none','Marker','square','Color',C3);
loglog(n, t_matlab,'LineStyle','none','Marker','diamond','Color',C4);

set(findall(gcf, 'Type', 'Line'), 'LineWidth', 1);
xlabel('Number of columns', 'FontSize', 10); 
ylabel('Runtime (sec)', 'FontSize', 10); 

% legend({'MP3JacobiSVD', ...
%     '\texttt{DGESVJ}', ...
%     '\texttt{DGEJSV}', ...
%     'MATLAB \texttt{svd}'}, ...
%     'Location', 'southeast', 'FontSize', 10);
% legend('boxoff')

%
figure(2)
loglog(n, mposj_disect_timing(:,2)./time_mposj_twoprec(:,2),'LineStyle','none','Marker', '<', 'Color', 'Black');
xlabel('Number of columns', 'FontSize', 10); 
ylabel('Runtime (sec)', 'FontSize', 10); 

%
figure(3)
loglog(n, mposj_disect_timing(:,3),'LineStyle','none','Marker','*','Color',C1); hold on; 
loglog(n, time_dgesvj,'LineStyle','none','Marker','pentagram','Color',C2);
xlabel('n');
ylabel('Runtime (sec)')

%
figure(4)
loglog(n, time_mposj_potential./time_dgejsv,'LineStyle','none','Marker', '>', 'Color', 'Black'); hold on;
loglog(n, time_dgesvj./time_dgejsv,'LineStyle','none','Marker','pentagram','Color',C2);
xlabel('n');
ylabel('Runtime (sec)')

fprintf("[END] svals_timing\n");
