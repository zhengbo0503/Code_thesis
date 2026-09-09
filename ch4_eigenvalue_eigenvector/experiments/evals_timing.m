clear; close all;
fprintf("[START] evals_timing\n");
addpath("../shared/");

% Load data from Julia
data_order = readmatrix("../julia/result/timing_order.csv");
data_kappa = readmatrix("../julia/result/timing_cnd.csv");

% Separate data_order
N = data_order(:,1);
order_tj = data_order(:,2);
order_tm2 = data_order(:,3);
order_tm3 = data_order(:,4);
order_tm2Prec = data_order(:,5);
order_tm2Apply = data_order(:,6);
order_tm2Jacobi = data_order(:,7);
order_tm2Else = data_order(:,8);
order_tm3Prec = data_order(:,9);
order_tm3Apply = data_order(:,10);
order_tm3Jacobi = data_order(:,11);
order_tm3Else = data_order(:,12);

% Separate data_kappa 
kappa = data_kappa(:,1);
kappa_tj = data_kappa(:,2);
kappa_tm2 = data_kappa(:,3);
kappa_tm3 = data_kappa(:,4);
kappa_tm2Prec = data_kappa(:,5);
kappa_tm2Apply = data_kappa(:,6);
kappa_tm2Jacobi = data_kappa(:,7);
kappa_tm2Else = data_kappa(:,8);
kappa_tm3Prec = data_kappa(:,9);
kappa_tm3Apply = data_kappa(:,10);
kappa_tm3Jacobi = data_kappa(:,11);
kappa_tm3Else = data_kappa(:,12);

%% Figure 1 for grand timings 
close all; figure(1); 
subplot(1,2,1);
plot(N, order_tj,"-o");
hold on;
plot(N, order_tm2, "-x");
plot(N, order_tm3, "->");
xlim([100,1e3]);
xlabel("$n$", Interpreter="latex")
xticks([100,300,500,700,1000]);
ylabel("Time (s)")
axis square
grid on

subplot(1,2,2);
semilogx(kappa, kappa_tj, "-o");
hold on;
semilogx(kappa, kappa_tm2, "-x");
semilogx(kappa, kappa_tm3, "->");
xlim([1e3,1e14]);
xlabel("$\kappa_2(A)$", Interpreter="latex")
xticks([1e3,1e6,1e9,1e12,1e14])
axis square
grid on

legend("Jacobi $\quad$", "MP2Jacobi $\quad$", "MP3Jacobi $\quad$", Location="north", NumColumns=3);

%% Figure 2 for Decomposing the timing for Orders

% Apply preconditioner
figure(2)
qdratio = 100; 

% Apply preconditioner 
subplot(2,2,1); 
plot(N,order_tm3Apply./order_tm2Apply,"->", "Color", "#0072BD");
hold on;
plot(N,(order_tm2Apply*qdratio)./order_tm2Apply, "--*", "Color", "#77AC30");
xlim([100,1e3]);
xlabel("$n$", Interpreter="latex")
xticks([100,300,500,700,1000]);
ylabel("Relative timings against MP2Jacobi")
title("(a) Applying the preconditioner","FontWeight","normal");
axis square
grid on
legend("MP3Jacobi", "MP3Jacobi (Potential)", "location", "best")
 
% Apply Jacobi 
subplot(2,2,2);
plot(N, order_tm2Jacobi,"-x", "Color", "#D95319");
hold on;
plot(N, order_tm3Jacobi, "->", "Color", "#0072BD");
xlim([100,1e3]);
xlabel("$n$", Interpreter="latex")
xticks([100,300,500,700,1000]);
ylabel("Timings")
title("(b) Applying the Jacobi algorithm","FontWeight","normal");
axis square
grid on
legend("MP2Jacobi", "MP3Jacobi", "location", "northwest") 

% Grand all 
subplot(2,2,3);
plot(N,order_tm3./order_tj,"->", "Color", "#0072BD");
hold on;
plot(N,(order_tm2Apply*(qdratio-1)+order_tm2)./order_tj, "--*", "Color", "#77AC30");
plot(N,order_tm2./order_tj,"-x", "Color", "#D95319");  
xlim([100,1e3]);
xlabel("$n$", Interpreter="latex")
xticks([100,300,500,700,1000]);
ylabel("Relative timings against Jacobi")
ylim([0.1,2]);
title("(c) Total time","FontWeight","normal");
axis square
grid on
legend("MP3Jacobi", "MP3Jacobi (Potential)", "MP2Jacobi", "location", "northwest") 


%% Export the plotted data as CSV for the thesis TikZ figures
% Consumers: figs/timing_union.tex (fig.ev-timing-union) and
%            figs/timings_2.tex   (fig.ev-timing-all)
% The TeX reads these columns directly; do not rename them.
% Routed through write_figure_csv so this shares the 17-digit, atomic
% write-and-publish path used by every other figure export.

assert(numel(N) == numel(kappa), 'Independent grids have different lengths.');
assert(all(order_tm2Apply > 0) && all(order_tj > 0), ...
    'Normalization denominators must be positive.');
assert(qdratio == 100, 'The fixed thesis caption requires qdratio = 100.');

union_header = {'n','kappa','t_jacobi_n_s','t_mp2_n_s','t_mp3_n_s', ...
    't_jacobi_kappa_s','t_mp2_kappa_s','t_mp3_kappa_s'};
union_values = [N, kappa, order_tj, order_tm2, order_tm3, ...
    kappa_tj, kappa_tm2, kappa_tm3];

breakdown_header = {'n','r_precond_mp3_mp2','r_precond_potential_mp2', ...
    't_jacobi_stage_mp2_s','t_jacobi_stage_mp3_s', ...
    'r_total_mp2_jacobi','r_total_mp3_jacobi','r_total_potential_jacobi'};
breakdown_values = [N, ...
    order_tm3Apply ./ order_tm2Apply, ...
    qdratio * ones(size(N)), ...
    order_tm2Jacobi, order_tm3Jacobi, ...
    order_tm2 ./ order_tj, ...
    order_tm3 ./ order_tj, ...
    (order_tm2Apply * (qdratio - 1) + order_tm2) ./ order_tj];

write_figure_csv("../output/timing_union.csv", union_values, union_header);
write_figure_csv("../output/timings_2.csv", breakdown_values, breakdown_header);

fprintf("[END] evals_timing\n");
