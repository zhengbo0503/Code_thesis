%TEST1 - forward accuracy test, changing condition number
addpath("../shared/");

close all; clear; rng(1);
fprintf("[START] svals_varying_cond\n");

m = 1000;
n = 800;
epsln = eps('double')/2;
kkappa = logspace(3,15,20);
f1 = zeros(length(m), 5);
f2 = zeros(length(m), 5);
f3 = zeros(length(m), 5);
f4 = zeros(length(m), 5);
bound1 = zeros(length(m),5);
bound2 = zeros(length(m),5);

for mode = 1:5

    for i = 1:length(kkappa)

        mm = m;
        nn = n;
        kappa = kkappa(i); 

        A = gallery('randsvd', [mm,nn], kappa, mode);
        sref = reference_singular_values(A);

        [U1,S1,V1,nos1,scnd] = mposj(A, 3, true);
        [f1(i,mode),~,~,~] = compute_error_svals(A, U1, S1, V1, sref);

        [U2,S2,V2,sva2,work2,info2] = dgesvj_mex(A,'G','U','V',nn,eye(nn),max(6,mm+nn));
        if info2 ~= 0
            fprintf("Error: DGESVJ does not converge.\n");
            break;
        end
        [f2(i,mode),~,~,~] = compute_error_svals(A, U2, S2, V2, sref);

        [U3,S3,V3,sva3,work3,iwork3,info3] = dgejsv_mex(A,'C','U','V','R','N','N');
        if info3 ~= 0
            fprintf("Error: DGEJSV does not converge.\n");
            break;
        end
        [f3(i,mode),~,~,~] = compute_error_svals(A, U3, S3, V3, sref);

        [U4,S4,V4] = svd(A,'econ');
        [f4(i,mode),~,~,~] = compute_error_svals(A, U4, S4, V4, sref);

        bound1(i,mode) = scond(A, 'C') * sqrt(mm * nn)* epsln;
        bound2(i,mode) = scnd * sqrt(mm * nn) * epsln;

        fprintf("Finished MODE = %d, %d of %d\n", mode, i, length(kkappa));

    end
end

savedata = 1;
if savedata == 1
    save("../data/svals_varying_cond.mat");
end

%% Export the plotted data for the thesis figure
% bound2 (preconditioned matrix) is the thesis reference curve; bound1 is
% drawn by the MATLAB figure only and is deliberately not exported.
csv_header = {'kappa', 'mp3_m1','dgesvj_m1','dgejsv_m1','matlab_m1','bound_m1', 'mp3_m2','dgesvj_m2','dgejsv_m2','matlab_m2','bound_m2', 'mp3_m3','dgesvj_m3','dgejsv_m3','matlab_m3','bound_m3', 'mp3_m4','dgesvj_m4','dgejsv_m4','matlab_m4','bound_m4', 'mp3_m5','dgesvj_m5','dgejsv_m5','matlab_m5','bound_m5'};
csv_values = [kkappa(:), f1(:,1), f2(:,1), f3(:,1), f4(:,1), bound2(:,1), f1(:,2), f2(:,2), f3(:,2), f4(:,2), bound2(:,2), f1(:,3), f2(:,3), f3(:,3), f4(:,3), bound2(:,3), f1(:,4), f2(:,4), f3(:,4), f4(:,4), bound2(:,4), f1(:,5), f2(:,5), f3(:,5), f4(:,5), bound2(:,5)];
write_figure_csv("../data/svals_varying_cond.csv", csv_values, csv_header);

%% 
close all; 
C1 = "#1171BE";
C2 = "#DD5400";
C3 = "#EDB120";
C4 = "#3BAA32";
    

for mode = 1:5
    figure(mode)
    
    % %%%%%%%%%%%%%%%%%%
    loglog(kkappa,f1(:,mode),'LineStyle','none','Marker','*', 'Color', C1); hold on;
    loglog(kkappa,f2(:,mode),'LineStyle','none','Marker','pentagram', 'Color', C2);
    loglog(kkappa,f3(:,mode),'LineStyle','none','Marker','square', 'Color', C3);
    loglog(kkappa,f4(:,mode),'LineStyle','none','Marker','diamond', 'Color', C4);
    loglog(kkappa,bound1(:,mode),'LineStyle','--','Marker', 'none', 'Color','k');
    loglog(kkappa,bound2(:,mode),'LineStyle',':','Marker', 'none', 'Color','k');
    set(findall(gcf, 'Type', 'Line'), 'LineWidth', 1);
    % %%%%%%%%%%%%%%%%%%

    axis square
    
    xlabel('Condition number $\kappa_2(A)$','FontSize',10);
    xlim([1e3,1e15]);
    xticks([1e3,1e6,1e9,1e12,1e15]);
    ylim([1e-16,1e2]);
    yticks([1e-16, 1e-13, 1e-10,1e-7,1e-4, 1e-1, 1e2]);
    
    % Only get ylabel if the subplot is on the left
    if mod(mode,2) == 1
        ylabel('$\mathrm{max}_k {\varepsilon}^{(k)}_{fwd}$','FontSize',10);
    end

    % Set the title 
    alphabet = ['a','b','c','d','e'];
    t = sprintf('(%s) MODE = %d', alphabet(mode), mode); 
    title(t, 'FontWeight', 'normal', 'FontSize', 10); 
end

fprintf("[END] svals_varying_cond\n");
