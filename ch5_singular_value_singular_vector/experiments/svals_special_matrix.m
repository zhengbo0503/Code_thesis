%TEST3 -- accuracy test for MPOSJ on matrices from the Anymatrix collection.

clear; close all; rng(1);
fprintf("[START] svals_special_matrix\n");

C1 = "#1171BE";
C2 = "#DD5400";
C3 = "#EDB120";
C4 = "#3BAA32";
epsln = eps('double')/2;
addpath("../shared/");

for i = 1:4
    A = get_testmatrix_anymatrix(i);
    [mm,nn] = size(A); 
    fprintf("[svals_special_matrix] Starting matrix %d of 4 (%d x %d).\n", i, mm, nn);
    n(i) = min(mm,nn); 
    m(i) = max(mm,nn); 
    
    % %%%
    Sref = reference_singular_values(A);
    
    [U1,S1,V1,nos1,scnd] = mposj(A, 3, true);
    S1 = reshape(sort(diag(S1), 'descend'), size(Sref,1), size(Sref,2));
    err1(1:min(mm,nn),i) = abs(S1 - Sref)./Sref;
    
    [U2,S2,V2,sva2,work2,info2] = dgesvj_mex(A,'G','U','V',nn,eye(nn),max(6,mm+nn));
    if info2 ~= 0
        fprintf("Error: DGESVJ does not converge.\n");
    end
    S2 = reshape(sort(diag(S2), 'descend'), size(Sref,1), size(Sref,2));
    err2(1:min(mm,nn),i) = abs(S2 - Sref)./Sref;
    
    [U3,S3,V3,sva3,work3,iwork3,info3] = dgejsv_mex(A,'C','U','V','R','N','N');
    if info3 ~= 0
        fprintf("Error: DGEJSV does not converge.\n");
    end
    S3 = reshape(sort(diag(S3), 'descend'), size(Sref,1), size(Sref,2));
    err3(1:min(mm,nn),i) = abs(S3 - Sref)./Sref;
    
    [U4,S4,V4] = svd(A,'econ');
    S4 = reshape(sort(diag(S4), 'descend'), size(Sref,1), size(Sref,2));
    err4(1:min(mm,nn),i) = abs(S4 - Sref)./Sref;
    
    bound1(1:min(mm,nn),i) = scond(A, 'C') * sqrt(mm * nn)* epsln * ones(min(mm,nn),1);
    bound2(1:min(mm,nn),i) = scnd * sqrt(mm * nn) * epsln * ones(min(mm,nn),1);

    scond_store(i) = scond(A, 'C'); 
    scnd_store(i) = scnd; 
    fprintf("[svals_special_matrix] Finished matrix %d of 4.\n", i);
end

%% Export the plotted data for the thesis figure
% Per-matrix subsampling matches the thesis panels. The whiskycorr stress
% test deliberately has no theoretical reference curve, so it has no bound
% column even though the script computes one.
csv_files = {"../data/svals_special_matrix_blur.csv", ...
             "../data/svals_special_matrix_kahan.csv", ...
             "../data/svals_special_matrix_whiskycorr.csv", ...
             "../data/svals_special_matrix_lauchli.csv"};
for i = 1:4
    if i == 4
        p = (1:10:n(i)).';
    else
        p = (1:n(i)).';
    end
    if i == 3
        csv_header = {'index','mp3','dgesvj','dgejsv','matlab'};
        csv_values = [p, err1(p,i), err2(p,i), err3(p,i), err4(p,i)];
    else
        csv_header = {'index','mp3','dgesvj','dgejsv','matlab','bound'};
        csv_values = [p, err1(p,i), err2(p,i), err3(p,i), err4(p,i), bound2(p,i)];
    end
    write_figure_csv(csv_files{i}, csv_values, csv_header);
end

%%
close all; 
title_candi = ["anymatrix('regtools/blur', 10, 5, 1.4)", 
    "anymatrix('gallery/kahan',50,1e-2)", 
    "anymatrix('nessie/whiskycorr')", 
    "anymatrix('gallery/lauchli', 500, 1e-3)"];
title_alphabet = ['a', 'b', 'c', 'd'];

for i = 1:4
    bound = scnd_store(i)*ones(n(i),1)*sqrt(n(i)*m(i))*epsln;
    
    figure(i);
    if i ~= 4
        semilogy(1:n(i), err1(1:n(i),i),'LineStyle','none','Marker','*','Color',C1); hold on; 
        semilogy(1:n(i), err2(1:n(i),i),'LineStyle','none','Marker','pentagram','Color',C2);
        semilogy(1:n(i), err3(1:n(i),i),'LineStyle','none','Marker','square','Color',C3);
        semilogy(1:n(i), err4(1:n(i),i),'LineStyle','none','Marker','diamond','Color',C4);
        semilogy(1:n(i), bound, 'LineStyle',':','Marker','none','Color','k'); 
    else
        semilogy(1:10:n(i), err1(1:10:n(i),i),'LineStyle','none','Marker','*','Color',C1); hold on; 
        semilogy(1:10:n(i), err2(1:10:n(i),i),'LineStyle','none','Marker','pentagram','Color',C2);
        semilogy(1:10:n(i), err3(1:10:n(i),i),'LineStyle','none','Marker','square','Color',C3);
        semilogy(1:10:n(i), err4(1:10:n(i),i),'LineStyle','none','Marker','diamond','Color',C4);
        semilogy(1:10:n(i), bound(1:10:n(i)), 'LineStyle',':','Marker','none','Color','k');
    end

    axis square 
    title_tmp = sprintf("(%s) %s", title_alphabet(i), title_candi(i));
    title(title_tmp, 'FontWeight', 'normal', 'FontSize', 10); 
    set(findall(gcf, 'Type', 'Line'), 'LineWidth', 1);
    xlabel('$k$th largest singular value', 'FontSize', 10, 'Interpreter', 'latex'); 
    if mod(i,2) == 1
        ylabel('$\varepsilon_{fwd}^{(k)}$', 'FontSize', 10, 'Interpreter', 'latex');
    end

    % if (i ~= 3)
    %     ylim([1e-18,1e-8]); 
    %     yticks([1e-18,1e-16,1e-14,1e-12,1e-10,1e-8]);
    % end
    % 
    if (i == 4)
        ylim([1e-18, 1e-8]);
    end

    
    % legend({'MP3JacobiSVD', '\texttt{DGESVJ}', '\texttt{DGEJSV}', 'MATLAB \texttt{svd}'}, 'FontSize', 10);

end

fprintf("[END] svals_special_matrix\n");
