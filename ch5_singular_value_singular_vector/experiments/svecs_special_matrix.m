%TEST5 - Access the accuracy of the computed singular vectors for special
addpath("../shared/");
%matrices

clear; close all;
fprintf("[START] svecs_special_matrix\n");

epsln = eps('double')/2;
writecsv = 1;
drawgraph = 0;

if writecsv && ~isfolder("../data")
    mkdir("../data");
end

for id = 1:2
    A = get_testmatrix_kms_lehmer(id); 
    [m,n] = size(A); 
    fprintf("[svecs_special_matrix] Starting matrix %d of 2 (%d x %d).\n", id, m, n);

    tol1 = sqrt(m*n)*epsln;

    [~,S1,V1,~,scnd] = mposj_svec(A);
    [Vref, Sref] = get_reference(A, 's');
    [err1, relgap1] = compute_error_svecs(A, V1, 's', Vref, Sref);

    bound = (1+scnd./relgap1)*tol1;

    [U2,S2,V2,sva2,work2,info2] = dgesvj_mex(A,'G','U','V',m,eye(n),max(6,m+n));
    if info2 ~= 0
        fprintf("Error: DGESVJ does not converge.\n");
        break;
    end
    [err2, ~ ] = compute_error_svecs(A, V2, 's', Vref, Sref);
    
    % preconditioned one-sided Jacobi
    [U3,S3,V3,sva3,work3,iwork3,info3] = dgejsv_mex(A,'C','U','V','R','N','N');
    if info3 ~= 0
        fprintf("Error: DGEJSV does not converge.\n");
        break;
    end
    [err3, ~ ] = compute_error_svecs(A, V3, 's', Vref, Sref);
    
    % MATLAB 
    [U4,S4,V4] = svd(A,'econ');
    [err4, ~] = compute_error_svecs(A, V4, 's', Vref, Sref);

    if writecsv
        if id == 1
            plotidx = 1:n;
            csvfile = "../data/svecs_special_matrix_kms.csv";
        else
            plotidx = 1:10:n;
            csvfile = "../data/svecs_special_matrix_lehmer.csv";
        end
        plotdata = table(plotidx(:), err1(plotidx).', err2(plotidx).', ...
            err3(plotidx).', err4(plotidx).', bound(plotidx).', ...
            'VariableNames', {'index', 'mp3', 'dgesvj', 'dgejsv', 'matlab', 'bound'});
        write_figure_csv(csvfile, plotdata);
    end

    if drawgraph
    figure(id); 
    C1 = "#1171BE";
    C2 = "#DD5400";
    C3 = "#EDB120";
    C4 = "#3BAA32";
    if id == 2


        semilogy(1:10:n,err1(1:10:n),'LineStyle','none','Marker','*', 'Color', C1);
        hold on;
        semilogy(1:10:n,err2(1:10:n),'LineStyle','none','Marker','pentagram', 'Color', C2);
        semilogy(1:10:n,err3(1:10:n),'LineStyle','none','Marker','square', 'Color', C3);
        semilogy(1:10:n,err4(1:10:n),'LineStyle','none','Marker','diamond', 'Color', C4);
        
        semilogy(1:10:n,bound(1:10:n),'LineStyle',':','Marker', 'none', 'Color','k');
    else
          
        semilogy(1:n,err1,'LineStyle','none','Marker','*', 'Color', C1);
        hold on;
        semilogy(1:n,err2,'LineStyle','none','Marker','pentagram', 'Color', C2);
        semilogy(1:n,err3,'LineStyle','none','Marker','square', 'Color', C3);
        semilogy(1:n,err4,'LineStyle','none','Marker','diamond', 'Color', C4);
        
        semilogy(1:n,bound,'LineStyle',':','Marker', 'none', 'Color','k');

    end
    set(findall(gcf, 'Type', 'Line'), 'LineWidth', 1);

    xlabel('Singular vector index (ordered by descending singular value)')
    ylabel('max svec err')
    
    legend('MP3JacobiSVD', 'DGESVJ', 'DGEJSV', 'MATLAB svd');
    end
        
    fprintf("[svecs_special_matrix] Finished matrix %d of 2.\n", id);
end

%%
printout = 0;
if printout 
    cleanfigure;
    % EDIT THIS PATH to your own output directory
    matlab2tikz('tmp.tex');
end

fprintf("[END] svecs_special_matrix\n");
