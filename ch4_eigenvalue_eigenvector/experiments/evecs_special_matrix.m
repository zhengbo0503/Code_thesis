%TEST6 - Access the accuracy of the computed eigenvectors for special
addpath("../shared/");
%matrices 
clear; close all;
fprintf("[START] evecs_special_matrix\n");

epsln = eps('double')/2;
writecsv = 1;
drawgraph = 0;

if writecsv && ~isfolder("../data")
    mkdir("../data");
end

for id = 1:2
    A = get_testmatrix(id); 
    [m,n] = size(A); 
    fprintf("[evecs_special_matrix] Starting matrix %d of 2 (%d x %d).\n", id, m, n);

    tol1 = 5*sqrt(n)*epsln;

    [V1,D1,~,~,~,scnd] = mp_pjacobi(A,'mp3');
    [Vref, Sref] = get_reference(A, 'E');
    [err1, relgap1] = compute_error_evecs(A, V1, 'E', Vref, Sref);

    bound = (1+scnd./relgap1)*tol1;

    [V2,D2] = cjacobi(A);
    [err2, ~ ] = compute_error_evecs(A, V2, 'E', Vref, Sref); 

    % Preconditioned two-sided Jacobi
    [V3,D3] = mp_pjacobi(A,'mp2');
    [err3, ~ ] = compute_error_evecs(A, V3, 'E', Vref, Sref); 

    % MATLAB
    [V4,D4] = eigsort(A);    
    [err4, ~ ] = compute_error_evecs(A, V4, 'E', Vref, Sref);

    if writecsv
        if id == 1
            plotidx = 1:n;
            csvfile = "../data/evecs_special_matrix_kms.csv";
        else
            plotidx = 1:10:n;
            csvfile = "../data/evecs_special_matrix_lehmer.csv";
        end
        plotdata = table(plotidx(:), err1(plotidx).', err2(plotidx).', ...
            err3(plotidx).', err4(plotidx).', bound(plotidx).', ...
            'VariableNames', {'index', 'mp3', 'jacobi', 'mp2', 'matlab', 'bound'});
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

    xlabel('Eigenvector index (ordered by descending eigenvalue)')
    ylabel('max evec err')

    legend('MP3Jacobi', 'Jacobi', 'MP2Jacobi', 'MATLAB eig');
    end
    
    fprintf("[evecs_special_matrix] Finished matrix %d of 2.\n", id);
end

%%
printout = 0;
if printout 
    cleanfigure;
    % EDIT THIS PATH to your own output directory
    matlab2tikz('tmp.tex');
end

fprintf("[END] evecs_special_matrix\n");
