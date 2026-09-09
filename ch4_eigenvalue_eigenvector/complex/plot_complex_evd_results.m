function files = plot_complex_evd_results(results,kind,plotDir)
%PLOT_COMPLEX_EVD_RESULTS - Export plots from complex-experiment results
%
%   Usage:
%       files = plot_complex_evd_results(results, kind)
%       files = plot_complex_evd_results(results, kind, plotDir)
%
%   Purpose:
%       PLOT_COMPLEX_EVD_RESULTS creates the eigenvalue accuracy
%       figure associated with a saved complex-experiment result structure.
%       Figures are rendered invisibly and exported in PNG and vector PDF
%       formats. The output directory is created when necessary.
%
%   Arguments:
%       (1) results - Result structure returned by a complex experiment.
%       (2) kind - 'varying_kappa' or 'varying_n'.
%       (3) plotDir - Output directory, and by default the plots directory
%           under complex_experiments.
%
%   Outputs:
%       (1) files - Two-element cell array containing the PNG and PDF paths.
%

kind = validatestring(kind,{'varying_kappa','varying_n'});
if nargin < 3 || isempty(plotDir)
    plotDir = fullfile(fileparts(mfilename('fullpath')),'plots');
end
if ~isfolder(plotDir), mkdir(plotDir); end

fileBase = fullfile(plotDir,['complex_evals_' kind '_' results.profile]);
if isfield(results,'boundAt')
    boundAt = results.boundAt;
else
    boundAt = [];
end
plainNotConverged = false(size(results.forward,1),numel(results.modes));
if isfield(results,'info')
    plainNotConverged = results.info(:,:,2) ~= 0;
end
switch kind
    case 'varying_kappa'
        finiteKappas = results.actualKappas(isfinite(results.actualKappas));
        kappaLimits = [0.95*min(finiteKappas),1.05*max(finiteKappas)];
        make_accuracy_plot(results.actualKappas,results.modes,results.forward, ...
            results.methods,'Condition number $\kappa_2(A)$', ...
            ['Complex two-sided Jacobi: varying $\kappa$ (' results.profile ')'], ...
            boundAt,plainNotConverged,kappaLimits,[1e3 1e6 1e9 1e12 1e15],fileBase);
    case 'varying_n'
        make_accuracy_plot(results.nvalues,results.modes,results.forward, ...
            results.methods,'$n$', ...
            ['Complex two-sided Jacobi: varying $n$ (' results.profile ')'], ...
            boundAt,plainNotConverged,[],[],fileBase);
end
files = {[fileBase '.png'],[fileBase '.pdf']};
end

function make_accuracy_plot(x,modes,forward,methods,xlabelText,superTitle,boundAt,plainNotConverged,xLimits,xTicks,fileBase)
%MAKE_ACCURACY_PLOT - Export the tiled varying-kappa or varying-n figure.
colours = [17 113 190;221 84 0;237 177 32;59 170 50]/255;
markers = {'*','p','s','d'};
fig = figure('Color','w','Visible','off','Position',[100 100 1120 650]);
tl = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
handles = gobjects(1,numel(methods));
boundHandle = gobjects(1);
nonconvergedHandle = gobjects(1);
for jm = 1:numel(modes)
    ax = nexttile(tl);
    set(ax,'XScale','log','YScale','log');
    hold(ax,'on');
    if isvector(x)
        xvalues = x;
    else
        xvalues = x(:,jm);
    end
    for k = 1:numel(methods)
        handles(k) = plot(ax,xvalues,forward(:,jm,k),'LineStyle','none', ...
            'Marker',markers{k},'Color',colours(k,:),'LineWidth',1.1);
    end
    % Overlay the retained baseline iterates that did not meet stopping.
    failed = plainNotConverged(:,jm);
    if any(failed)
        jacobiErrors = forward(:,jm,2);
        nonconvergedHandle = plot(ax,xvalues(failed),jacobiErrors(failed), ...
            'LineStyle','none','Marker','x','Color','k','MarkerSize',10, ...
            'LineWidth',1.5);
    end
    if ~isempty(boundAt)
        boundHandle = plot(ax,xvalues,boundAt(:,jm),'LineStyle',':', ...
            'Marker','none','Color','k','LineWidth',1.1);
    end
    if ~isempty(xLimits)
        xlim(ax,xLimits);
        xticks(ax,xTicks);
    end
    grid(ax,'on'); axis(ax,'square');
    xlabel(ax,xlabelText,'Interpreter','latex');
    ylabel(ax,'$\max_k\varepsilon_{\rm fwd}^{(k)}$','Interpreter','latex');
    title(ax,sprintf('MODE = %d',modes(jm)),'FontWeight','normal');
end
legendHandles = handles;
legendLabels = methods;
if ~isempty(boundAt)
    legendHandles(end+1) = boundHandle;
    legendLabels{end+1} = '$7nu\,\kappa_D(\widetilde A)$ (real-case reference)';
end
if isgraphics(nonconvergedHandle)
    legendHandles(end+1) = nonconvergedHandle;
    legendLabels{end+1} = 'Jacobi: not converged (sweep limit)';
end
lgd = legend(legendHandles,legendLabels,'Location','northwest', ...
    'Interpreter','latex');
lgd.Layout.Tile = 6;
sgtitle(tl,superTitle,'Interpreter','latex');
exportgraphics(fig,[fileBase '.png'],'Resolution',220);
exportgraphics(fig,[fileBase '.pdf'],'ContentType','vector');
close(fig);
end
