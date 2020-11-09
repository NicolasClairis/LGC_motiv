function [figH] = do_boxplot_cell(Xcell, xlabels, boxColor, xTitle, yTitle, genTitle, enableLatex, fullScreen)
% 
% [figH] = do_boxplot(Xcell, xlabels, boxColor, xTitle, yTitle, genTitle, enableLatex, fullScreen)
%
% This function plots Xcell cell array (1xp size) as p boxplots
% and returns the figure handle figH
% -------------------------------------------------------------------------

    figH = figure;
    X = cell2mat_pad(Xcell);
   
  % Set display options
  % ---
    if nargin > 6 && enableLatex
      % Enable LateX interpreter  
        set(gca, 'TickLabelInterpreter', 'latex');
    end

    if nargin > 7  && fullScreen
      % Enlarge figure to full-size
        set(gcf, 'Position', get(0, 'Screensize'));
    end
    
    
  % Draw '0' baseline (dashed)
  % ---
    plot([0 length(xlabels)+1],[0 0], '--', 'Color', [.7 .7 .7], 'LineWidth', 1.1)
    hold on
    
    
  % Set titles  
  % ---
    if nargin > 3 && ~isempty(xTitle)
        xlabel(xTitle)
    end
    
    if nargin > 4 && ~isempty(yTitle)
        ylabel(yTitle);
    end
    
    if nargin > 5 && ~isempty(genTitle)
        title(genTitle, 'Interpreter', 'none')
    end
    
    
    
  % Draw box plot
  % ---
  
    if nargin < 3 || isempty(boxColor)
  	    boxColor = [.2 .9 .6]; % default color
    end
  
    iosr.statistics.boxPlot(xlabels, X, ...
        'notch', false, ...
        'notchDepth', .2, ...
        'lineColor', [.6 .6 .6], ...
        'lineStyle', ':', ...
        'lineWidth', .8, ...
        'medianColor', [.4 .4 .4], ...
        'medianWidth', 2, ....
        'boxColor', boxColor, ...
        'boxAlpha', .2, ...
        'boxWidth', .3, ....
        'showScatter', true, ...
        'scatterAlpha', .1,...
        'scatterColor', boxColor, ...
        'scatterMarker', '.', ...
        'scatterSize', 300, ...
        'symbolColor', boxColor, ...
        'outlierSize', 12);
    
    box off
  
    
end