 function [stats, figH] = ttestAndBoxPlot(B, saveFilename, labels, titleGraph, color)

% function [stats, figH] = ttestAndPlot(B, saveFilename, labels, titleGraph, color)

% This function performs a ttest on B array of betas (1 row = 1 subject, 1
% column = 1 regressor), saves stats to 'saveFilename' (/!\ do not include .mat extension in filename), 
% then draws a box plot with labels 'labels' (cell array of char strings),
% and optionally specifies regression options 'optionNames' (cell array of strings) and 'optionValues' (array of doubles)
%  OUT: stats = some useful statistics about Betas
%       figH = figure window handle 
%       bH = bar graph handle
%       axH = axis handle
% EB 03/18

%%
 % t-test beta Group level
 % ---------------------------

     numRegX = length(labels);
     stats.pGroup = nan(1, numRegX);

     stats.meanBetasGroup = nanmean(B, 1);
     stats.semBetasGroup = sem(B, 1); 

     for iRegX = 1:numRegX
       [~, stats.pGroup(iRegX)] = ttest(B(:,iRegX));
     end
     
   labelsComplete = cell(1,length(labels));
   for i = 1:length(labels)
      pVal = num2str(stats.pGroup(i), '%.2g');
      labelsComplete{i} = ['\begin{tabular}{c}', labels{i}, '\\', pVal, '\end{tabular}'];
   end

     
 % Draw boxplot
 % ---
 
    figH = figure;
  % Draw '0' baseline (dashed)
    plot([0 numRegX + 1],[0 0], '--', 'Color', [.7 .7 .7], 'LineWidth', 1.1)
    ylabel('Beta estimates');
    hold on
    
  % Enable LateX interpreter  
    set(gca, 'TickLabelInterpreter', 'latex');%,  'HorizontalAlignment', 'center');
    
  % Enlarge figure to full-size
    set(gcf, 'Position', get(0, 'Screensize'));
    
  % Draw box plot
  
    [~, outlierIdx] = min(B(1,:));
    B(1,outlierIdx) = NaN;
    iosr.statistics.boxPlot(labelsComplete, B, ...
        'notch', false, ...
        'notchDepth', .2, ...
        'lineColor', [.6 .6 .6], ...
        'lineStyle', ':', ...
        'lineWidth', .8, ...
        'medianColor', [.4 .4 .4], ...
        'medianWidth', 2, ....
        'boxColor', color, ...
        'boxAlpha', .2, ...
        'boxWidth', .3, ....
        'showScatter', true, ...
        'scatterColor', color, ...
        'scatterMarker', '.', ...
        'scatterSize', 100, ...
        'symbolColor', color, ...
        'outlierSize', 12);
    
    box off
    
    title(titleGraph);
  
  % Make it wide  
    pbaspect([3 1 1])
    
  % Save plot   
    print(saveFilename, '-dpng', '-r600')
    
    
 end
