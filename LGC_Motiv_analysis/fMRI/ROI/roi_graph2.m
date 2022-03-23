function[roi_fig] = roi_graph2(selectedContrastIndex,...
            con_avg1, con_errorbar1,...
            con_avg2, con_errorbar2,...
            figConName, ttest_pval_1_vs_2, xname1, xname2)
% [roi_fig] = roi_graph2(selectedContrastIndextrastIndex,...
%             con_avg1, con_sem1,...
%             con_avg2, con_sem2,...
%             figConName, ttest_pval_1_vs_2, name1, name2)
% roi_graph2 will display a bar graph showing 2 columns with a star if the
% difference between the two is significant.
%
% INPUTS
% selectedContrastIndextrastIndex: index for the contrast
%
% con_avg1: average for the first group of data
%
% con_errorbar1: error bar (SEM or SD) for the first group of data
%
% con_avg2: average for the second group of data
%
% con_errorbar2: error bar (SEM or SD) for the second group of data
%
% figConName: name for the contrast (display name for X)
%
% ttest_pval_1_vs_2: p.value for a t.test comparing the two groups
%
% xname1: name for first column
%
% xname2: name for second column
%
% OUTPUTS
%
% roi_fig: handle for the figure
%

%% extract general informations
langage = 'engl';
lSize = 3;
ftSize = 18;
pSize = 50;
roi_color = ([143, 0, 0])./255;
% roi_color = ([143 143 143])./255;

% default yscale values
minYvalue1 = min(min(min(min(con_avg1(selectedContrastIndex) - con_errorbar1(selectedContrastIndex) )))); % extract min beta value across subs and ROI
maxYvalue1 = max(max(max(max(con_avg1(selectedContrastIndex) + con_errorbar1(selectedContrastIndex) )))); % extract max beta value across subs and ROI
minYvalue2 = min(min(min(min(con_avg2(selectedContrastIndex) - con_errorbar2(selectedContrastIndex) )))); % extract min beta value across subs and ROI
maxYvalue2 = max(max(max(max(con_avg2(selectedContrastIndex) + con_errorbar2(selectedContrastIndex) )))); % extract max beta value across subs and ROI
minYvalue = min(minYvalue1, minYvalue2);
maxYvalue = max(maxYvalue1, maxYvalue2);
dist_min_max_Yvalue = abs(maxYvalue - minYvalue);
% define various possible default yscale values depending on minBeta and
% maxBeta
% stop scale at 0 if all contrasts are positive or negative
if minYvalue > 0 
    minYscale = 0;
else
    minYscale = minYvalue - dist_min_max_Yvalue/3;
end
if maxYvalue < 0 
    maxYscale = 0;
else
    maxYscale = maxYvalue + dist_min_max_Yvalue/3;
end
yscale = [minYscale,  maxYscale];
fig_size = abs(yscale(2) - yscale(1));
% position of line and stars for each column
if (con_avg1(selectedContrastIndex) > 0) || (con_avg2(selectedContrastIndex) > 0)
    dYpos_yscale = abs(yscale(2) - maxYvalue);
    pval_line_Ypos = maxYvalue + (1/3)*dYpos_yscale;
    pval_Ypos = maxYvalue + (4/5)*dYpos_yscale;
else
    dYpos_yscale = abs(yscale(1) - minYvalue);
    pval_line_Ypos = minYvalue - (1/3)*dYpos_yscale;
    pval_Ypos = minYvalue - (1/3)*dYpos_yscale;
end

%% display figure
scale_ok_idx = 0;
while scale_ok_idx == 0
    roi_fig = fig();

    %% loop through sequences
    xpos = (1:2)';

    bar_hdl = bar(xpos,...
        [con_avg1(selectedContrastIndex), con_avg2(selectedContrastIndex)]);
    bar_hdl.EdgeColor = roi_color;
    bar_hdl.FaceColor = roi_color;
    hold on;
    error_hdl = errorbar(xpos,...
        [con_avg1(selectedContrastIndex), con_avg2(selectedContrastIndex)],...
        [con_errorbar1(selectedContrastIndex), con_errorbar2(selectedContrastIndex)],...
        'LineStyle','none',...
        'LineWidth',lSize,'Color','k');

    % adapt p.value position depending on number of contrasts and number of sequences
    pval_xpos = 1.5;
    % add line on top (or below) bars and below (or above) stars
    if ttest_pval_1_vs_2(selectedContrastIndex) <= 0.05
        line(xpos, [pval_line_Ypos, pval_line_Ypos],...
            'Color','k','LineStyle','-','LineWidth',lSize);

        % add stars above (or below) bars and line
        if ttest_pval_1_vs_2(selectedContrastIndex) <= 0.05 && ttest_pval_1_vs_2(selectedContrastIndex) > 0.01
            text(pval_xpos, pval_Ypos, '*',...
                'HorizontalAlignment','center','VerticalAlignment', 'top', 'FontSize', ftSize)
        elseif ttest_pval_1_vs_2(selectedContrastIndex) <= 0.01 && ttest_pval_1_vs_2(selectedContrastIndex) > 0.005
            text(pval_xpos, pval_Ypos, '**',...
                'HorizontalAlignment', 'center','VerticalAlignment', 'top', 'FontSize', ftSize)
        elseif ttest_pval_1_vs_2(selectedContrastIndex) <= 0.005 && ttest_pval_1_vs_2(selectedContrastIndex) > 0.001
            text(pval_xpos, pval_Ypos, '***',...
                'HorizontalAlignment', 'center','VerticalAlignment', 'top', 'FontSize', ftSize)
        elseif ttest_pval_1_vs_2(selectedContrastIndex) <= 0.001
            text(pval_xpos, pval_Ypos, '****',...
                'HorizontalAlignment', 'center','VerticalAlignment', 'top', 'FontSize', ftSize)
        end % how many stars for p.value?
    end % is p.value significant?

    % define graph limits
    xlim([0 3]);
    ylim(yscale);

    % draw a line to show where the zero stands (if scale not thresholded
    % at zero)
    if (minYscale < 0) && (maxYscale > 0)
        line(xlim(),[0 0],'Color','k','LineStyle','-','LineWidth',lSize);
    end

    xticks(xpos);
    set(gca,'xtick',xpos,'XTickLabel',{xname1,xname2});
    title(figConName);
    switch langage
        case 'engl'
            ylabel('Regression estimate');
        case 'fr'
            ylabel('coefficient de régression');
    end

    %
    legend_size(pSize);
    
    %% check if graph scale is fine for proper display or not
    scale_fig = figure;
    set(scale_fig,'Position',[0,0,250,100]); % reduce new figure size
    scale_ok = spm_input('Y scale ok ?',1,'b','yes | no',[1 2]);
    if scale_ok == 1 % if scale is fine
        scale_ok_idx = 1;
    elseif scale_ok == 2
        yscale = spm_input('New scale values [y1 (low) y2 (high)]?',1,'e');
        fig_size = abs(yscale(2) - yscale(1));
        close(roi_fig); % close graph with previous values if not ok
    end
    close(scale_fig);
    
end % scale

end % function