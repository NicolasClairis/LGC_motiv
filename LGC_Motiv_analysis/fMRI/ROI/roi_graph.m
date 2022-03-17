function[roi_fig] = roi_graph(selectedCon, n_ROIs,...
        con_avg, con_sem, con_names, ttest_pval)
% [roi_fig] = roi_graph(selectedCon, n_ROIs, con_avg, con_sem, con_names, ttest_pval);
% roi_graph will represent the data for the selected ROI
%
% INPUTS
% selectedCon: index for selected contrasts
%
% n_ROIs: number of regions of interest (ROIs)
%
% con_vec_all: matrix with contrast*subjects*ROIs
%
% con_avg: contrasts*ROIs (average across subjects)
%
% con_sem: contrasts*ROIs (sem across subjects)
%
% con_names: contrast name
%
% ttest_pval:  contrasts*ROIs (t.test p.value across subjects)
%
% ROI_nm: ROI name for each ROI
%
% OUTPUTS
% roi_fig: handle for figure

%% extract general informations
n_cons = length(selectedCon);
langage = 'engl';
lSize = 3;
pSize = 50;
roi_color = ([143 0 0])./255;
% roi_color = ([143 143 143])./255;

% no more than 4 subplots per line:
nb_ROI_graph_lines = 1*(n_ROIs <= 4) + 2*(n_ROIs > 4)*(n_ROIs <= 8) + 3*(n_ROIs >8);
% if more than 4 ROIs studied, then max nber of graphs per line = 4 (more = hard to read)
max_graph_per_line = n_ROIs*(n_ROIs < 4) + 4*(n_ROIs >=4);

% default yscale values
minYvalue = min(min(min(min(con_avg(selectedCon,:) - con_sem(selectedCon,:) )))); % extract min beta value across subs and ROI
maxYvalue = max(max(max(max(con_avg(selectedCon,:) + con_sem(selectedCon,:) )))); % extract max beta value across subs and ROI
dist_min_max_Yvalue = abs(maxYvalue - minYvalue);
% define various possible default yscale values depending on minBeta and
% maxBeta
% stop scale at 0 if all contrasts are positive or negative
if (minYvalue - dist_min_max_Yvalue/3) > 0 
    minYscale = 0;
else
    minYscale = minYvalue - dist_min_max_Yvalue/3;
end
if (maxYvalue + dist_min_max_Yvalue/3) < 0 
    maxYscale = 0;
else
    maxYscale = maxYvalue + dist_min_max_Yvalue/3;
end
yscale = [minYscale,  maxYscale];
fig_size = abs(yscale(2) - yscale(1));
% position of stars for each column
pval_Ypos_pos = yscale(2) - fig_size/6;
pval_Ypos_neg = yscale(1) + fig_size/6;

%% display figure
scale_ok_idx = 0;
while scale_ok_idx == 0
    roi_fig = fig();
    
    %% loop through ROIs
    for iROI = 1:n_ROIs
        
        % create subplot containing data
        ROI_plot = subplot(nb_ROI_graph_lines, max_graph_per_line, iROI);
        
        %% loop through sequences
        xpos = (1:1:n_cons)';
        
        bar_hdl = bar(1:n_cons, con_avg(selectedCon, iROI)' );
        bar_hdl.EdgeColor = roi_color;
        bar_hdl.FaceColor = roi_color;
        hold on;
        error_hdl = errorbar(1:n_cons,...
            con_avg(selectedCon)',...
            con_sem(selectedCon)',...
            'LineStyle','none',...
            'LineWidth',lSize,'Color','k');
        
        % add stars for significant differences
        % put the p.value close to top of screen
        for iCon = 1:n_cons
            % adapt p.value position depending on number of contrasts and number of sequences
            pval_xpos = xpos(iCon);
            % disp p.value (above bar when positive contrast and below when
            % negative)
            if con_avg(selectedCon(iCon)) >= 0
                pval_Ypos = pval_Ypos_pos;
            else
                pval_Ypos = pval_Ypos_neg;
            end
            if ttest_pval(selectedCon(iCon), iROI) <= 0.05 && ttest_pval(selectedCon(iCon), iROI) > 0.01
                text(pval_xpos, pval_Ypos, '*', 'HorizontalAlignment','center','VerticalAlignment', 'top', 'FontSize', 18)
            elseif ttest_pval(selectedCon(iCon), iROI) <= 0.01 && ttest_pval(selectedCon(iCon), iROI) > 0.005
                text(pval_xpos, pval_Ypos, '**', 'HorizontalAlignment', 'center','VerticalAlignment', 'top', 'FontSize', 18)
            elseif ttest_pval(selectedCon(iCon), iROI) <= 0.005 && ttest_pval(selectedCon(iCon), iROI) > 0.001
                text(pval_xpos, pval_Ypos, '***', 'HorizontalAlignment', 'center','VerticalAlignment', 'top', 'FontSize', 18)
            elseif ttest_pval(selectedCon(iCon), iROI) <= 0.001
                text(pval_xpos, pval_Ypos, '****', 'HorizontalAlignment', 'center','VerticalAlignment', 'top', 'FontSize', 18)
            end
        end % contrast loop
        
        % draw a line to show where the zero stands
        line(xlim(),[0 0],'Color','k','LineStyle','-','LineWidth',lSize);
        
        % define graph limits
        xlim([0 n_cons+1]);
        ylim(yscale);
        
        xticks(1:n_cons);
        set(gca,'xtick',xpos,'XTickLabel',con_names);
        % ylabel('t.value');
        switch langage
            case 'engl'
                ylabel('Regression estimate');
            case 'fr'
                ylabel('coefficient de régression');
        end
        
        %
        legend_size(pSize);
        
    end % ROI loop
    
    %% check if graph scale is fine for proper display or not
    scale_fig = figure;
    set(scale_fig,'Position',[0,0,250,100]); % reduce new figure size
    scale_ok = spm_input('Y scale ok ?',1,'b','yes | no',[1 2]);
    if scale_ok == 1 % if scale is fine
        scale_ok_idx = 1;
    elseif scale_ok == 2
        yscale = spm_input('New scale values [y1 (low) y2 (high)]?',1,'e');
        fig_size = abs(yscale(2) - yscale(1));
        % position of stars for each column
        column_pval_Ypos = yscale(2) - fig_size/8;
        close(roi_fig); % close graph with previous values if not ok
    end
    close(scale_fig);
    
end % scale
    
end % function