function[roi_fig] = roi_graph(selectedCon, n_ROIs,...
        con_vec_all, con_avg, con_sem, con_names, ttest_pval)
% [roi_fig] = roi_graph(selectedCon, n_ROIs, con_vec_all, con_avg, con_sem, con_names, ttest_pval);
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
%
% Script requires violinplot function which can be downloaded here: https://github.com/bastibe/Violinplot-Matlab/


%% extract general informations
n_cons = length(selectedCon);
NS = size(con_vec_all, 2);
langage = 'engl';
lSize = 3;
pSize = 50;
% roi_color = ([143 0 0])./255;
% roi_color = ([143 143 143])./255;

% no more than 4 subplots per line:
nb_ROI_graph_lines = 1*(n_ROIs <= 4) + 2*(n_ROIs > 4)*(n_ROIs <= 8) + 3*(n_ROIs >8);
% if more than 4 ROIs studied, then max nber of graphs per line = 4 (more = hard to read)
max_graph_per_line = n_ROIs*(n_ROIs < 4) + 4*(n_ROIs >=4);

%% show mean +/- SEM or violinplot?
graphType = 'violin'; % 'violin'/'meanSem'

%% y.scale range
switch graphType
    case 'meanSem'
        % default yscale values if displaying mean +/- SEM
        minYvalue = min(min(min(min(con_avg(selectedCon,:) - con_sem(selectedCon,:) )))); % extract min beta value across subs and ROI
        maxYvalue = max(max(max(max(con_avg(selectedCon,:) + con_sem(selectedCon,:) )))); % extract max beta value across subs and ROI
    case 'violin'
        % default yscale values if displaying violinplot
        minYvalue = min(min(min(min(con_vec_all(selectedCon,:,:))))); % extract min beta value across subs and ROI
        maxYvalue = max(max(max(max(con_vec_all(selectedCon,:,:))))); % extract max beta value across subs and ROI
end
dist_min_max_Yvalue = abs(maxYvalue - minYvalue);
% define various possible default yscale values depending on minBeta and
% maxBeta
% stop scale at 0 if all contrasts are positive or negative
marginY = dist_min_max_Yvalue/8;
if (minYvalue - marginY) > 0 
    minYscale = 0;
else
    minYscale = minYvalue - marginY;
end
if (maxYvalue + marginY) < 0 
    maxYscale = 0;
else
    maxYscale = maxYvalue + marginY;
end
yscale = [minYscale,  maxYscale];

%% ratio x/y/z axis
ax_ratio = [1 1 1]; % all axis equal

%% display figure
scale_ok_idx = 0;
while scale_ok_idx == 0
    roi_fig = fig();
    
    %% loop through ROIs
    for iROI = 1:n_ROIs
        
        % create subplot containing data
        ROI_plot = subplot(nb_ROI_graph_lines, max_graph_per_line, iROI);
        pbaspect(ax_ratio);
        
        %% loop through sequences
        xpos = (1:1:n_cons)';
        
        switch graphType
            case 'meanSem'
                % show mean +/- SEM
                bar_hdl = bar(1:n_cons, con_avg(selectedCon, iROI)' );
                bar_hdl.EdgeColor = roi_color;
                bar_hdl.FaceColor = roi_color;
                hold on;
                error_hdl = errorbar(1:n_cons,...
                    con_avg(selectedCon, iROI)',...
                    con_sem(selectedCon, iROI)',...
                    'LineStyle','none',...
                    'LineWidth',lSize,'Color','k');
                
            case 'violin'
                % show violinplot
                nCon = length(selectedCon);
                con_values = NaN(NS, nCon);
                for iC = 1:nCon
                    jCon = selectedCon(iC);
                    con_values(:, iC) = con_vec_all(jCon,:, iROI);
                end
                violinplot(con_values);
        end
        
        %% add stars for significant effects
        stars_signif_display(iROI, n_cons, xpos, yscale,...
                con_avg, con_sem, con_vec_all, ttest_pval,...
                selectedCon, graphType);
        
        %% other display parameters
        % define graph limits
        xlim([0.5 n_cons+0.5]);
        ylim(yscale);

        % draw a line to show where the zero stands
        line(xlim(),[0 0],'Color','k','LineStyle','-','LineWidth',lSize);
        
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
        close(roi_fig); % close graph with previous values if not ok
    end
    close(scale_fig);
    
    %% replace stars for significant effects if scale was modified
    if scale_ok == 2
        for iROI = 1:n_ROIs
            stars_signif_display(iROI, n_cons, xpos, yscale,...
                con_avg, con_sem, con_vec_all, ttest_pval,...
                selectedCon, graphType);
        end % ROI loop
    end % scale modified
    
end % scale
    
end % function

%% function to display stars above significant effects
function[] = stars_signif_display(iROI, n_cons, xpos, yscale,...
    con_avg, con_sem, con_vec_all, ttest_pval, selectedCon, graphType)
% stars_signif_display will display the stars above the significant
% regressors.
% INPUTS: all required information to know where to place the stars (x
% location (xpos), yscale information (yscale), contrast values (con_avg, con_sem,
% con_vec_all), significant or not (ttest_pval), index of contrasts to
% check (selectedCon), number of total contrasts displayed, index of the
% current ROI (iROI) and information regarding figure type
% (graphType: meanSem= mean and SEM bars; while violin = violinplots).

% put the p.value close to top of screen
for iCon = 1:n_cons
    % adapt p.value position depending on number of contrasts and number of sequences
    pval_xpos = xpos(iCon);
    % disp p.value (above bar when positive contrast and below when
    % negative)
    switch graphType
        case 'meanSem'
            maxYvalue_tmp = con_avg(selectedCon(iCon), iROI) + con_sem(selectedCon(iCon), iROI);
            minYvalue_tmp = con_avg(selectedCon(iCon), iROI) - con_sem(selectedCon(iCon), iROI);
        case 'violin'
            maxYvalue_tmp = min([max(max(max(con_vec_all(selectedCon,:,iROI)))),yscale(2)]);
            minYvalue_tmp = max([min(min(min(con_vec_all(selectedCon,:,iROI)))),yscale(1)]);
    end
    if maxYvalue_tmp > 0
        dYpos_yscale_tmp = abs(yscale(2) - maxYvalue_tmp);
        pval_Ypos = maxYvalue_tmp + (1/3)*dYpos_yscale_tmp;
    else
        dYpos_yscale_tmp = abs(yscale(1) - minYvalue_tmp);
        pval_Ypos = minYvalue_tmp - (1/3)*dYpos_yscale_tmp;
    end
    
    if ttest_pval(selectedCon(iCon), iROI) <= 0.05
        if ttest_pval(selectedCon(iCon), iROI) <= 0.05 && ttest_pval(selectedCon(iCon), iROI) > 0.01
            text(pval_xpos, pval_Ypos, '*',...
                'HorizontalAlignment','center','VerticalAlignment', 'top', 'FontSize', 18)
        elseif ttest_pval(selectedCon(iCon), iROI) <= 0.01 && ttest_pval(selectedCon(iCon), iROI) > 0.005
            text(pval_xpos, pval_Ypos, '**',...
                'HorizontalAlignment', 'center','VerticalAlignment', 'top', 'FontSize', 18)
        elseif ttest_pval(selectedCon(iCon), iROI) <= 0.005 && ttest_pval(selectedCon(iCon), iROI) > 0.001
            text(pval_xpos, pval_Ypos, '***',...
                'HorizontalAlignment', 'center','VerticalAlignment', 'top', 'FontSize', 18)
        elseif ttest_pval(selectedCon(iCon), iROI) <= 0.001
            text(pval_xpos, pval_Ypos, '****',...
                'HorizontalAlignment', 'center','VerticalAlignment', 'top', 'FontSize', 18)
        end % how many stars for p.value?
    end % is p.value significant?
end % contrast loop
end % sub-function