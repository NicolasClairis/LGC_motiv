function[l_hdl, star_hdl] = add_pval_comparison(x1_vals, x2_vals, pval, x1_pos, x2_pos, NS_choice)
% [l_hdl, star_hdl] = add_pval_comparison(x1_vals, x2_vals, pval, x1_pos, x2_pos)
% add_pval_comparison will add a trait and a star above two violin or
% scatter plots which are significantly different
%
% INPUTS
% x1_vals: 1*nx1 vector with values for x1 variable 
%
% x2_vals: 1*nx2 vector with values for x2 variable
%
% pval: p.value for the comparison between x1 and x2
%
% x1_pos: x coordinate for variable x1
%
% x2_pos: x coordinate for variable x2
%
% NS_choice: indicate whether you want to keep empty if no significant
% difference or if you want to add a NS (non-significant) on top
% 'NS': add abbreviation
% '': nothing (default)
%
% OUTPUTS
% l_hdl: line handle for the trait
%
% star_hdl: star handle for the stars above the trait

%% extract size of figure in y.axis
yscale = ylim;
size_y = abs(yscale(2) - yscale(1));

%% define trait and star properties
lW = 3; % line width for trait
ftSize = 18; % star font size

%% what to do by default if NS_choice empty?
if ~exist('NS_choice','var') || isempty(NS_choice) || ~ismember(NS_choice,{'NS',''})
    NS_choice = ''; % empty by default
end

%% add trait and star if difference is significant
% if difference between x1 and x2 is significant, add a bar and a certain
% number of stars
if pval < 0.05
    
    %% define location of trait
    y_val = max( [x1_vals, x2_vals], [], 2,'omitnan') + size_y/20;
    l_hdl = line([x1_pos x2_pos],...
        [y_val y_val],...
        'LineWidth',lW,'LineStyle','-','Color','k');
    pval_xpos = mean([x1_pos x2_pos]);
    
    %% add stars on top of the trait
    pval_Ypos = y_val + size_y/20;
    
    if pval < 0.05 && pval > 0.01
        star_hdl = text(pval_xpos, pval_Ypos, '*',...
            'HorizontalAlignment','center','VerticalAlignment', 'top', 'FontSize', ftSize);
    elseif pval < 0.01 && pval > 0.001
        star_hdl = text(pval_xpos, pval_Ypos, '**',...
            'HorizontalAlignment','center','VerticalAlignment', 'top', 'FontSize', ftSize);
    elseif pval < 0.001
        star_hdl = text(pval_xpos, pval_Ypos, '***',...
            'HorizontalAlignment','center','VerticalAlignment', 'top', 'FontSize', ftSize);
    end % pval threshold
    
else % no significant difference
    switch NS_choice
        case ''
            % nothing
            l_hdl = [];
            star_hdl = [];
        case 'NS'
            %% define location of trait (same as for when difference is significant)
            y_val = max( [x1_vals, x2_vals], [], 2,'omitnan') + size_y/20;
            l_hdl = line([x1_pos x2_pos],...
                [y_val y_val],...
                'LineWidth',lW,'LineStyle','-','Color','k');
            pval_xpos = mean([x1_pos x2_pos]);
            
            %% add stars on top of the trait
            pval_Ypos = y_val + size_y/80;
            star_hdl = text(pval_xpos, pval_Ypos, 'NS',...
            'HorizontalAlignment','center','VerticalAlignment', 'bottom', 'FontSize', ftSize);
    end
    
end % significant p.value

end % function