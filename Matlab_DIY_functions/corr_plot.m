function[fig_hdl, subplot_hdl] = corr_plot(corr_mtrx, pval_mtrx,...
    corr_range, xlabels, ylabels, xlabel_nm, ylabel_nm,...
    apply_pval_threshold, pval_threshold, disp_stars,...
    fig_hdl, subplot_hdl)
% [fig_hdl, subplot_hdl] = corr_plot(corr_mtrx, pval_mtrx,...
%     corr_range, xlabels, ylabels, xlabel_nm, ylabel_nm,...
%     apply_pval_threshold, pval_threshold, disp_stars,...
%     fig_hdl, subplot_hdl)
% corr_plot will display a correlation plot with a color scale based on
% correlation coefficients and p.values entered in inputs
%
% INPUTS
% corr_mtrx: nY*nX matrix with correlation coefficients
%
% pval_mtrx: nY*nX matrix with corresponding p.values for each correlation
% test
%
% corr_range: [r1 r2] binary vector indicating the range to use for the
% color scale ([-1 1] by default if left empty)
%
% xlabels: 1*nX cell with label for each column of the matrix (ascending
% numbers by default if left empty)
%
% ylabels: 1*nY cell with label for each line of the matrix (ascending
% numbers by default if left empty)
%
% xlabel_nm, ylabel_nm: global x/y label (on top of xlabels and ylabels
% used for each tick). Can be left empty if you don't want to have a global
% label
%
% apply_pval_threshold:
% (false): display all correlations (by default if left empty)
% (true): filter to only have a display of the significant elements in
% color (all the rest will be forced to be displayed as if r=0)
%
% pval_threshold: value defining the p.value threshold if
% apply_pval_threshold is true (otherwise it is basically useless)
% pval_threshold = 0.05 by default so that any box where p>0.05 is not
% shown if apply_pval_threshold = true.
%
% disp_stars:
% (false): no stars on significant boxes
% (true): add stars depending on how significant each box is.
% Thresholds being: * for p <= 0.05, **p<0.01, ***p<0.001
% You may want to decide between pval_threshold and disp_stars what to use
% as both might be partly redundant (both serve to highlight what is significant)
%
% fig_hdl: figure handle (will be created if left empty)
%
% subplot_hdl: subplot handle (will be created if left empty)
%
% OUTPUTS
% fig_hdl: figure handle
%
% subplot_hdl: subplot handle

%% determine number of x and y elements in the matrix
n_x_vars = size(corr_mtrx, 2);
n_y_vars = size(corr_mtrx, 1);

%% default inputs
% x labels by default = numbers
if ~exist('xlabels','var')
    xlabels = cell(1,n_x_vars);
    for iX = 1:n_x_vars
        xlabels{iX} = num2str(iX);
    end
end

% y labels by default = numbers
if ~exist('ylabels','var')
    ylabels = cell(1,n_y_vars);
    for iY = 1:n_y_vars
        ylabels{iY} = num2str(iY);
    end
end

% correlation range by default
if ~exist('corr_range','var') || isempty(corr_range)
    corr_range = [-1 1];
end

% pvalue threshold default
if ~exist('apply_pval_threshold','var') || isempty(apply_pval_threshold)
    apply_pval_threshold = false;
elseif apply_pval_threshold == true
    % define default p.value threshold if apply_pval_threshold = true but
    % no p.value threshold has been applied
    if ~exist('pval_threshold','var') || isempty(pval_threshold)
        pval_threshold = 0.05;
    end
end

% display stars (by default)
if ~exist('disp_stars','var') || isempty(disp_stars)
    disp_stars = true;
end

%% mask non-significant boxes if asked in input
if apply_pval_threshold == true
    corr_mtrx(pval_mtrx > pval_threshold) = 0;
end

%% mark r=NaN as r=0 because otherwise Matlab uses the lowest corr_range value by default which gives very misleading conclusions
if sum(sum(isnan(corr_mtrx))) > 0
    corr_mtrx(isnan(corr_mtrx)) = 0;
end

%% figure parameters
[~, ~, col] = general_fig_prm;

% define which colormap you want to use (see full list here if you are not
% happy with the selection:
% https://ch.mathworks.com/help/matlab/ref/colormap.html)
% color_range_choices = 'hot';
% color_range_choices = 'turbo';
% color_range_choices = 'jet';
color_range_choices = redblue(45);

% p.value stars font size
pval_star_color = col.black;
n_totalBoxes = n_x_vars*n_y_vars;
% threshold font size to avoid making it completely invisible
pval_star_min_fontSize = 10;
pval_star_fontSize = pval_star_min_fontSize + exp(-log(n_totalBoxes/1000)); % adjust star size to the number of correlations displayed (smaller stars when boxes are bigger)

%% define figure
% create figure (if not entered in input)
if ~exist('fig_hdl','var') || isempty(fig_hdl)
    fig_hdl = fig;
end
% create subplot (otherwise display looks strange) (unless subplot already
% defined in input)
if ~exist('subplot_hdl','var') || isempty(subplot_hdl)
    subplot_hdl = subplot(1,1,1);
end
% display correlation matrix
imagesc(corr_mtrx, corr_range);
colormap(subplot_hdl, color_range_choices);

% display scale
cbar = colorbar;
cbar.Label.String = 'r';

% x legends (columns) from left to right
if ~isempty(xlabels)
    xticks(1:n_x_vars);
    xticklabels(xlabels);
end
% add global xlabel if entered in inputs (otherwise leave empty)
if exist('xlabel_nm','var') && ~isempty(xlabel_nm)
    xlabel(xlabel_nm);
end

% y legends (lines) from top to bottom
if ~isempty(ylabels)
    yticks(1:n_y_vars);
    yticklabels(ylabels);
end
% add global ylabel if entered in inputs (otherwise leave empty)
if exist('ylabel_nm','var') && ~isempty(ylabel_nm)
    ylabel(ylabel_nm);
end

%% add stars in the graph for significant correlations
if disp_stars == true
    for iX_column_var = 1:n_x_vars % loop through columns (from left to right)
        for jY_line_var = 1:n_y_vars % loop through lines (from top to bottom)
            if pval_mtrx(jY_line_var, iX_column_var) <= 0.05 % filter significant p.values
                
                % add stars if correlation is significant
                % Note: for text x and y are inverted because the X coordinate
                % for text corresponds to the columns and the Y coordinate for
                % text corresponds to the lines but reference point (0,0) is
                % top-left in both cases due to imagesc
                if pval_mtrx(jY_line_var, iX_column_var) > 0.01 && pval_mtrx(jY_line_var, iX_column_var) <= 0.05
                    pval_hdl = text(iX_column_var, jY_line_var, '*');
                elseif pval_mtrx(jY_line_var, iX_column_var) > 0.001 && pval_mtrx(jY_line_var, iX_column_var) <= 0.01
                    pval_hdl = text(iX_column_var, jY_line_var, '**');
                elseif pval_mtrx(jY_line_var, iX_column_var) <= 0.001
                    pval_hdl = text(iX_column_var, jY_line_var, '***');
                end % p.value
                
                % adjust star parameters (font size, color, etc.)
                pval_hdl.Color = pval_star_color;
                pval_hdl.FontSize = pval_star_fontSize;
                pval_hdl.FontWeight = 'bold';
                pval_hdl.HorizontalAlignment = 'center'; % center text on x-axis
                pval_hdl.VerticalAlignment = 'middle'; % center text on y-axis
            end % when p.value is significant
        end % loop over Y variables (lines)
    end % loop over X variables (columns)
end % display significant p.values on the graph

end % function