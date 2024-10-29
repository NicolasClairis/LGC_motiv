function[sortedContributions, sorted_var_names,...
    fig_PC1_vs_PC2, fig_PC1_vs_PC2_and_mSplit,...
    fig_cumulativeExplVariance, fig_PC_contrib] = pca_figures(n_PC_loading_details, coeff, score, explained,...
    var_names,...
    mSplit_disp, low_mSplit, high_mSplit, low_mSplit_nm, high_mSplit_nm)
% [fig_PC1_vs_PC2, fig_PC1_vs_PC2_and_mSplit,...
%     fig_cumulativeExplVariance, fig_PC_contrib] = pca_figures(n_PC_loading_details, coeff, score, explained,...
%     var_names,...
%     mSplit_disp, low_mSplit, high_mSplit, low_mSplit_nm, high_mSplit_nm)
% pca_figures will take the outputs of the pca.m function and will create
% the figures you normally want to visualize after performing a PCA.
%
% INPUTS
% n_PC_loading_details: number of PC for which you want to see the detail
% of the variables loading upon them.
%
% coeff: nVariables*nPCs PCA loadings/coefficients for each variable
% (lines) and each PC (columns)
%
% score: nObservations*nPCs principal component scores
%
% explained: nPCs*1 percentage of explained variance for each PC
%
% var_names: 1*nVariables cell with name for each variable included in the
% PCA
%
% mSplit_disp: highlight the spreading of a given variable of interest based on
% a median split
% (0) nothing
% (1) will highlight datapoints with different colours depending on if they
% belong to the low or the high group
%
% low_mSplit/high_mSplit: nObservations*1 binary vector indicating which
% subjects belongs to which group (low/high) (can be left empty if mSplit =
% 0)
%
% low_mSplit_nm/high_mSplit_nm: name for low and high median split variable
% to display as a legend on the figure
%
% OUTPUTS
% sortedContributions: nVariables*nPCs matrix containing the percentage of
% contribution of each variable for each PC
%
% sorted_var_names: corresponding variable names for sortedContributions
%
% fig_PC1_vs_PC2: figure handle of figure displaying PC1 against PC2
% (figure aimed at showing how well the PCs allow to split the data)
%
% fig_PC1_vs_PC2_and_mSplit: figure handle of the figure displaying PC1
% against PC2 (like fig_PC1_vs_PC2) but this time, the data is split
% according to one variable whose low values are shown in blue and high
% values are shown in red
%
% fig_cumulativeExplVariance: figure handle of the figure displaying the
% cumulative variance explained depending on the number of PCs considered
%
% fig_PC_contrib: structure containing the figure handle for the figures 

%% general figure parameters
pSize = 10;
n_vars = size(coeff,1);
n_PCs = size(coeff,2);

%% PC1 vs PC2
fig_PC1_vs_PC2 = fig;
scat_hdl = scatter(score(:,1), score(:,2));
scat_hdl_upgrade(scat_hdl);
xlabel('Principal Component 1');
ylabel('Principal Component 2');

%% highlight datapoints that are low/high on the graph if mSplit = 1
if mSplit_disp == 1
    fig_PC1_vs_PC2_and_mSplit = fig;
    scat_low_hdl = scatter(score(low_mSplit,1), score(low_mSplit,2));
    scat_high_hdl = scatter(score(high_mSplit,1), score(high_mSplit,2));
    scat_hdl_upgrade(scat_low_hdl);
    scat_hdl_upgrade(scat_high_hdl);
    scat_low_hdl.MarkerFaceColor = [171 217 233]./255; % blue
    scat_high_hdl.MarkerFaceColor = [215 25 28]./255; % red
    xlabel('Principal Component 1');
    ylabel('Principal Component 2');
    legend([scat_high_hdl,scat_low_hdl],{high_mSplit_nm,low_mSplit_nm});
    legend('boxoff');
else
    fig_PC1_vs_PC2_and_mSplit = [];
end % median split figure

%% show cumulative of explained variance depending on the number of PCs included
cumulativeExplained = cumsum(explained);
fig_cumulativeExplVariance = fig;
plot(cumulativeExplained, '-o');
xlabel('Number of Principal Components');
ylabel('Cumulative Explained Variance (%)');
title('Explained Variance by Principal Components');

%% Loadings on each PC
% convert loadings into percentage
squaredLoadings = coeff.^2; % Note: the sum of the squared coefficients/loadings 
% on each PC should be equal to 1 (i.e. sum(squaredLoadings,1) = 1)
percentageContribution = 100*squaredLoadings;

% general heatmap of the weights on each PC
fig;
corr_range = [-1 1];
color_range_choices = redblue(45);
subplot_hdl = subplot(1,1,1);
% display correlation matrix
imagesc(coeff, corr_range);
colormap(subplot_hdl, color_range_choices);
xlabel('Principal Components');
yticks(1:n_vars);
yticklabels(var_names);
ylabel('Variables');
% display scale
colorbar;
legend_size(pSize);

% For each PC, sort variables depending on their respective contribution to
% the corresponding PC
sortedContributions = NaN(n_vars, n_PCs);
sorted_var_names = cell(n_vars,n_PCs);
for iPC = 1:n_PCs
    % sort contributions independently for each PC in descending order
    [sortedContributions(:,iPC), sortedIndices_tmp] = sort(percentageContribution(:, iPC), 'descend');
    sorted_var_names(:,iPC) = strrep(var_names(sortedIndices_tmp),'_',' ');
end % loop over PCS

% if n_PC_loading_details not initiated => ask which value to use
if ~exist('n_PC_loading_details','var') || isempty(n_PC_loading_details)
    n_PC_loading_details_cell = inputdlg('How many PCs you want to decompose?');
    n_PC_loading_details = str2double(n_PC_loading_details_cell);
end

% plot weight for the first n_PC_loading_details PCs selected
% fig; % one figure where all plots will be displayed
for iPC = 1:n_PC_loading_details  % Loop through first selected PCs
    
    % plot corresponding results for contribution
    fig_PC_contrib.weight.(['PC',num2str(iPC)]) = fig; % one new figure for each plot
%     subplot(1, n_PC_loading_details, iPC);  % Create subplot for each PC
    bar(1:n_vars, coeff(:, iPC));  % Bar plot of sorted contributions
    xlabel('Variable');
    xticks(1:n_vars);
    xticklabels(var_names);
    ylabel('Weight');
    title(['Variable Weights on Principal Component ', num2str(iPC)]);
    legend_size(10);
end

% plot percentage contributions for the first n_PC_loading_details PCs
% selected
% fig; % one figure where all plots will be displayed
for iPC = 1:n_PC_loading_details  % Loop through first selected PCs
    
    % plot corresponding results for contribution
    fig_PC_contrib.sortedContribution.(['PC',num2str(iPC)]) = fig; % one new figure for each plot
%     subplot(1, n_PC_loading_details, iPC);  % Create subplot for each PC
    bar(sortedContributions(:, iPC));  % Bar plot of sorted contributions
    xticks(1:n_vars);
    xticklabels(sorted_var_names(:,iPC));
    xtickangle(45);
    ylabel('Percentage Contribution (%)');
    title(['Contributions to PC', num2str(iPC)]);
    legend_size(10);
end

end % function