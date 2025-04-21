function[fig_hdl, circ_graph] = circ_disp(corr_mtrx, pval_mtrx, pval_thresh, var_labels)
% [fig_hdl, circ_graph] = circ_disp(corr_mtrx, pval_mtrx, pval_thresh, var_labels)
% circ_disp will display circular correlations.
%
% INPUTS
% corr_mtrx: nV*nV matrix with correlation coefficients
%
% pval_mtrx: nV*nV matrix with corresponding p.values
%
% pval_thresh: p.value threshold to use for display
%
% var_labels: nV cell with labels of variables
%
% OUTPUTS
% fig_hdl: figure handle
%
% circ_graph: graph handle
%

%% check which correlations are positive and which are negative, ideally
% to mark them with a different color but currently not implemented
% sign_corr_mtrx = (corr_mtrx >= 0) - (corr_mtrx < 0);

%% extract absolute correlation matrix to have weights (removing sign)
abs_corr_mtrx = abs(corr_mtrx);

%% remove non-significant correlations (if pval_mtrx entered in inputs)
if exist('pval_mtrx','var')
    if ~exist('pval_thresh','var') || isempty(pval_thresh) ||...
            pval_thresh < 0 || pval_thresh > 1
        pval_thresh = 0.05; % p.value threshold
    end
    abs_corr_mtrx(pval_mtrx > pval_thresh) = 0; % set r to 0 if not significant to remove
end

%% define color map
colormap = lines(length(abs_corr_mtrx));

% ideal: colour associated to sign of correlation to distinguish positive
% and negative correlations, but quite complex because current
% circularGraph script associates a color to each node and attributes the
% same color to all the connections launched from that node, instead of
% associating it to the connection => one needs to manipulate the original
% code quite substantially to make that work accordingly => left out for
% now)
% pos_col = [33 102 172];
% neg_col = [178 24 43];

%% display figure
fig_hdl = fig;
if exist('var_labels','var') && ~isempty(var_labels)
    circ_graph = circularGraph(abs_corr_mtrx,...
        'Colormap',colormap,...
        'Label',var_labels);
else
    circ_graph = circularGraph(abs_corr_mtrx,...
        'Colormap',colormap);
end

% change legend size to keep them small and visible
legend_size(10);
end % function