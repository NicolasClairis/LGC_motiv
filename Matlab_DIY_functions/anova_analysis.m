function[pval] = anova_analysis(datamatrix, column_labels)
% anova_analysis will analyze the data entered in datamatrix with an ANOVA
% including a post-hoc Bonferroni test.
%
% INPUTS
% datamatrix: nSamples*nVar matrix with datat to analyze
%
% column_labels: cell with label for each column (important to create the
% post-hoc resulting table).
% 
% See also:
% https://sites.duke.edu/adhoc_boss/files/2013/02/Detailed-Multicompare-MATLAB-Tutorial.pdf
% for more details
% Requires the Statistics and Machine Learning matlab Toolbox

%% extract number of variables
nVars = size(datamatrix, 2);

%% perform ANOVA
[pval.global, tbl, stats] = anova1(datamatrix);
[comparison,means,h,gnames] = multcompare(stats,'alpha',0.05,'ctype','bonferroni');

%% perform post-hoc Bonferroni test
for iVar1 = 1:nVars
    column_label1 = column_labels{iVar1};
    for iVar2 = (iVar1+1):nVars
        column_label2 = column_labels{iVar2};
        comparison_nm = [column_label1,'_vs_',column_label2];
        idx_posthoc.(comparison_nm) = comparison(:,1) == iVar1 & comparison(:,2) == iVar2;
        pval.(comparison_nm) = comparison(idx_posthoc.(comparison_nm), 6);
    end % loop over variable 2
end % loop over variable 1

end % function