% This script will allow to test whether the BOLD activity in reaction to a
% given contrast of interest in a given GLM actually correlates with the
% parameter sensitivities of the different behavioural models.

%% define subjects
condition = subject_condition();
[subject_id, NS] = LGCM_subject_selection('study1',condition);

%% define GLM number
fMRI_GLM = spm_input('fMRI GLM number',1,'e');

%% define behavioural model to check
[prm] = prm_extraction(subject_id);

%% extract list of all potential parameters
allPrm_names = fieldnames(prm);
nPrm = length(allPrm_names);

%% extract ROI contrasts
[con_vec_all,...
    ~, ~, ~,...
    con_names,...
    ROI_coords, ttest_ROI] = ROI_extraction_group('study1', fMRI_GLM,...
    subject_id, condition, 0);
n_cons = size(con_vec_all, 1);
n_ROIs = size(con_vec_all,3);
if n_ROIs > 1
    error('script not ready yet for more than 1 ROI');
end
ROI_BOLD_nm = ROI_coords.ROI_nm.ROI_1_shortName;

%% select contrast of interest
% prepare contrast names for question
list_con_question = con_names{1};
for iCon = 2:n_cons
    list_con_question = [list_con_question,' | ',con_names{iCon}];
end
% select contrast to test
which_con = zeros(1, n_cons); % column with 0/1 to see which contrasts to display at the end
% select which contrasts you want to display for each figure
selectedContrast = spm_input('What contrast to test ?',1,'m',...
    list_con_question, ...
    1:n_cons, 0);
which_con(1,selectedContrast) = 1; % 1 for selected contrasts
ROI_beta_values = NaN(1,NS);
ROI_beta_values(:) = con_vec_all(selectedContrast,:,1);
con_nm = con_names{selectedContrast};

%% test all potential correlations
pSize = 30;
lSize = 3;
for iPrm = 1:nPrm
    prm_nm = allPrm_names{iPrm};
    prm_tmp = prm.(prm_nm);
    goodSubs = ~isnan(prm_tmp);
    [betas_tmp, ~, stats_tmp] = glmfit(ROI_beta_values(goodSubs), prm_tmp(goodSubs), 'normal');
    betas.(prm_nm) = betas_tmp;
    pval.(prm_nm) = stats_tmp.p;
    fitted_prm_tmp = glmval(betas_tmp, ROI_beta_values(goodSubs), 'identity');
    
    disp([prm_nm,'=f(',ROI_BOLD_nm,' ',con_names{selectedContrast},') ;',...
        'p = ',num2str(stats_tmp.p(2))]);
    
    % display figure with correlation data
    fig;
    hold on;
    scatter(ROI_beta_values(goodSubs), prm_tmp(goodSubs),...
        'LineWidth',3);
    plot(ROI_beta_values(goodSubs), fitted_prm_tmp,...
        'LineStyle','--','LineWidth',lSize,'Color','k');
    xlabel([ROI_BOLD_nm,' ',con_nm]);
    ylabel(prm_nm);
    legend_size(pSize);
end % parameter loop