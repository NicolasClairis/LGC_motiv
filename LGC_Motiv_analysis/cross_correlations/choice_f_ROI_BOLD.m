% This script will allow to test whether the BOLD activity in reaction to a
% given contrast of interest in a given GLM actually correlates with the
% proportion of effortful choices.

%% define subjects
study_nm = 'study1';
condition = subject_condition();
[subject_id, NS] = LGCM_subject_selection(study_nm,condition);

%% define GLM number
fMRI_GLM = spm_input('fMRI GLM number',1,'e');

%% define behavioural model to check
[choice_hE] = choice_hE_proportion(study_nm, condition, subject_id, 0);
allTask_names = fieldnames(choice_hE);
allTask_names(ismember(allTask_names,{'mean','sem'})) = [];
nTasks = length(allTask_names);

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
ROI_BOLD_short_nm1 = inputdlg('ROI BOLD contrast short name?');
ROI_BOLD_short_nm = ROI_BOLD_short_nm1{1};

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
grey = [143 143 143]./255;
subsIncluded = {'allSubs','withoutOutliers'};

for iUncorrCorr = 1:length(subsIncluded)
    uncCorr_nm = subsIncluded{iUncorrCorr};
    switch iUncorrCorr
        case 1
            disp('Uncorrected data:');
        case 2
            disp(' ');
            disp('After outlier removal:');
    end
    for iT = 1:nTasks
        task_nm = allTask_names{iT};
        choice_hE_tmp = choice_hE.(task_nm);
        switch iUncorrCorr
            case 1 % all subjects included as long as no-NaN values
                goodSubs = ~isnan(ROI_beta_values.*choice_hE_tmp);
            case 2 % remove any outlier in any measure
                [~,~,ROI_beta_values_bis] = rmv_outliers_3sd(ROI_beta_values);
                [~,~,choice_hE_tmp] = rmv_outliers_3sd(choice_hE_tmp);
                goodSubs = ~isnan(ROI_beta_values_bis.*choice_hE_tmp);
        end
        [betas_tmp, ~, stats_tmp] = glmfit(ROI_beta_values(goodSubs), choice_hE_tmp(goodSubs), 'normal');
        [rho.(uncCorr_nm).(task_nm), pval_rho.(uncCorr_nm).(task_nm)] = corr(ROI_beta_values(goodSubs)', choice_hE_tmp(goodSubs)');
        betas.(uncCorr_nm).(task_nm) = betas_tmp;
        pval.(uncCorr_nm).(task_nm) = stats_tmp.p;
        ROI_b_ascOrder = sort(ROI_beta_values(goodSubs));
        fitted_prm_tmp = glmval(betas_tmp, ROI_b_ascOrder, 'identity');
        
        disp([task_nm,'=f(',ROI_BOLD_short_nm,' ',con_names{selectedContrast},'); ',...
            'N = ',num2str(sum(goodSubs)),'; ',...
            'r = ',num2str(round(rho.(uncCorr_nm).(task_nm),3)),'; ',...
            'p = ',num2str(stats_tmp.p(2))]);
        
        % display figure with correlation data
        fig; subplot(1,3,1:2);
        hold on;
        scat_hdl = scatter(ROI_beta_values(goodSubs), choice_hE_tmp(goodSubs));
        scat_hdl_upgrade(scat_hdl);
        fit_hdl = plot(ROI_b_ascOrder, fitted_prm_tmp);
        fit_hdl_upgrade(fit_hdl);
        %         xlabel([ROI_BOLD_short_nm,' ',con_nm]);
        %         xlabel([ROI_BOLD_short_nm,' slope']);
        xlabel([ROI_BOLD_short_nm,' regression estimate']);
        %         ylabel(['Choices % ',task_nm]);
        ylabel('Choices (%)');
        %     % if you want to check the bad subject id
        %     for iS = 1:length(goodSubs)
        %         if goodSubs(iS) == 1
        %             sub_nm = subject_id{iS};
        %             text(ROI_beta_values(iS), prm_tmp(iS),...
        %                 sub_nm, 'HorizontalAlignment','center',...
        %                 'VerticalAlignment', 'top', 'FontSize', 18);
        %         end
        %     end
        [txt_hdl, txtSize] = place_r_and_pval(rho.(uncCorr_nm).(task_nm),...
            pval.(uncCorr_nm).(task_nm)(2));
        legend_size(pSize);
    end % task loop
end % uncorr/corr for outliers loop