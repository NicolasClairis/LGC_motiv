% script aimed at extracting ROI in a given GLM based on the level of a
% metabolite of interest (to select) in one of the two brain areas of study
% that can be checked (dmPFC/aIns). Then, you can select some specific
% contrast of interest and look at the linear correlation between the level
% of metabolite in the given brain area and the BOLD activity.

%% figure display?
fig_disp = 1;

%% define GLM number
GLM = spm_input('GLM number',1,'e');

%% define all subjects
study_nm = 'study1';
condition = subject_condition();
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% define metabolite and ROI you want to focus on
[metabolite_allSubs, ROI_nm, metabolite_nm] = metabolite_extraction(study_nm, subject_id);

%% extract ROI from fRMI
[con_vec_all,...
    ~, ~, ~,...
    con_names,...
    ROI_coords, ttest_ROI] = ROI_extraction_group('study1', GLM,...
    subject_id, condition, 0);
n_cons = size(con_vec_all, 1);
n_ROIs = size(con_vec_all,3);

%% how many figures do you want to plot
if fig_disp == 1
    pSize = 50;
    lSize  = 3;
    n_figs = spm_input('How many figures?',1,'e');

    % prepare contrast names for question
    list_con_question = con_names{1};
    for iCon = 2:n_cons
        list_con_question = [list_con_question,' | ',con_names{iCon}];
    end

    %% select contrast of interest (ie those you would like to display in each figure)
    for iFig = 1:n_figs
        % select contrasts to display
        which_con = zeros(n_figs, n_cons); % column with 0/1 to see which contrasts to display at the end
        % select which contrasts you want to display for each figure
        selectedContrast = spm_input(['What contrast for fig.',num2str(iFig),' ?'],1,'m',...
            list_con_question, ...
            1:n_cons, 0);
        which_con(iFig,selectedContrast) = 1; % 1 for selected contrasts
        % give name to the contrast for the figure
        figConName= spm_input([con_names{selectedContrast},...
            ' short name:'],1,'s');
        
        %% figures
        for iROI = 1:n_ROIs
            ROI_BOLD_nm = ROI_coords.ROI_nm.(['ROI_',num2str(iROI),'_shortName']);
            
            %% perform the correlation
            ROI_beta_values = NaN(1,NS);
            ROI_beta_values(:) = con_vec_all(selectedContrast,:,iROI);
            good_subs = ~isnan(metabolite_allSubs);
            [betas_tmp, ~, stats_tmp] = glmfit(metabolite_allSubs(good_subs), ROI_beta_values(good_subs), 'normal');
            betas.(['MRS_',ROI_nm,'_',metabolite_nm]).(['fMRI_',ROI_BOLD_nm]) = betas_tmp;
            pval.(['MRS_',ROI_nm,'_',metabolite_nm]).(['fMRI_',ROI_BOLD_nm]) = stats_tmp.p;
            fitted_ROIbetas_tmp = glmval(betas_tmp, metabolite_allSubs(good_subs), 'identity');
            
            % display figure with correlation data
            fig;
            hold on;
            scatter(metabolite_allSubs(good_subs), ROI_beta_values(good_subs),...
                'LineWidth',3);
            plot(metabolite_allSubs(good_subs), fitted_ROIbetas_tmp,...
                'LineStyle','--','LineWidth',lSize,'Color','k');
            xlabel([ROI_nm,' ',metabolite_nm]);
            ylabel([ROI_BOLD_nm,' regression estimate']);
            title(figConName);
            legend_size(pSize);
        end % ROI loop
    end % figure loop
end % figure display