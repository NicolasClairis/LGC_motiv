% script aimed at extracting ROI in a given GLM based on the level of a
% metabolite of interest (to select) in one of the two brain areas of study
% that can be checked (dmPFC/aIns). Then, you can select some specific
% contrast of interest and look at the difference in BOLD activity between
% the two groups

%% figure display?
fig_disp = 1;

%% define GLM number
GLM = spm_input('GLM number',1,'e');

%% define all subjects
condition = subject_condition();
[subject_id, NS] = LGCM_subject_selection('study1', condition);

%% extract metabolite levels
[low_met_subs, high_met_subs, metabolite_nm] = medSplit_metabolites(subject_id);

%% extract ROI
[con_vec_all,...
    ~, ~, ~,...
    con_names,...
    ROI_coords, ttest_ROI] = ROI_extraction_group('study1', GLM,...
    subject_id, condition, 0);
n_cons = size(con_vec_all, 1);
n_ROIs = size(con_vec_all,3);

%% extract ROI data split according to the metabolite levels
[con_avg_lowMet, con_avg_highMet,...
    con_sem_lowMet, con_sem_highMet,...
    con_std_lowMet, con_std_highMet,...
    ttest_pval_lowMet_vs_highMet, ttest_tval_lowMet_vs_highMet] = deal(NaN(n_ROIs, n_cons));
for iROI = 1:n_ROIs
    for iCon = 1:n_cons
        [con_avg_lowMet(iROI, iCon),...
            con_sem_lowMet(iROI, iCon),...
            con_std_lowMet(iROI, iCon)] = mean_sem_sd(con_vec_all(iCon, low_met_subs, iROI),2);
        [con_avg_highMet(iROI, iCon),...
            con_sem_highMet(iROI, iCon),...
            con_std_highMet(iROI, iCon)] = mean_sem_sd(con_vec_all(iCon, high_met_subs, iROI),2);
        [~,ttest_pval_lowMet_vs_highMet(iROI, iCon),~,stats_tmp] = ttest2(con_vec_all(iCon, low_met_subs, iROI), con_vec_all(iCon, high_met_subs, iROI));
        ttest_tval_lowMet_vs_highMet(iROI, iCon) = stats_tmp.tstat;
    end % contrast loop
end % ROI loop

%% how many figures do you want to plot
if fig_disp == 1
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
            % first figure: simple bar graph
            [roi_fig1] = roi_graph2(selectedContrast,...
                con_avg_lowMet(iROI,:), con_sem_lowMet(iROI,:),...
                con_avg_highMet(iROI,:), con_sem_highMet(iROI,:),...
                figConName, ttest_pval_lowMet_vs_highMet(iROI,:),...
                ['low ',metabolite_nm],['high ',metabolite_nm],...
                ROI_BOLD_nm);
            % second figure: same but with violin plots
            [roi_fig2] = roi_graph3(selectedContrast,...
                con_vec_all(:,:, iROI),...
                low_met_subs, high_met_subs,...
                figConName, ttest_pval_lowMet_vs_highMet(iROI,:),...
                ['low ',metabolite_nm],['high ',metabolite_nm],...
                ROI_BOLD_nm);
        end % ROI loop
    end % figure loop
end % figure display