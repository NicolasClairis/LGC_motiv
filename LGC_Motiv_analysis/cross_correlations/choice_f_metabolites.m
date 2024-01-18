%% script to compare metabolites concentrations vs proportion of high effort choices across individuals

%% subject selection
study_nm = 'study1';
condition = subject_condition();
[subject_id, NS] = LGCM_subject_selection(study_nm, condition, 'all');

%% load metabolites
[metabolite_allSubs, MRS_ROI_nm, mb_nm] = metabolite_extraction(study_nm, subject_id);

%% load proportion of choices
fig_disp = 0;
[choice_hE] = choice_hE_proportion(study_nm, condition, subject_id, fig_disp);
task_names = {'Ep','Em','EpEm'};
nTasks = length(task_names);

%% perform correlations
[pSize, lWidth, col, mSize] = general_fig_prm;
for iT = 1:nTasks
    task_nm = task_names{iT};
    goodSubs.(task_nm) = ~isnan(choice_hE.(task_nm).*metabolite_allSubs);
    NS_goodSubs.(task_nm) = sum(goodSubs.(task_nm));
    [r_corr.(task_nm),pval_corr.(task_nm)] = corr(choice_hE.(task_nm)(goodSubs.(task_nm))',...
        metabolite_allSubs(goodSubs.(task_nm))');
    [~, betas.(task_nm),...
        pval.(task_nm), ~,...
        mb_sorted.(task_nm),...
        choice_fit_mbSorted.(task_nm)] = glm_package(metabolite_allSubs(goodSubs.(task_nm))',...
        choice_hE.(task_nm)(goodSubs.(task_nm))', 'normal', 'on');
    
    %% figure
    fig;
    scat_hdl = scatter(metabolite_allSubs(goodSubs.(task_nm))',...
        choice_hE.(task_nm)(goodSubs.(task_nm))');
    fit_hdl = plot(mb_sorted.(task_nm),...
        choice_fit_mbSorted.(task_nm));
    scat_hdl_upgrade(scat_hdl);
    fit_hdl_upgrade(fit_hdl);
    xlabel([MRS_ROI_nm,' - ',mb_nm,' (mM)']);
    ylabel(['Choices (%) - ',task_nm]);
    place_r_and_pval(r_corr.(task_nm), pval_corr.(task_nm));
    ylim([0 100]);
    legend_size(pSize);
end % task loop