%% script to compare metabolites concentrations in the plasma
% vs proportion of high effort choices across individuals

%% subject selection
study_nm = 'study1';
condition = subject_condition();
[subject_id, NS] = LGCM_subject_selection(study_nm, condition, 'all');

%% load metabolites
[plasmaM, mb_names, n_mb] = load_plasma_metabolites(subject_id);

%% load proportion of choices
fig_disp = 0;
[choice_hE] = choice_hE_proportion(study_nm, condition, subject_id, fig_disp);
task_names = {'Ep','Em','EpEm'};
nTasks = length(task_names);

%% perform correlations
[pSize, lWidth, col, mSize] = general_fig_prm;
for iT = 1:nTasks
    task_nm = task_names{iT};
    
    for iM = 1:n_mb
        mb_nm = mb_names{iM};
        
        for iRaw3SD = 1:2
            switch iRaw3SD
                case 1 % raw data
                    dataType_nm = 'raw';
                    metabolite_allSubs = plasmaM.(mb_nm);
                case 2 % filtered
                    dataType_nm = 'mb_3sd_filtered';
                    metabolite_allSubs = plasmaM.filtered.(mb_nm).(mb_nm);
            end
            
            %% correlation
            GS_tmp = ~isnan(choice_hE.(task_nm).*metabolite_allSubs);
            goodSubs.(task_nm).(mb_nm).(dataType_nm) = GS_tmp;
            NS_goodSubs.(task_nm).(mb_nm).(dataType_nm) = sum(goodSubs.(task_nm).(mb_nm).(dataType_nm));
            [r_corr.(task_nm).(mb_nm).(dataType_nm),pval_corr.(task_nm).(mb_nm).(dataType_nm)] = corr(choice_hE.(task_nm)(GS_tmp)',...
                metabolite_allSubs(GS_tmp)');
            [~, betas.(task_nm).(mb_nm).(dataType_nm),...
                pval.(task_nm).(mb_nm).(dataType_nm), ~,...
                mb_sorted.(task_nm).(mb_nm).(dataType_nm),...
                choice_fit_mbSorted.(task_nm).(mb_nm).(dataType_nm)] = glm_package(metabolite_allSubs(GS_tmp)',...
                choice_hE.(task_nm)(GS_tmp)', 'normal', 'on');
            
            %% figure
            if pval_corr.(task_nm).(mb_nm).(dataType_nm) < 0.05
                fig;
                scat_hdl = scatter(metabolite_allSubs(GS_tmp)',...
                    choice_hE.(task_nm)(GS_tmp)');
                fit_hdl = plot(mb_sorted.(task_nm).(mb_nm).(dataType_nm),...
                    choice_fit_mbSorted.(task_nm).(mb_nm).(dataType_nm));
                scat_hdl_upgrade(scat_hdl);
                fit_hdl_upgrade(fit_hdl);
                xlabel(['plasma ',mb_nm,' (μM)']);
                ylabel(['Choices (%) - ',task_nm]);
                place_r_and_pval(r_corr.(task_nm).(mb_nm).(dataType_nm), pval_corr.(task_nm).(mb_nm).(dataType_nm));
                ylim([0 100]);
                legend_size(pSize);
            end % filter significant p.value
            
        end % loop over raw/mean+/-3*SD filter
    end % metabolite loop
end % task loop