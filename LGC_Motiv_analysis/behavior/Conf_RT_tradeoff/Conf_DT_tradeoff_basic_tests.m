function[] = Conf_DT_tradeoff_basic_tests(fig_disp)


%% subject selection
[study_nm, condition, gender, subject_id, NS] = sub_id;

%% working dir
computerRoot = LGCM_root_paths;

%% general parameters
mdl_nm = 'mdl_4';
nTrialsPerRun = 54;
n_totalRuns = 4;
n_runsPerTask = 2;
[]
n_bins = 6;
n_dV_bins = 3;
task_names = {'EpEm','Ep','Em'}; nTasks = length(task_names);
for iT = 1:nTasks
    task_nm = task_names{iT};
    switch task_nm
        case 'EpEm'
            n_runs0 = n_totalRuns;
        otherwise
            n_runs0 = n_runsPerTask;
    end
    [choice_f_dV.low_Conf.(task_nm), choice_f_dV.high_Conf.(task_nm),...
        DT_f_dV.low_Conf.(task_nm), DT_f_dV.high_Conf.(task_nm),...
        dV_f_dV.low_Conf.(task_nm), dV_f_dV.high_Conf.(task_nm),...
        Conf_f_dV.consistent.(task_nm), Conf_f_dV.inconsistent.(task_nm),...
        dV_f_dV.consistent.(task_nm), dV_f_dV.inconsistent.(task_nm)] = deal(NaN(n_bins, NS, n_runs0));
    
    for i_dV_bins = 1:n_dV_bins
        dV_nm = ['dV',num2str(i_dV_bins)];
        [Conf_f_DT.(dV_nm), DT_f_DT.(dV_nm)] = deal((NaN(n_bins,NS, n_runs0)));
    end % dV loop
end % task loop

%% loop over subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    subBehaviorFolder = [fullfile(computerRoot,study_nm,['CID',sub_nm],'behavior'),filesep];
    [runs, n_runs] = runs_definition(study_nm, sub_nm, condition);
    
    for iR = 1:n_runs
        run_nm = runs.runsToKeep(iR);
        run_task_nm = runs.tasks{iR};
        [task_fullName] = task_fullName_extraction(run_task_nm);
        % extract predicted confidence + subjective value
        [NV_chosen_min_unch, deltaNV_hE_min_lE, confidence_highE, pChoice,...
            deltaNV_hE_min_lE_plus_bias, NV_ch_min_unch_with_bias] = extract_bayesian_mdl(resultsFolder, subBehaviorFolder,...
            sub_nm, run_nm, task_fullName, mdl_nm);
        % extract RT
        [DT] = extract_RT(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        % extract actual confidence rating
        [conf_rtg] = extract_confidence_rating(subBehaviorFolder, sub_nm, run_nm, task_fullName);
    end
end % subject loop

%% figures
if fig_disp == 1
    
end % figure display

end % function