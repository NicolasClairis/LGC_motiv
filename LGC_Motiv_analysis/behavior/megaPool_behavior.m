function[mega_mtrx] = megaPool_behavior()
%[] = megaPool_behavior()
% megaPool_behavior will pool all subjects, all sessions together to create
% one big structure pooling all the data together.
%
% OUTPUTS
% mega_mtrx: structure with all the relevant information. Each subject will
% have one subfield of mega_mtrx. Then in each subject subfield, you will
% find one field indicating the task type (physical/mental), the session
% number (1-4), the incentive of the low and the high effort option, the
% nature of each trial (reward or punishment), the high effort level, the 
% sum of accumulated effort, the previous efficiency, the choice (0/1: low/high)
% and the confidence (0/1: low/high) of the trial.
%

%% subject selection
study_nm = 'study1';
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);
% store informations
mega_mtrx.subject_selection_info.subject_id = subject_id;
mega_mtrx.subject_selection_info.NS = NS;
mega_mtrx.subject_selection_info.condition = condition;

%% working directories
rootPath = fullfile('E:',study_nm);

%% main parameters
n_runs = 4;
n_trialsPerRun = 54;
n_totalTrials = n_trialsPerRun.*n_runs;

%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    sub_fullnm = ['CID',sub_nm];
    subBehaviorFolder = [fullfile(rootPath,sub_fullnm,'behavior'), filesep];
    
    % initialize variable of interest
    mega_mtrx.(sub_fullnm).taskType = cell(n_totalTrials,1);
    [mega_mtrx.(sub_fullnm).session_nb,...
        mega_mtrx.(sub_fullnm).low_I,...
        mega_mtrx.(sub_fullnm).high_I,...
        mega_mtrx.(sub_fullnm).dI,...
        mega_mtrx.(sub_fullnm).RP,...
        mega_mtrx.(sub_fullnm).high_E,...
        mega_mtrx.(sub_fullnm).choice_hE,...
        mega_mtrx.(sub_fullnm).confidence,...
        mega_mtrx.(sub_fullnm).RT_choice,...
        mega_mtrx.(sub_fullnm).Fp,...
        mega_mtrx.(sub_fullnm).prevAccuracySpeed] = deal(NaN(n_totalTrials,1));
    
    % load which runs to use
    [runs, n_runs] = runs_definition(study_nm, sub_nm, condition);
    
    % loop through runs
    for iR = 1:n_runs
        jR = runs.runsToKeep(iR);
        trials_idx = (1:n_trialsPerRun) + n_trialsPerRun.*(jR - 1);
        task_nm_tmp = runs.tasks{iR};
        [task_fullName] = task_fullName_extraction(task_nm_tmp);
        run_nm = num2str(jR);
        
        % extract data
        [deltaI_money_tmp, low_I_money_tmp, high_I_money_tmp] = extract_monetary_incentive(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        [RP_trial_tmp] = extract_RP(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        [hE_level_tmp] = extract_hE_level(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        [choice_high_E_tmp] = extract_choice_hE(subBehaviorFolder,...
            sub_nm, run_nm, task_fullName);
        [conf_rtg_tmp] = extract_confidence_rating(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        [RT_tmp] = extract_RT(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        switch task_nm_tmp
            case 'Ep'
                [~, sumPrevAUC] = extract_physical_fatigue(subBehaviorFolder, sub_nm, run_nm, task_fullName);
            case 'Em'
                [~,...
                    ~,...
                    ~,...
                    ~,...
                    prevEfficacy_ter_with2first_tmp,...
                    ~] = extract_mental_previous_efficacy(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        end % physical/mental
        
        % load data
        mega_mtrx.(sub_fullnm).taskType(trials_idx) = repmat({task_nm_tmp},n_trialsPerRun,1);
        mega_mtrx.(sub_fullnm).session_nb(trials_idx) = jR;
        mega_mtrx.(sub_fullnm).low_I(trials_idx) = low_I_money_tmp;
        mega_mtrx.(sub_fullnm).high_I(trials_idx) = high_I_money_tmp;
        mega_mtrx.(sub_fullnm).dI(trials_idx) = deltaI_money_tmp;
        mega_mtrx.(sub_fullnm).RP(trials_idx) = RP_trial_tmp;
        mega_mtrx.(sub_fullnm).high_E(trials_idx) = hE_level_tmp;
        mega_mtrx.(sub_fullnm).choice_hE(trials_idx) = choice_high_E_tmp;
        mega_mtrx.(sub_fullnm).confidence(trials_idx) = conf_rtg_tmp;
        mega_mtrx.(sub_fullnm).RT_choice(trials_idx) = RT_tmp;
        switch task_nm_tmp
            case 'Ep'
                mega_mtrx.(sub_fullnm).Fp(trials_idx) = sumPrevAUC;
            case 'Em'
                mega_mtrx.(sub_fullnm).prevAccuracySpeed(trials_idx) = prevEfficacy_ter_with2first_tmp;
        end % physical/mental
        
    end % run loop
end % subject loop


end % function