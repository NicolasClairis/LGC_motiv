function[matlabbatch] = First_level_loadRegressors_DCM(matlabbatch, GLMprm, study_nm, sub_nm, sub_idx, runs, n_runs,...
    subBehaviorFolder, computerRoot, n_scansPerRun, TR, DCM_mode)
% [matlabbatch] = First_level_loadRegressors_DCM(matlabbatch, GLMprm, study_nm, sub_nm, iS, runs, n_runs,...
%         subj_behavior_folder, computerRoot, n_scansPerRun, TR, DCM_mode);
% First_level_loadRegressors_DCM will load the regressors of interest for each
% task. The difference with First_level_loadRegressors is that First_level_loadRegressors_DCM
% will pool all the sessions together into one single big session while
% First_level_loadRegressors models each session independently
%
% INPUTS
% matlabbatch: structure with the First level batch
%
% GLMprm: structure with GLM parameters
%
% study_nm: study name
%
% sub_nm: subject identification number (as a string)
%
% sub_idx: index for matlabbatch
%
% runs: structure with information relative to the runs to include or not
%
% n_runs: number of runs
%
% subBehaviorFolder: subject behavioral folder name
%
% computerRoot: root of where the data is (required to be able to load net
% value if the GLM asks for it)
%
% n_scansPerRun: variable indicating the number of scans contained in
% previous runs (and henceforth the indication about the duration we have
% to move the onsets)
%
% TR: repetition time of fMRI acquisition (to multiply n_scansPerRun when
% computing the delay to add to the onsets for sessions > 1)
%
% DCM_mode:
% (1) all sessions modeled independently like in a classic univariate GLM
% => hard to manipulate for DCM but could be useful for testing
% session-specific effects or comparing sessions
% (2) sessions pooled within each task (ex: session 1 and 3 of physical
% effort will be concatenated into one single regressor) but each task will
% be modeled separately
% (3) all sessions pooled together
% (4) all trial periods are pooled together across sessions except for
% choice and effort which are modeled independently for each task (but
% pooled across sessions of the same task)
% (5) all trial periods are pooled together across sessions except for
% the effort period which is modeled independently for each task (but
% pooled across sessions of the same task)
%
% OUTPUTS
% matlabbatch: structure updated with the regressors of interest depending
% on GLMprm values
%
% See also First_level_batch_concat_DCM.m

%% is it a GLM for onsets-only?
onsets_only_GLM = GLMprm.gal.onsets_only;

%% load onsets and regressors of interest/session
[onsets, durations, regressors] = First_level_DCM_prepare_regressors_to_load(study_nm, sub_nm, runs, n_runs,...
    subBehaviorFolder, computerRoot, n_scansPerRun, TR, GLMprm);
% extract main regressors
choice_RT = regressors.choice_RT;
choice_hE = regressors.choice_hE;
choice_hE_bis = regressors.choice_hE_bis;
E_chosen = regressors.E_chosen;
E_chosen_bis = regressors.E_chosen_bis;
E_unchosen = regressors.E_unchosen;
E_varOption = regressors.E_varOption;
RP_var_binary = regressors.RP_var_binary;
NV_ch_min_unch = regressors.NV_ch_min_unch;
NV_ch_min_unch_with_bias = regressors.NV_ch_min_unch_with_bias;
pChosen = regressors.pChosen;
R_amount_varOption = regressors.R_amount_varOption;
R_level_varOption = regressors.R_level_varOption;
z_R_amount_varOption = regressors.z_R_amount_varOption;
z_R_level_varOption = regressors.z_R_level_varOption;
R_amount_chosen = regressors.R_amount_chosen;
R_level_chosen = regressors.R_level_chosen;
z_R_amount_chosen = regressors.z_R_amount_chosen;
z_R_level_chosen = regressors.z_R_level_chosen;
P_amount_varOption = regressors.P_amount_varOption;
P_level_varOption = regressors.P_level_varOption;
z_P_amount_varOption = regressors.z_P_amount_varOption;
z_P_level_varOption = regressors.z_P_level_varOption;
P_amount_chosen = regressors.P_amount_chosen;
P_level_chosen = regressors.P_level_chosen;
z_P_amount_chosen = regressors.z_P_amount_chosen;
z_P_level_chosen = regressors.z_P_level_chosen;
money_amount_left = regressors.money_amount_left;
money_level_left = regressors.money_level_left;
money_amount_right = regressors.money_amount_right;
money_level_right = regressors.money_level_right;
money_amount_chosen = regressors.money_amount_chosen;
abs_money_amount_chosen = regressors.abs_money_amount_chosen;
money_level_chosen = regressors.money_level_chosen;
abs_money_level_chosen = regressors.abs_money_level_chosen;
money_amount_unchosen = regressors.money_amount_unchosen;
abs_money_amount_unchosen = regressors.abs_money_amount_unchosen;
money_level_unchosen = regressors.money_level_unchosen;
abs_money_level_unchosen = regressors.abs_money_level_unchosen;
money_amount_varOption = regressors.money_amount_varOption;
abs_money_amount_varOption = regressors.abs_money_amount_varOption;
money_level_varOption = regressors.money_level_varOption;
abs_money_level_varOption = regressors.abs_money_level_varOption;
moneyChosen_min_moneyUnchosen_amount = regressors.moneyChosen_min_moneyUnchosen_amount;
moneyChosen_min_moneyFixed_amount = regressors.moneyChosen_min_moneyFixed_amount;
moneyChosen_min_moneyFixed_level = regressors.moneyChosen_min_moneyFixed_level;
money_amount_sum = regressors.money_amount_sum;
money_amount_obtained = regressors.money_amount_obtained;
win_vs_loss_fbk = regressors.win_vs_loss_fbk;
E_left = regressors.E_left;
E_right = regressors.E_right;
E_chosen_min_E_unchosen = regressors.E_chosen_min_E_unchosen;
Ech_min_Efixed = regressors.Ech_min_Efixed;
E_sum = regressors.E_sum;
money_level_x_E_varOption = regressors.money_level_x_E_varOption;
money_level_x_E_chosen = regressors.money_level_x_E_chosen;
R_level_x_E_varOption = regressors.R_level_x_E_varOption;
R_level_x_E_chosen = regressors.R_level_x_E_chosen;
P_level_x_E_varOption = regressors.P_level_x_E_varOption;
P_level_x_E_chosen = regressors.P_level_x_E_chosen;
NV_varOption = regressors.NV_varOption;
NV_varOption_plus_bias = regressors.NV_varOption_plus_bias;
pChoice_hE = regressors.pChoice_hE;
AUC = regressors.AUC;
AUC_overshoot = regressors.AUC_overshoot;
AUC_N = regressors.AUC_N;
AUC_overshoot_N = regressors.AUC_overshoot_N;
fatigue = regressors.fatigue;
efficacy_with2first = regressors.efficacy_with2first;
efficacy_pureNback = regressors.efficacy_pureNback;
efficacy_bis_with2first = regressors.efficacy_bis_with2first;
efficacy_bis_pureNback = regressors.efficacy_bis_pureNback;
prevEfficacy_with2first = regressors.prevEfficacy_with2first;
prevEfficacy_pureNback = regressors.prevEfficacy_pureNback;
prevEfficacy_bis_with2first = regressors.prevEfficacy_bis_with2first;
prevEfficacy_bis_pureNback = regressors.prevEfficacy_bis_pureNback;
trialN = regressors.trialN;
trialN_dEch = regressors.trialN_dEch;
trialN_dEnonDef_min_Edef = regressors.trialN_dEnonDef_min_Edef;
trialN_dEnonDef = regressors.trialN_dEnonDef;
confidence = regressors.confidence;
latency = regressors.latency;
n_correct = regressors.n_correct;
n_errors = regressors.n_errors;
RT_avg = regressors.RT_avg;
forcePeak = regressors.forcePeak;
forcePeak_N = regressors.forcePeak_N;

%% load general GLM parameters
orth_vars = GLMprm.gal.orth_vars;
zPerRun = GLMprm.gal.zPerRun;
switch zPerRun
    case 0 % raw data = no change
        raw_or_z = @(x) x;
    case 1 % zscore all variables per run
        raw_or_z = @(x) zscore(x);
end

% GLM parameters for choice Reward/Punishment pool
choice_RPpool = GLMprm.choice.(task_id).RPpool;
switch choice_RPpool
    case 0
        RPchoiceCond = {'R','P'};
    case 1
        RPchoiceCond = {'RP'};
end
chosen_RPpool = GLMprm.chosen.(task_id).RPpool;
switch chosen_RPpool
    case 0
        RPchosenCond = {'R','P'};
    case 1
        RPchosenCond = {'RP'};
end
preEcross_RPpool = GLMprm.preEffortCross.(task_id).RPpool;
switch preEcross_RPpool
    case 0
        RPpreEcrossCond = {'R','P'};
    case 1
        RPpreEcrossCond = {'RP'};
end
Eperf_RPpool = GLMprm.Eperf.(task_id).RPpool;
switch Eperf_RPpool
    case 0
        RPperfCond = {'R','P'};
    case 1
        RPperfCond = {'RP'};
end
fbk_RPpool = GLMprm.fbk.(task_id).RPpool;
switch fbk_RPpool
    case 0
        RPfbkCond = {'R','P'};
    case 1
        RPfbkCond = {'RP'};
end

% GLM parameters for choice effort pool
choice_splitPerE = GLMprm.choice.(task_id).splitPerE;
switch choice_splitPerE
    case 0
        EsplitchoiceCond = {'E'};
    case 1
        EsplitchoiceCond = {'E1','E2','E3'};
    case 2
        EsplitchoiceCond = {'Ech0','Ech1','Ech2','Ech3'};
    case 3
        EsplitchoiceCond = {'lEch','hEch'};
end
chosen_splitPerE = GLMprm.chosen.(task_id).splitPerE;
switch chosen_splitPerE
    case 0
        EsplitchosenCond = {'E'};
    case 1
        EsplitchosenCond = {'E1','E2','E3'};
    case 2
        EsplitchosenCond = {'Ech0','Ech1','Ech2','Ech3'};
    case 3
        EsplitchosenCond = {'lEch','hEch'};
end
preEcross_splitPerE = GLMprm.preEffortCross.(task_id).splitPerE;
switch preEcross_splitPerE
    case 0
        EsplitpreEcrossCond = {'E'};
    case 1
        EsplitpreEcrossCond = {'E1','E2','E3'};
    case 2
        EsplitpreEcrossCond = {'Ech0','Ech1','Ech2','Ech3'};
    case 3
        EsplitpreEcrossCond = {'lEch','hEch'};
end
Eperf_splitPerE = GLMprm.Eperf.(task_id).splitPerE;
switch Eperf_splitPerE
    case 0
        EsplitEperfCond = {'E'};
    case 1
        EsplitEperfCond = {'E1','E2','E3'};
    case 2
        EsplitEperfCond = {'Ech0','Ech1','Ech2','Ech3'};
    case 3
        EsplitEperfCond = {'lEch','hEch'};
end
fbk_splitPerE = GLMprm.fbk.(task_id).splitPerE;
switch fbk_splitPerE
    case 0
        EsplitFbkCond = {'E'};
    case 1
        EsplitFbkCond = {'E1','E2','E3'};
    case 2
        EsplitFbkCond = {'Ech0','Ech1','Ech2','Ech3'};
    case 3
        EsplitFbkCond = {'lEch','hEch'};
end

%% load matlabbatch
% initialize conditions
iCond = 0;
tasks = {'Ep','Em'};
nTasks = length(tasks);

%% all fixation cross
allCrossesModel = GLMprm.model_onset.(task_id).allCrosses;
if ismember(allCrossesModel,{'stick','boxcar'})
    
    %% adapt depending on DCM_mode
    switch DCM_mode
        case 1 % all sessions independent
            for iRun = 1:n_runs
                run_full_nm = ['run',num2str(iRun)];
                iCond = iCond + 1; % change for each run
                switch allCrossesModel
                    case 'stick'
                        modelAllCrossesdur = 0;
                    case 'boxcar'
                        modelAllCrossesdur = durations.allCrossesDur.(run_full_nm);
                end
                
                %% all cross modulators: none
                n_allCrossesMods = 0;
                allCrosses_modNames = cell(1,1);
                allCrosses_modVals = [];
                
                %% load all
                [matlabbatch] = First_level_loadEachCondition_DCM(matlabbatch, sub_idx, iCond,...
                    'fixation cross', onsets.allCrossesOnsets.(run_full_nm), modelAllCrossesdur,...
                    n_allCrossesMods, allCrosses_modNames, allCrosses_modVals,...
                    orth_vars, onsets_only_GLM);
            end % run loop
        case 2 % all tasks independent but sessions pooled
            for iTask = 1:nTasks
                task_nm = tasks{iTask};
                iCond = iCond + 1; % change for each task
                switch allCrossesModel
                    case 'stick'
                        modelAllCrossesdur = 0;
                    case 'boxcar'
                        modelAllCrossesdur = durations.allCrossesDur.(task_nm);
                end
                
                %% all cross modulators: none
                n_allCrossesMods = 0;
                allCrosses_modNames = cell(1,1);
                allCrosses_modVals = [];
                
                %% load all
                [matlabbatch] = First_level_loadEachCondition_DCM(matlabbatch, sub_idx, iCond,...
                    'fixation cross', onsets.allCrossesOnsets.(task_nm), modelAllCrossesdur,...
                    n_allCrossesMods, allCrosses_modNames, allCrosses_modVals,...
                    orth_vars, onsets_only_GLM);
            end % task loop
        case {3,4,5} % all fixation crosses pooled across sessions
            iCond = iCond + 1;
            switch allCrossesModel
                case 'stick'
                    modelAllCrossesdur = 0;
                case 'boxcar'
                    modelAllCrossesdur = durations.allCrossesDur.allTrials;
            end
            
            %% all cross modulators: none
            n_allCrossesMods = 0;
            allCrosses_modNames = cell(1,1);
            allCrosses_modVals = [];
            
            %% load all
            [matlabbatch] = First_level_loadEachCondition_DCM(matlabbatch, sub_idx, iCond,...
                'fixation cross', onsets.allCrossesOnsets.allTrials, modelAllCrossesdur,...
                n_allCrossesMods, allCrosses_modNames, allCrosses_modVals,...
                orth_vars, onsets_only_GLM);
    end % DCM mode: info about how regressors concatenated or not across tasks or sessions
end % filter if all fixation cross included or not

%% fixation cross before choice
preChoiceCrossModel = GLMprm.model_onset.(task_id).preChoiceCross;
preChoiceCrossModel_RT = GLMprm.preChoiceCross.(task_id).RT;
preChoiceCrossModel_choiceHighE = GLMprm.preChoiceCross.(task_id).choiceHighE;
preChoiceCrossModel_E_chosen    = GLMprm.preChoiceCross.(task_id).E_chosen;
if ismember(preChoiceCrossModel,{'stick','boxcar','boxcar_bis'})
    %% adapt depending on how sessions are pooled
    switch DCM_mode
        case 1 % all sessions independent
            for iRun = 1:n_runs
                run_full_nm = ['run',num2str(iRun)];
                
                iCond = iCond + 1;
                preChoiceCrossOnsets = onsets.preChoiceCrossOnsets.(run_full_nm);
                switch preChoiceCrossModel
                    case 'stick'
                        modelPreChoiceCrossdur = 0;
                    case 'boxcar'
                        modelPreChoiceCrossdur = durations.preChoiceCrossDur.(run_full_nm);
                    case 'boxcar_bis' % go from fixation cross display until choice is made
                        modelPreChoiceCrossdur = durations.preChoiceCrossDur.(run_full_nm) + choice_RT.(run_full_nm);
                end
                
                %% pre-choice cross modulators
                n_preChoiceCrossMods = 0;
                preChoiceCross_modNames = cell(1,1);
                preChoiceCross_modVals = [];
                
                % RT (first regressor)
                if preChoiceCrossModel_RT > 0 && ismember(preChoiceCrossModel_RT,[4,5,6])
                    n_preChoiceCrossMods = n_preChoiceCrossMods + 1;
                    preChoiceCross_modNames{n_preChoiceCrossMods} = 'choice RT';
                    switch preChoiceCrossModel_RT
                        case 4
                            preChoiceCross_modVals(n_preChoiceCrossMods,:) = raw_or_z(choice_RT.(run_full_nm));
                        case 5
                            preChoiceCross_modVals(n_preChoiceCrossMods,:) = zscore(choice_RT.(run_full_nm));
                        otherwise
                            error('not ready yet');
                    end
                end
                
                % choice = high effort
                switch preChoiceCrossModel_choiceHighE
                    case 0
                    case 1
                        n_preChoiceCrossMods = n_preChoiceCrossMods + 1;
                        preChoiceCross_modNames{n_preChoiceCrossMods} = 'choice = high effort';
                        preChoiceCross_modVals(n_preChoiceCrossMods,:) = choice_hE.(run_full_nm); % binary variable => no zscore
                    case 2
                        n_preChoiceCrossMods = n_preChoiceCrossMods + 1;
                        preChoiceCross_modNames{n_preChoiceCrossMods} = 'choice = high effort';
                        preChoiceCross_modVals(n_preChoiceCrossMods,:) = choice_hE_bis.(run_full_nm); % binary variable => no zscore
                    otherwise
                        error('not ready yet');
                end
                
                % effort chosen
                if preChoiceCrossModel_E_chosen > 0
                    n_preChoiceCrossMods = n_preChoiceCrossMods + 1;
                    preChoiceCross_modNames{n_preChoiceCrossMods} = 'effort chosen';
                    switch preChoiceCrossModel_E_chosen
                        case 0
                        case 1
                            preChoiceCross_modVals(n_preChoiceCrossMods,:) = raw_or_z(E_chosen.(run_full_nm));
                        case 3
                            preChoiceCross_modVals(n_preChoiceCrossMods,:) = raw_or_z(E_chosen_bis.(run_full_nm));
                        case 4
                            preChoiceCross_modVals(n_preChoiceCrossMods,:) = zscore(E_chosen.(run_full_nm));
                        otherwise
                            error('not ready yet');
                    end
                end
                
                % RT (last regressor)
                if preChoiceCrossModel_RT > 0 && ismember(preChoiceCrossModel_RT,[1,2,3])
                    n_preChoiceCrossMods = n_preChoiceCrossMods + 1;
                    preChoiceCross_modNames{n_preChoiceCrossMods} = 'choice RT';
                    switch preChoiceCrossModel_RT
                        case 1
                            preChoiceCross_modVals(n_preChoiceCrossMods,:) = raw_or_z(choice_RT.(run_full_nm));
                        case 2
                            preChoiceCross_modVals(n_preChoiceCrossMods,:) = zscore(choice_RT.(run_full_nm));
                        otherwise
                            error('not ready yet');
                    end
                end
                
                %% load all
                [matlabbatch] = First_level_loadEachCondition_DCM(matlabbatch, sub_idx, iCond,...
                    'preChoice fixation cross', preChoiceCrossOnsets, modelPreChoiceCrossdur,...
                    n_preChoiceCrossMods, preChoiceCross_modNames, preChoiceCross_modVals,...
                    orth_vars, onsets_only_GLM);
            end % run loop
        case 2 % all tasks independent but sessions pooled
            for iTask = 1:nTasks
                task_nm = tasks{iTask};
                iCond = iCond + 1;
                preChoiceCrossOnsets = onsets.preChoiceCrossOnsets.(task_nm);
                switch preChoiceCrossModel
                    case 'stick'
                        modelPreChoiceCrossdur = 0;
                    case 'boxcar'
                        modelPreChoiceCrossdur = durations.preChoiceCrossDur.(task_nm);
                    case 'boxcar_bis' % go from fixation cross display until choice is made
                        modelPreChoiceCrossdur = durations.preChoiceCrossDur.(task_nm) + choice_RT.(task_nm);
                end
                
                %% pre-choice cross modulators
                n_preChoiceCrossMods = 0;
                preChoiceCross_modNames = cell(1,1);
                preChoiceCross_modVals = [];
                
                % RT (first regressor)
                if preChoiceCrossModel_RT > 0 && ismember(preChoiceCrossModel_RT,[4,5,6])
                    n_preChoiceCrossMods = n_preChoiceCrossMods + 1;
                    preChoiceCross_modNames{n_preChoiceCrossMods} = 'choice RT';
                    switch preChoiceCrossModel_RT
                        case 4
                            preChoiceCross_modVals(n_preChoiceCrossMods,:) = raw_or_z(choice_RT.(task_nm));
                        case 5
                            preChoiceCross_modVals(n_preChoiceCrossMods,:) = zscore(choice_RT.(task_nm));
                        otherwise
                            error('not ready yet');
                    end
                end
                
                % choice = high effort
                switch preChoiceCrossModel_choiceHighE
                    case 0
                    case 1
                        n_preChoiceCrossMods = n_preChoiceCrossMods + 1;
                        preChoiceCross_modNames{n_preChoiceCrossMods} = 'choice = high effort';
                        preChoiceCross_modVals(n_preChoiceCrossMods,:) = choice_hE.(task_nm); % binary variable => no zscore
                    case 2
                        n_preChoiceCrossMods = n_preChoiceCrossMods + 1;
                        preChoiceCross_modNames{n_preChoiceCrossMods} = 'choice = high effort';
                        preChoiceCross_modVals(n_preChoiceCrossMods,:) = choice_hE_bis.(task_nm); % binary variable => no zscore
                    otherwise
                        error('not ready yet');
                end
                
                % effort chosen
                if preChoiceCrossModel_E_chosen > 0
                    n_preChoiceCrossMods = n_preChoiceCrossMods + 1;
                    preChoiceCross_modNames{n_preChoiceCrossMods} = 'effort chosen';
                    switch preChoiceCrossModel_E_chosen
                        case 0
                        case 1
                            preChoiceCross_modVals(n_preChoiceCrossMods,:) = raw_or_z(E_chosen.(task_nm));
                        case 3
                            preChoiceCross_modVals(n_preChoiceCrossMods,:) = raw_or_z(E_chosen_bis.(task_nm));
                        case 4
                            preChoiceCross_modVals(n_preChoiceCrossMods,:) = zscore(E_chosen.(task_nm));
                        otherwise
                            error('not ready yet');
                    end
                end
                
                % RT (last regressor)
                if preChoiceCrossModel_RT > 0 && ismember(preChoiceCrossModel_RT,[1,2,3])
                    n_preChoiceCrossMods = n_preChoiceCrossMods + 1;
                    preChoiceCross_modNames{n_preChoiceCrossMods} = 'choice RT';
                    switch preChoiceCrossModel_RT
                        case 1
                            preChoiceCross_modVals(n_preChoiceCrossMods,:) = raw_or_z(choice_RT.(task_nm));
                        case 2
                            preChoiceCross_modVals(n_preChoiceCrossMods,:) = zscore(choice_RT.(task_nm));
                        otherwise
                            error('not ready yet');
                    end
                end
                
                %% load all
                [matlabbatch] = First_level_loadEachCondition_DCM(matlabbatch, sub_idx, iCond,...
                    'preChoice fixation cross', preChoiceCrossOnsets, modelPreChoiceCrossdur,...
                    n_preChoiceCrossMods, preChoiceCross_modNames, preChoiceCross_modVals,...
                    orth_vars, onsets_only_GLM);
            end % task loop
        case {3,4,5} % all crosses pooled across sessions
            iCond = iCond + 1;
            preChoiceCrossOnsets = onsets.preChoiceCrossOnsets.allTrials;
            switch preChoiceCrossModel
                case 'stick'
                    modelPreChoiceCrossdur = 0;
                case 'boxcar'
                    modelPreChoiceCrossdur = preChoiceCrossDur.allTrials;
                case 'boxcar_bis' % go from fixation cross display until choice is made
                    modelPreChoiceCrossdur = preChoiceCrossDur.allTrials + choice_RT.allTrials;
            end
            
            %% pre-choice cross modulators
            n_preChoiceCrossMods = 0;
            preChoiceCross_modNames = cell(1,1);
            preChoiceCross_modVals = [];
            
            % RT (first regressor)
            if preChoiceCrossModel_RT > 0 && ismember(preChoiceCrossModel_RT,[4,5,6])
                n_preChoiceCrossMods = n_preChoiceCrossMods + 1;
                preChoiceCross_modNames{n_preChoiceCrossMods} = 'choice RT';
                switch preChoiceCrossModel_RT
                    case 4
                        preChoiceCross_modVals(n_preChoiceCrossMods,:) = raw_or_z(choice_RT.allTrials);
                    case 5
                        preChoiceCross_modVals(n_preChoiceCrossMods,:) = zscore(choice_RT.allTrials);
                    otherwise
                        error('not ready yet');
                end
            end
            
            % choice = high effort
            switch preChoiceCrossModel_choiceHighE
                case 0
                case 1
                    n_preChoiceCrossMods = n_preChoiceCrossMods + 1;
                    preChoiceCross_modNames{n_preChoiceCrossMods} = 'choice = high effort';
                    preChoiceCross_modVals(n_preChoiceCrossMods,:) = choice_hE.allTrials; % binary variable => no zscore
                case 2
                    n_preChoiceCrossMods = n_preChoiceCrossMods + 1;
                    preChoiceCross_modNames{n_preChoiceCrossMods} = 'choice = high effort';
                    preChoiceCross_modVals(n_preChoiceCrossMods,:) = choice_hE_bis.allTrials; % binary variable => no zscore
                otherwise
                    error('not ready yet');
            end
            
            % effort chosen
            if preChoiceCrossModel_E_chosen > 0
                n_preChoiceCrossMods = n_preChoiceCrossMods + 1;
                preChoiceCross_modNames{n_preChoiceCrossMods} = 'effort chosen';
                switch preChoiceCrossModel_E_chosen
                    case 0
                    case 1
                        preChoiceCross_modVals(n_preChoiceCrossMods,:) = raw_or_z(E_chosen.allTrials);
                    case 3
                        preChoiceCross_modVals(n_preChoiceCrossMods,:) = raw_or_z(E_chosen_bis.allTrials);
                    case 4
                        preChoiceCross_modVals(n_preChoiceCrossMods,:) = zscore(E_chosen.allTrials);
                    otherwise
                        error('not ready yet');
                end
            end
            
            % RT (last regressor)
            if preChoiceCrossModel_RT > 0 && ismember(preChoiceCrossModel_RT,[1,2,3])
                n_preChoiceCrossMods = n_preChoiceCrossMods + 1;
                preChoiceCross_modNames{n_preChoiceCrossMods} = 'choice RT';
                switch preChoiceCrossModel_RT
                    case 1
                        preChoiceCross_modVals(n_preChoiceCrossMods,:) = raw_or_z(choice_RT.allTrials);
                    case 2
                        preChoiceCross_modVals(n_preChoiceCrossMods,:) = zscore(choice_RT.allTrials);
                    otherwise
                        error('not ready yet');
                end
            end
            
            %% load all
            [matlabbatch] = First_level_loadEachCondition_DCM(matlabbatch, sub_idx, iCond,...
                'preChoice fixation cross', preChoiceCrossOnsets, modelPreChoiceCrossdur,...
                n_preChoiceCrossMods, preChoiceCross_modNames, preChoiceCross_modVals,...
                orth_vars, onsets_only_GLM);
    end % DCM mode: info about how regressors concatenated or not across tasks or sessions
end % filter pre-choice cross


%% choice period
choiceModel = GLMprm.model_onset.(task_id).choice;
if ismember(choiceModel,{'stick','boxcar','boxcar_bis'})
    
    for iRP_choice = 1:length(RPchoiceCond)
        RP_dispChoice_nm = RPchoiceCond{iRP_choice};
        for iEsplit_dispChoice = 1:length(EsplitchoiceCond)
            splitE_dispChoice_nm = EsplitchoiceCond{iEsplit_dispChoice};
            choiceModel_RP                              = GLMprm.choice.(task_id).(RP_dispChoice_nm).(splitE_dispChoice_nm).R_vs_P;
            choiceModel_choicehE                        = GLMprm.choice.(task_id).(RP_dispChoice_nm).(splitE_dispChoice_nm).choiceHighE;
            choiceModel_R_varOption                     = GLMprm.choice.(task_id).(RP_dispChoice_nm).(splitE_dispChoice_nm).R_varOption;
            choiceModel_R_chosen                        = GLMprm.choice.(task_id).(RP_dispChoice_nm).(splitE_dispChoice_nm).R_chosen;
            choiceModel_P_varOption                     = GLMprm.choice.(task_id).(RP_dispChoice_nm).(splitE_dispChoice_nm).P_varOption;
            choiceModel_P_chosen                        = GLMprm.choice.(task_id).(RP_dispChoice_nm).(splitE_dispChoice_nm).P_chosen;
            choiceModel_moneyLeft                       = GLMprm.choice.(task_id).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_left;
            choiceModel_moneyRight                      = GLMprm.choice.(task_id).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_right;
            choiceModel_moneyChosen                     = GLMprm.choice.(task_id).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_chosen;
            choiceModel_moneyUnchosen                   = GLMprm.choice.(task_id).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_unchosen;
            choiceModel_moneyNonDefault                 = GLMprm.choice.(task_id).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_varOption;
            choiceModel_moneyChosen_min_moneyUnchosen   = GLMprm.choice.(task_id).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_ch_min_unch;
            choiceModel_moneyChosen_min_moneyDefault    = GLMprm.choice.(task_id).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_ch_min_fixOption;
            choiceModel_moneySum                        = GLMprm.choice.(task_id).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_sum;
            choiceModel_E_left                          = GLMprm.choice.(task_id).(RP_dispChoice_nm).(splitE_dispChoice_nm).E_left;
            choiceModel_E_right                         = GLMprm.choice.(task_id).(RP_dispChoice_nm).(splitE_dispChoice_nm).E_right;
            choiceModel_E_chosen                        = GLMprm.choice.(task_id).(RP_dispChoice_nm).(splitE_dispChoice_nm).E_chosen;
            choiceModel_E_unchosen                      = GLMprm.choice.(task_id).(RP_dispChoice_nm).(splitE_dispChoice_nm).E_unchosen;
            choiceModel_E_nonDefault                    = GLMprm.choice.(task_id).(RP_dispChoice_nm).(splitE_dispChoice_nm).E_varOption;
            choiceModel_Ech_min_Eunch                   = GLMprm.choice.(task_id).(RP_dispChoice_nm).(splitE_dispChoice_nm).E_ch_min_unch;
            choiceModel_Ech_min_Efixed                  = GLMprm.choice.(task_id).(RP_dispChoice_nm).(splitE_dispChoice_nm).E_ch_min_fixOption;
            choiceModel_E_sum                           = GLMprm.choice.(task_id).(RP_dispChoice_nm).(splitE_dispChoice_nm).E_sum;
            choiceModel_money_level_x_E_varOption       = GLMprm.choice.(task_id).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_level_x_E_varOption;
            choiceModel_money_level_x_E_chosen          = GLMprm.choice.(task_id).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_level_x_E_chosen;
            choiceModel_R_level_x_E_varOption           = GLMprm.choice.(task_id).(RP_dispChoice_nm).(splitE_dispChoice_nm).R_level_x_E_varOption;
            choiceModel_R_level_x_E_chosen              = GLMprm.choice.(task_id).(RP_dispChoice_nm).(splitE_dispChoice_nm).R_level_x_E_chosen;
            choiceModel_P_level_x_E_varOption           = GLMprm.choice.(task_id).(RP_dispChoice_nm).(splitE_dispChoice_nm).P_level_x_E_varOption;
            choiceModel_P_level_x_E_chosen              = GLMprm.choice.(task_id).(RP_dispChoice_nm).(splitE_dispChoice_nm).P_level_x_E_chosen;
            choiceModel_NV_chosen                       = GLMprm.choice.(task_id).(RP_dispChoice_nm).(splitE_dispChoice_nm).NV_chosen;
            choiceModel_NV_varOption                    = GLMprm.choice.(task_id).(RP_dispChoice_nm).(splitE_dispChoice_nm).NV_varOption;
            choiceModel_NV_varOption_bis                = GLMprm.choice.(task_id).(RP_dispChoice_nm).(splitE_dispChoice_nm).NV_varOption_bis;
            switch task_id
                case 'Ep'
                    choiceModel_F_integral = GLMprm.choice.(task_id).(RP_dispChoice_nm).(splitE_dispChoice_nm).F_integral;
                    choiceModel_fatigue = GLMprm.choice.(task_id).(RP_dispChoice_nm).(splitE_dispChoice_nm).fatigue;
                    choiceModel_Ech_x_fatigue = GLMprm.choice.(task_id).(RP_dispChoice_nm).(splitE_dispChoice_nm).Ech_x_fatigue;
                case 'Em'
                    choiceModel_efficacy = GLMprm.choice.(task_id).(RP_dispChoice_nm).(splitE_dispChoice_nm).efficacy;
                    choiceModel_prevEfficacy = GLMprm.choice.(task_id).(RP_dispChoice_nm).(splitE_dispChoice_nm).prevEfficacy;
                    choiceModel_Ech_x_prevEfficacy = GLMprm.choice.(task_id).(RP_dispChoice_nm).(splitE_dispChoice_nm).Ech_x_prevEfficacy;
            end
            choiceModel_conf                            = GLMprm.choice.(task_id).(RP_dispChoice_nm).(splitE_dispChoice_nm).confidence;
            choiceModel_RT                              = GLMprm.choice.(task_id).(RP_dispChoice_nm).(splitE_dispChoice_nm).RT;
            choiceModel_trialN                          = GLMprm.choice.(task_id).(RP_dispChoice_nm).(splitE_dispChoice_nm).trialN;
            
            %% adapt depending on how sessions are pooled
            switch DCM_mode
                case 1 % all sessions independent
                    for iRun = 1:n_runs
                        run_full_nm = ['run',num2str(iRun)];
                        
                        % extract trial index for the current loop
                        switch RP_dispChoice_nm
                            case 'RP'
                                RPfilter_dispChoice = true(1,length(double(RP_var_binary.(run_full_nm))));
                            case 'R'
                                RPfilter_dispChoice = (RP_var_binary.(run_full_nm) == 1);
                            case 'P'
                                RPfilter_dispChoice = (RP_var_binary.(run_full_nm) == 0);
                        end
                        switch splitE_dispChoice_nm
                            case 'E'
                                Efilter_dispChoice = true(1,length(double(RP_var_binary.(run_full_nm))));
                            case 'E1'
                                Efilter_dispChoice = (E_varOption.(run_full_nm) == 1);
                            case 'E2'
                                Efilter_dispChoice = (E_varOption.(run_full_nm) == 2);
                            case 'E3'
                                Efilter_dispChoice = (E_varOption.(run_full_nm) == 3);
                            case 'Ech0'
                                Efilter_dispChoice = (E_chosen.(run_full_nm) == 0);
                            case 'Ech1'
                                Efilter_dispChoice = (E_chosen.(run_full_nm) == 1);
                            case 'Ech2'
                                Efilter_dispChoice = (E_chosen.(run_full_nm) == 2);
                            case 'Ech3'
                                Efilter_dispChoice = (E_chosen.(run_full_nm) == 3);
                            case 'lEch'
                                Efilter_dispChoice = (choice_hE.(run_full_nm) == 0);
                            case 'hEch'
                                Efilter_dispChoice = (choice_hE.(run_full_nm) == 1);
                        end
                        choice_trial_idx = (RPfilter_dispChoice.*Efilter_dispChoice) == 1; % NEED to transform it into logical or will just focus on the first trial
                        
                        %% choice onset
                        iCond = iCond + 1;
                        modelChoiceOnset = onsets.dispChoiceOptionOnsets.(run_full_nm)(choice_trial_idx);
                        % duration
                        switch choiceModel
                            case 'stick'
                                modelChoiceDur = 0;
                            case 'boxcar'
                                modelChoiceDur = durations.dispChoiceOptionsDur.(run_full_nm)(choice_trial_idx);
                            case 'boxcar_bis'
                                modelChoiceDur = durations.dispChoiceOptionsDur.(run_full_nm)(choice_trial_idx) + durations.dispChosenDur.(run_full_nm)(choice_trial_idx);
                        end
                        
                        %% choice modulators
                        n_choiceMods = 0;
                        choice_modNames = cell(1,1);
                        choice_modVals = [];
                        
                        % RT (first regressor)
                        if choiceModel_RT > 0 && ~ismember(choiceModel_RT,[1,2,3])
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'choice RT';
                            switch choiceModel_RT
                                case 4
                                    choice_modVals(n_choiceMods,:) = raw_or_z(choice_RT.(run_full_nm)(choice_trial_idx));
                                case 5
                                    choice_modVals(n_choiceMods,:) = zscore(choice_RT.(run_full_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % net value chosen (first regressor)
                        if choiceModel_NV_chosen > 0 && ~ismember(choiceModel_NV_chosen,[1,2,3])
                            n_choiceMods = n_choiceMods + 1;
                            switch choiceModel_NV_chosen
                                case 4
                                    choice_modNames{n_choiceMods} = 'NVch-NVunch';
                                    choice_modVals(n_choiceMods,:) = raw_or_z(NV_ch_min_unch.(run_full_nm)(choice_trial_idx));
                                case 5
                                    choice_modNames{n_choiceMods} = 'p(chosen)';
                                    choice_modVals(n_choiceMods,:) = raw_or_z(pChosen.(run_full_nm)(choice_trial_idx));
                                case 6
                                    choice_modNames{n_choiceMods} = 'NVch-NVunch';
                                    choice_modVals(n_choiceMods,:) = raw_or_z(NV_ch_min_unch_with_bias.(run_full_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % reward vs punishments
                        if choiceModel_RP == 1
                            if strcmp(RP_dispChoice_nm,'RP')
                                n_choiceMods = n_choiceMods + 1;
                                choice_modNames{n_choiceMods} = 'R vs P';
                                choice_modVals(n_choiceMods,:) = RP_var_binary.(run_full_nm)(choice_trial_idx); % binary variable => no zscore
                            else
                                error('cannot split R and P trials and add a variable representing R/P trials') ;
                            end
                        end
                        
                        % choice = high effort
                        switch choiceModel_choicehE
                            case 0
                            case 1
                                n_choiceMods = n_choiceMods + 1;
                                choice_modNames{n_choiceMods} = 'choice = high effort';
                                choice_modVals(n_choiceMods,:) = choice_hE.(run_full_nm)(choice_trial_idx); % binary variable => no zscore
                            case 2
                                n_choiceMods = n_choiceMods + 1;
                                choice_modNames{n_choiceMods} = 'choice = high effort';
                                choice_modVals(n_choiceMods,:) = choice_hE_bis.(run_full_nm)(choice_trial_idx); % binary variable => no zscore
                            otherwise
                                error('case not ready yet');
                        end
                        
                        % reward for the non-default option
                        if choiceModel_R_varOption > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'R non-default';
                            switch choiceModel_R_varOption
                                case 1 % R high effort option (amount)
                                    choice_modVals(n_choiceMods,:) = raw_or_z(R_amount_varOption.(run_full_nm)(choice_trial_idx));
                                case 2 % R high effort option (level)
                                    choice_modVals(n_choiceMods,:) = raw_or_z(R_level_varOption.(run_full_nm)(choice_trial_idx));
                                case 3 % z(R high effort option) (amount)
                                    choice_modVals(n_choiceMods,:) = raw_or_z(z_R_amount_varOption.(run_full_nm)(choice_trial_idx));
                                case 4 % z(R high effort option) (level)
                                    choice_modVals(n_choiceMods,:) = raw_or_z(z_R_level_varOption.(run_full_nm)(choice_trial_idx));
                            end
                        end
                        
                        % reward for the chosen option
                        if choiceModel_R_chosen > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'R chosen';
                            switch choiceModel_R_chosen
                                case 1 % R chosen option amount
                                    choice_modVals(n_choiceMods,:) = raw_or_z(R_amount_chosen.(run_full_nm)(choice_trial_idx));
                                case 2 % R chosen option level (0/1/2/3)
                                    choice_modVals(n_choiceMods,:) = raw_or_z(R_level_chosen.(run_full_nm)(choice_trial_idx));
                                case 3 % z(R chosen option) amount
                                    choice_modVals(n_choiceMods,:) = raw_or_z(z_R_amount_chosen.(run_full_nm)(choice_trial_idx));
                                case 4 % z(R chosen option) level (0/1/2/3)
                                    choice_modVals(n_choiceMods,:) = raw_or_z(z_R_level_chosen.(run_full_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % punishment for the non-default option
                        if choiceModel_P_varOption > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'P non-default';
                            switch choiceModel_P_varOption
                                case 1 % R high effort option (amount)
                                    choice_modVals(n_choiceMods,:) = raw_or_z(P_amount_varOption.(run_full_nm)(choice_trial_idx));
                                case 2 % R high effort option (level)
                                    choice_modVals(n_choiceMods,:) = raw_or_z(P_level_varOption.(run_full_nm)(choice_trial_idx));
                                case 3 % z(R high effort option) (amount)
                                    choice_modVals(n_choiceMods,:) = raw_or_z(z_P_amount_varOption.(run_full_nm)(choice_trial_idx));
                                case 4 % z(R high effort option) (level)
                                    choice_modVals(n_choiceMods,:) = raw_or_z(z_P_level_varOption.(run_full_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % punishment for the chosen option
                        if choiceModel_P_chosen > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'P chosen';
                            switch choiceModel_P_chosen
                                case 1 % R high effort option (amount)
                                    choice_modVals(n_choiceMods,:) = raw_or_z(P_amount_chosen.(run_full_nm)(choice_trial_idx));
                                case 2 % R high effort option (level)
                                    choice_modVals(n_choiceMods,:) = raw_or_z(P_level_chosen.(run_full_nm)(choice_trial_idx));
                                case 3 % z(R high effort option) (amount)
                                    choice_modVals(n_choiceMods,:) = raw_or_z(z_P_amount_chosen.(run_full_nm)(choice_trial_idx));
                                case 4 % z(R high effort option) (level)
                                    choice_modVals(n_choiceMods,:) = raw_or_z(z_P_level_chosen.(run_full_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % money left
                        if choiceModel_moneyLeft > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'money left';
                            switch choiceModel_moneyLeft
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(money_amount_left.(run_full_nm)(choice_trial_idx));
                                case 2
                                    choice_modVals(n_choiceMods,:) = raw_or_z(abs(money_amount_left.(run_full_nm)(choice_trial_idx)));
                                case 3
                                    choice_modVals(n_choiceMods,:) = raw_or_z(money_level_left.(run_full_nm)(choice_trial_idx));
                                case 4
                                    choice_modVals(n_choiceMods,:) = raw_or_z(abs(money_level_left.(run_full_nm)(choice_trial_idx)));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % money right
                        if choiceModel_moneyRight > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'money right';
                            switch choiceModel_moneyRight
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(money_amount_right.(run_full_nm)(choice_trial_idx));
                                case 2
                                    choice_modVals(n_choiceMods,:) = raw_or_z(abs(money_amount_right.(run_full_nm)(choice_trial_idx)));
                                case 3
                                    choice_modVals(n_choiceMods,:) = raw_or_z(money_level_right.(run_full_nm)(choice_trial_idx));
                                case 4
                                    choice_modVals(n_choiceMods,:) = raw_or_z(abs(money_level_right.(run_full_nm)(choice_trial_idx)));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % money chosen
                        if choiceModel_moneyChosen > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'money chosen';
                            switch choiceModel_moneyChosen
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(money_amount_chosen.(run_full_nm)(choice_trial_idx));
                                case 2
                                    choice_modVals(n_choiceMods,:) = raw_or_z(abs_money_amount_chosen.(run_full_nm)(choice_trial_idx));
                                case 3
                                    choice_modVals(n_choiceMods,:) = raw_or_z(money_level_chosen.(run_full_nm)(choice_trial_idx));
                                case 4
                                    choice_modVals(n_choiceMods,:) = raw_or_z(abs_money_level_chosen.(run_full_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % money unchosen
                        if choiceModel_moneyUnchosen > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'money unchosen';
                            switch choiceModel_moneyUnchosen
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(money_amount_unchosen.(run_full_nm)(choice_trial_idx));
                                case 2
                                    choice_modVals(n_choiceMods,:) = raw_or_z(abs_money_amount_unchosen.(run_full_nm)(choice_trial_idx));
                                case 3
                                    choice_modVals(n_choiceMods,:) = raw_or_z(money_level_unchosen.(run_full_nm)(choice_trial_idx));
                                case 4
                                    choice_modVals(n_choiceMods,:) = raw_or_z(abs_money_level_unchosen.(run_full_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % money associated to the non-default option
                        if choiceModel_moneyNonDefault > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'money non-default';
                            switch choiceModel_moneyNonDefault
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(money_amount_varOption.(run_full_nm)(choice_trial_idx));
                                case 2
                                    choice_modVals(n_choiceMods,:) = raw_or_z(abs_money_amount_varOption.(run_full_nm)(choice_trial_idx));
                                case 3
                                    choice_modVals(n_choiceMods,:) = raw_or_z(money_level_varOption.(run_full_nm)(choice_trial_idx));
                                case 4
                                    choice_modVals(n_choiceMods,:) = raw_or_z(abs_money_level_varOption.(run_full_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % money chosen - money unchosen
                        if choiceModel_moneyChosen_min_moneyUnchosen > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'money chosen - unchosen';
                            switch choiceModel_moneyChosen_min_moneyUnchosen
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(moneyChosen_min_moneyUnchosen_amount.(run_full_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % money chosen - money fixed option (low R low E)
                        if choiceModel_moneyChosen_min_moneyDefault > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'money chosen - money fixed';
                            switch choiceModel_moneyChosen_min_moneyDefault
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(moneyChosen_min_moneyFixed_amount.(run_full_nm)(choice_trial_idx));
                                case 2
                                    choice_modVals(n_choiceMods,:) = raw_or_z(moneyChosen_min_moneyFixed_level.(run_full_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % sum of money
                        if choiceModel_moneySum > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'money sum';
                            switch choiceModel_moneySum
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(money_amount_sum.(run_full_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % effort left
                        if choiceModel_E_left > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'E left';
                            switch choiceModel_E_left
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(E_left.(run_full_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % effort right
                        if choiceModel_E_right > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'E right';
                            switch choiceModel_E_right
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(E_right.(run_full_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % effort chosen
                        if choiceModel_E_chosen > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'effort chosen';
                            switch choiceModel_E_chosen
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(E_chosen.(run_full_nm)(choice_trial_idx));
                                case 3
                                    choice_modVals(n_choiceMods,:) = raw_or_z(E_chosen_bis.(run_full_nm)(choice_trial_idx));
                                case 4
                                    choice_modVals(n_choiceMods,:) = zscore(E_chosen.(run_full_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % effort unchosen
                        if choiceModel_E_unchosen > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'effort unchosen';
                            switch choiceModel_E_unchosen
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(E_unchosen.(run_full_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % effort non-default option
                        if choiceModel_E_nonDefault > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'effort non-default option';
                            switch choiceModel_E_nonDefault
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(E_varOption.(run_full_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % effort chosen - unchosen
                        if choiceModel_Ech_min_Eunch > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'effort chosen - unchosen';
                            switch choiceModel_Ech_min_Eunch
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(E_chosen_min_E_unchosen.(run_full_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % effort chosen - fixed option
                        if choiceModel_Ech_min_Efixed > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'effort chosen - fixed';
                            switch choiceModel_Ech_min_Efixed
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(Ech_min_Efixed.(run_full_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % sum of efforts
                        if choiceModel_E_sum > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'effort sum';
                            switch choiceModel_E_sum
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(E_sum.(run_full_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % money level*effort level high E option
                        if choiceModel_money_level_x_E_varOption > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'money x effort non-default option';
                            switch choiceModel_money_level_x_E_varOption
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(money_level_x_E_varOption.(run_full_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % money level * effort chosen
                        if choiceModel_money_level_x_E_chosen > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'money x effort chosen';
                            switch choiceModel_money_level_x_E_chosen
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(money_level_x_E_chosen.(run_full_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % R level * effort high E option
                        if choiceModel_R_level_x_E_varOption > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'R x E non-default option';
                            switch choiceModel_R_level_x_E_varOption
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(R_level_x_E_varOption.(run_full_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % R level * effort chosen
                        if choiceModel_R_level_x_E_chosen > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'R x E chosen';
                            switch choiceModel_R_level_x_E_chosen
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(R_level_x_E_chosen.(run_full_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % P level * effort high E option
                        if choiceModel_P_level_x_E_varOption > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'P x E non-default option';
                            switch choiceModel_P_level_x_E_varOption
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(P_level_x_E_varOption.(run_full_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % P level * effort chosen
                        if choiceModel_P_level_x_E_chosen > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'P x E chosen';
                            switch choiceModel_P_level_x_E_chosen
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(P_level_x_E_chosen.(run_full_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % net value chosen
                        if choiceModel_NV_chosen > 0 && ~ismember(choiceModel_NV_chosen,[4,5,6])
                            n_choiceMods = n_choiceMods + 1;
                            switch choiceModel_NV_chosen
                                case 1
                                    choice_modNames{n_choiceMods} = 'NVch-NVunch';
                                    choice_modVals(n_choiceMods,:) = raw_or_z(NV_ch_min_unch.(run_full_nm)(choice_trial_idx));
                                case 2
                                    choice_modNames{n_choiceMods} = 'p(chosen)';
                                    choice_modVals(n_choiceMods,:) = raw_or_z(pChosen.(run_full_nm)(choice_trial_idx));
                                case 3
                                    choice_modNames{n_choiceMods} = 'NVch-NVunch';
                                    choice_modVals(n_choiceMods,:) = raw_or_z(NV_ch_min_unch_with_bias.(run_full_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % net value non-default option
                        if choiceModel_NV_varOption > 0
                            n_choiceMods = n_choiceMods + 1;
                            switch choiceModel_NV_varOption
                                case 1
                                    choice_modNames{n_choiceMods} = 'delta NV high E - low E';
                                    choice_modVals(n_choiceMods,:) = raw_or_z(NV_varOption.(run_full_nm)(choice_trial_idx));
                                case 2
                                    choice_modNames{n_choiceMods} = '|delta NV high E - low E|';
                                    choice_modVals(n_choiceMods,:) = raw_or_z(abs(NV_varOption.(run_full_nm)(choice_trial_idx)));
                                case 3
                                    choice_modNames{n_choiceMods} = 'p(choice=hE)';
                                    choice_modVals(n_choiceMods,:) = raw_or_z(pChoice_hE.(run_full_nm)(choice_trial_idx));
                                case 4
                                    choice_modNames{n_choiceMods} = 'delta NV high E - low E + bias';
                                    choice_modVals(n_choiceMods,:) = raw_or_z(NV_varOption_plus_bias.(run_full_nm)(choice_trial_idx));
                                case 5
                                    choice_modNames{n_choiceMods} = '|delta NV high E - low E + bias|';
                                    choice_modVals(n_choiceMods,:) = raw_or_z(abs(NV_varOption_plus_bias.(run_full_nm)(choice_trial_idx)));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % net value non-default option bis
                        if choiceModel_NV_varOption_bis > 0
                            n_choiceMods = n_choiceMods + 1;
                            switch choiceModel_NV_varOption_bis
                                case 1
                                    choice_modNames{n_choiceMods} = 'delta NV high E - low E';
                                    choice_modVals(n_choiceMods,:) = raw_or_z(NV_varOption.(run_full_nm)(choice_trial_idx));
                                case 2
                                    choice_modNames{n_choiceMods} = '|delta NV high E - low E|';
                                    choice_modVals(n_choiceMods,:) = raw_or_z(abs(NV_varOption.(run_full_nm)(choice_trial_idx)));
                                case 3
                                    choice_modNames{n_choiceMods} = 'p(choice=hE)';
                                    choice_modVals(n_choiceMods,:) = raw_or_z(pChoice_hE.(run_full_nm)(choice_trial_idx));
                                case 4
                                    choice_modNames{n_choiceMods} = 'delta NV high E - low E + bias';
                                    choice_modVals(n_choiceMods,:) = raw_or_z(NV_varOption_plus_bias.(run_full_nm)(choice_trial_idx));
                                case 5
                                    choice_modNames{n_choiceMods} = '|delta NV high E - low E + bias|';
                                    choice_modVals(n_choiceMods,:) = raw_or_z(abs(NV_varOption_plus_bias.(run_full_nm)(choice_trial_idx)));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        if strcmp(task_id,'Ep')
                            % force integral
                            if choiceModel_F_integral > 0
                                n_choiceMods = n_choiceMods + 1;
                                choice_modNames{n_choiceMods} = 'F integral';
                                switch choiceModel_F_integral
                                    case 1
                                        choice_modVals(n_choiceMods,:) = raw_or_z(AUC.(run_full_nm)(choice_trial_idx));
                                    case 2
                                        choice_modVals(n_choiceMods,:) = raw_or_z(AUC_overshoot.(run_full_nm)(choice_trial_idx));
                                    case 3
                                        choice_modVals(n_choiceMods,:) = raw_or_z(AUC_N.(run_full_nm)(choice_trial_idx));
                                    case 4
                                        choice_modVals(n_choiceMods,:) = raw_or_z(AUC_overshoot_N.(run_full_nm)(choice_trial_idx));
                                    case 5
                                        choice_modVals(n_choiceMods,:) = zscore(AUC.(run_full_nm)(choice_trial_idx));
                                    case 6
                                        choice_modVals(n_choiceMods,:) = zscore(AUC_overshoot.(run_full_nm)(choice_trial_idx));
                                    case 7
                                        choice_modVals(n_choiceMods,:) = zscore(AUC_N.(run_full_nm)(choice_trial_idx));
                                    case 8
                                        choice_modVals(n_choiceMods,:) = zscore(AUC_overshoot_N.(run_full_nm)(choice_trial_idx));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                            
                            % fatigue
                            if choiceModel_fatigue > 0
                                n_choiceMods = n_choiceMods + 1;
                                choice_modNames{n_choiceMods} = 'fatigue';
                                switch choiceModel_fatigue
                                    case 1
                                        choice_modVals(n_choiceMods,:) = raw_or_z(fatigue.(run_full_nm)(choice_trial_idx));
                                    case 2
                                        choice_modVals(n_choiceMods,:) = zscore(fatigue.(run_full_nm)(choice_trial_idx));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                            
                            % Ech*fatigue
                            if choiceModel_Ech_x_fatigue > 0
                                n_choiceMods = n_choiceMods + 1;
                                choice_modNames{n_choiceMods} = 'Ech_x_fatigue';
                                switch choiceModel_Ech_x_fatigue
                                    case 1
                                        choice_modVals(n_choiceMods,:) = raw_or_z(E_chosen.(run_full_nm)(choice_trial_idx).*fatigue.(run_full_nm)(choice_trial_idx));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                        end % physical effort filter
                        
                        if strcmp(task_id,'Em')
                            % efficacy
                            if choiceModel_efficacy > 0
                                n_choiceMods = n_choiceMods + 1;
                                choice_modNames{n_choiceMods} = 'efficacy';
                                switch choiceModel_efficacy
                                    case 1
                                        choice_modVals(n_choiceMods,:) = raw_or_z(efficacy_with2first.(run_full_nm)(choice_trial_idx));
                                    case 2
                                        choice_modVals(n_choiceMods,:) = raw_or_z(efficacy_pureNback.(run_full_nm)(choice_trial_idx));
                                    case 3
                                        choice_modVals(n_choiceMods,:) = raw_or_z(efficacy_bis_with2first.(run_full_nm)(choice_trial_idx));
                                    case 4
                                        choice_modVals(n_choiceMods,:) = raw_or_z(efficacy_bis_pureNback.(run_full_nm)(choice_trial_idx));
                                    case 5
                                        choice_modVals(n_choiceMods,:) = zscore(efficacy_with2first.(run_full_nm)(choice_trial_idx));
                                    case 6
                                        choice_modVals(n_choiceMods,:) = zscore(efficacy_pureNback.(run_full_nm)(choice_trial_idx));
                                    case 7
                                        choice_modVals(n_choiceMods,:) = zscore(efficacy_bis_with2first.(run_full_nm)(choice_trial_idx));
                                    case 8
                                        choice_modVals(n_choiceMods,:) = zscore(efficacy_bis_pureNback.(run_full_nm)(choice_trial_idx));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                            
                            % previous trial efficacy
                            if choiceModel_prevEfficacy > 0
                                n_choiceMods = n_choiceMods + 1;
                                choice_modNames{n_choiceMods} = 'previous trial efficacy';
                                switch choiceModel_prevEfficacy
                                    case 1
                                        choice_modVals(n_choiceMods,:) = raw_or_z(prevEfficacy_with2first.(run_full_nm)(choice_trial_idx));
                                    case 2
                                        choice_modVals(n_choiceMods,:) = raw_or_z(prevEfficacy_pureNback.(run_full_nm)(choice_trial_idx));
                                    case 3
                                        choice_modVals(n_choiceMods,:) = raw_or_z(prevEfficacy_bis_with2first.(run_full_nm)(choice_trial_idx));
                                    case 4
                                        choice_modVals(n_choiceMods,:) = raw_or_z(prevEfficacy_bis_pureNback.(run_full_nm)(choice_trial_idx));
                                    case 5
                                        choice_modVals(n_choiceMods,:) = zscore(prevEfficacy_with2first.(run_full_nm)(choice_trial_idx));
                                    case 6
                                        choice_modVals(n_choiceMods,:) = zscore(prevEfficacy_pureNback.(run_full_nm)(choice_trial_idx));
                                    case 7
                                        choice_modVals(n_choiceMods,:) = zscore(prevEfficacy_bis_with2first.(run_full_nm)(choice_trial_idx));
                                    case 8
                                        choice_modVals(n_choiceMods,:) = zscore(prevEfficacy_bis_pureNback.(run_full_nm)(choice_trial_idx));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                            
                            % Ech*(previous trial efficacy)
                            if choiceModel_Ech_x_prevEfficacy > 0
                                n_choiceMods = n_choiceMods + 1;
                                choice_modNames{n_choiceMods} = 'Ech_x_previous trial efficacy';
                                switch choiceModel_Ech_x_prevEfficacy
                                    case 1
                                        choice_modVals(n_choiceMods,:) = raw_or_z(E_chosen.(run_full_nm)(choice_trial_idx).*prevEfficacy_with2first.(run_full_nm)(choice_trial_idx));
                                    case 2
                                        choice_modVals(n_choiceMods,:) = raw_or_z(E_chosen.(run_full_nm)(choice_trial_idx).*prevEfficacy_pureNback.(run_full_nm)(choice_trial_idx));
                                    case 3
                                        choice_modVals(n_choiceMods,:) = raw_or_z(E_chosen.(run_full_nm)(choice_trial_idx).*prevEfficacy_bis_with2first.(run_full_nm)(choice_trial_idx));
                                    case 4
                                        choice_modVals(n_choiceMods,:) = raw_or_z(E_chosen.(run_full_nm)(choice_trial_idx).*prevEfficacy_bis_pureNback.(run_full_nm)(choice_trial_idx));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                        end % mental effort filter
                        
                        % trial number
                        if choiceModel_trialN > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'trial number';
                            switch choiceModel_trialN
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(trialN.(run_full_nm)(choice_trial_idx));
                                case 2
                                    choice_modVals(n_choiceMods,:) = raw_or_z(trialN_dEch.(run_full_nm)(choice_trial_idx));
                                case 3
                                    choice_modVals(n_choiceMods,:) = raw_or_z(trialN_dEnonDef_min_Edef.(run_full_nm)(choice_trial_idx));
                                case 4
                                    choice_modVals(n_choiceMods,:) = raw_or_z(trialN_dEnonDef.(run_full_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % choice confidence
                        if choiceModel_conf > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'confidence';
                            switch choiceModel_conf
                                case 1 % binary variable => no zscore
                                    choice_modVals(n_choiceMods,:) = confidence.(run_full_nm)(choice_trial_idx);
                                case {2,3,4} % confidence inferred by the model => ok to zscore
                                    choice_modVals(n_choiceMods,:) = raw_or_z(confidence.(run_full_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % RT (last regressor)
                        if choiceModel_RT > 0 && ~ismember(choiceModel_RT,[4,5,6])
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'choice RT';
                            switch choiceModel_RT
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(choice_RT.(run_full_nm)(choice_trial_idx));
                                case 2
                                    choice_modVals(n_choiceMods,:) = zscore(choice_RT.(run_full_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        [matlabbatch] = First_level_loadEachCondition_DCM(matlabbatch, sub_idx, iCond,...
                            ['choice_',RP_dispChoice_nm,'_',splitE_dispChoice_nm], modelChoiceOnset, modelChoiceDur,...
                            n_choiceMods, choice_modNames, choice_modVals,...
                            orth_vars, onsets_only_GLM);
                    end % run loop
                case {2,4} % all tasks independent but sessions pooled
                    for iTask = 1:nTasks
                        task_nm = tasks{iTask};
                        
                        % extract trial index for the current loop
                        switch RP_dispChoice_nm
                            case 'RP'
                                RPfilter_dispChoice = true(1,length(double(RP_var_binary.(task_nm))));
                            case 'R'
                                RPfilter_dispChoice = (RP_var_binary.(task_nm) == 1);
                            case 'P'
                                RPfilter_dispChoice = (RP_var_binary.(task_nm) == 0);
                        end
                        switch splitE_dispChoice_nm
                            case 'E'
                                Efilter_dispChoice = true(1,length(double(RP_var_binary.(task_nm))));
                            case 'E1'
                                Efilter_dispChoice = (E_varOption.(task_nm) == 1);
                            case 'E2'
                                Efilter_dispChoice = (E_varOption.(task_nm) == 2);
                            case 'E3'
                                Efilter_dispChoice = (E_varOption.(task_nm) == 3);
                            case 'Ech0'
                                Efilter_dispChoice = (E_chosen.(task_nm) == 0);
                            case 'Ech1'
                                Efilter_dispChoice = (E_chosen.(task_nm) == 1);
                            case 'Ech2'
                                Efilter_dispChoice = (E_chosen.(task_nm) == 2);
                            case 'Ech3'
                                Efilter_dispChoice = (E_chosen.(task_nm) == 3);
                            case 'lEch'
                                Efilter_dispChoice = (choice_hE.(task_nm) == 0);
                            case 'hEch'
                                Efilter_dispChoice = (choice_hE.(task_nm) == 1);
                        end
                        choice_trial_idx = (RPfilter_dispChoice.*Efilter_dispChoice) == 1; % NEED to transform it into logical or will just focus on the first trial
                        
                        %% choice onset
                        iCond = iCond + 1;
                        modelChoiceOnset = onsets.dispChoiceOptionOnsets.(task_nm)(choice_trial_idx);
                        % duration
                        switch choiceModel
                            case 'stick'
                                modelChoiceDur = 0;
                            case 'boxcar'
                                modelChoiceDur = durations.dispChoiceOptionsDur.(task_nm)(choice_trial_idx);
                            case 'boxcar_bis'
                                modelChoiceDur = durations.dispChoiceOptionsDur.(task_nm)(choice_trial_idx) + durations.dispChosenDur.(task_nm)(choice_trial_idx);
                        end
                        
                        %% choice modulators
                        n_choiceMods = 0;
                        choice_modNames = cell(1,1);
                        choice_modVals = [];
                        
                        % RT (first regressor)
                        if choiceModel_RT > 0 && ~ismember(choiceModel_RT,[1,2,3])
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'choice RT';
                            switch choiceModel_RT
                                case 4
                                    choice_modVals(n_choiceMods,:) = raw_or_z(choice_RT.(task_nm)(choice_trial_idx));
                                case 5
                                    choice_modVals(n_choiceMods,:) = zscore(choice_RT.(task_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % net value chosen (first regressor)
                        if choiceModel_NV_chosen > 0 && ~ismember(choiceModel_NV_chosen,[1,2,3])
                            n_choiceMods = n_choiceMods + 1;
                            switch choiceModel_NV_chosen
                                case 4
                                    choice_modNames{n_choiceMods} = 'NVch-NVunch';
                                    choice_modVals(n_choiceMods,:) = raw_or_z(NV_ch_min_unch.(task_nm)(choice_trial_idx));
                                case 5
                                    choice_modNames{n_choiceMods} = 'p(chosen)';
                                    choice_modVals(n_choiceMods,:) = raw_or_z(pChosen.(task_nm)(choice_trial_idx));
                                case 6
                                    choice_modNames{n_choiceMods} = 'NVch-NVunch';
                                    choice_modVals(n_choiceMods,:) = raw_or_z(NV_ch_min_unch_with_bias.(task_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % reward vs punishments
                        if choiceModel_RP == 1
                            if strcmp(RP_dispChoice_nm,'RP')
                                n_choiceMods = n_choiceMods + 1;
                                choice_modNames{n_choiceMods} = 'R vs P';
                                choice_modVals(n_choiceMods,:) = RP_var_binary.(task_nm)(choice_trial_idx); % binary variable => no zscore
                            else
                                error('cannot split R and P trials and add a variable representing R/P trials') ;
                            end
                        end
                        
                        % choice = high effort
                        switch choiceModel_choicehE
                            case 0
                            case 1
                                n_choiceMods = n_choiceMods + 1;
                                choice_modNames{n_choiceMods} = 'choice = high effort';
                                choice_modVals(n_choiceMods,:) = choice_hE.(task_nm)(choice_trial_idx); % binary variable => no zscore
                            case 2
                                n_choiceMods = n_choiceMods + 1;
                                choice_modNames{n_choiceMods} = 'choice = high effort';
                                choice_modVals(n_choiceMods,:) = choice_hE_bis.(task_nm)(choice_trial_idx); % binary variable => no zscore
                            otherwise
                                error('case not ready yet');
                        end
                        
                        % reward for the non-default option
                        if choiceModel_R_varOption > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'R non-default';
                            switch choiceModel_R_varOption
                                case 1 % R high effort option (amount)
                                    choice_modVals(n_choiceMods,:) = raw_or_z(R_amount_varOption.(task_nm)(choice_trial_idx));
                                case 2 % R high effort option (level)
                                    choice_modVals(n_choiceMods,:) = raw_or_z(R_level_varOption.(task_nm)(choice_trial_idx));
                                case 3 % z(R high effort option) (amount)
                                    choice_modVals(n_choiceMods,:) = raw_or_z(z_R_amount_varOption.(task_nm)(choice_trial_idx));
                                case 4 % z(R high effort option) (level)
                                    choice_modVals(n_choiceMods,:) = raw_or_z(z_R_level_varOption.(task_nm)(choice_trial_idx));
                            end
                        end
                        
                        % reward for the chosen option
                        if choiceModel_R_chosen > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'R chosen';
                            switch choiceModel_R_chosen
                                case 1 % R chosen option amount
                                    choice_modVals(n_choiceMods,:) = raw_or_z(R_amount_chosen.(task_nm)(choice_trial_idx));
                                case 2 % R chosen option level (0/1/2/3)
                                    choice_modVals(n_choiceMods,:) = raw_or_z(R_level_chosen.(task_nm)(choice_trial_idx));
                                case 3 % z(R chosen option) amount
                                    choice_modVals(n_choiceMods,:) = raw_or_z(z_R_amount_chosen.(task_nm)(choice_trial_idx));
                                case 4 % z(R chosen option) level (0/1/2/3)
                                    choice_modVals(n_choiceMods,:) = raw_or_z(z_R_level_chosen.(task_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % punishment for the non-default option
                        if choiceModel_P_varOption > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'P non-default';
                            switch choiceModel_P_varOption
                                case 1 % R high effort option (amount)
                                    choice_modVals(n_choiceMods,:) = raw_or_z(P_amount_varOption.(task_nm)(choice_trial_idx));
                                case 2 % R high effort option (level)
                                    choice_modVals(n_choiceMods,:) = raw_or_z(P_level_varOption.(task_nm)(choice_trial_idx));
                                case 3 % z(R high effort option) (amount)
                                    choice_modVals(n_choiceMods,:) = raw_or_z(z_P_amount_varOption.(task_nm)(choice_trial_idx));
                                case 4 % z(R high effort option) (level)
                                    choice_modVals(n_choiceMods,:) = raw_or_z(z_P_level_varOption.(task_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % punishment for the chosen option
                        if choiceModel_P_chosen > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'P chosen';
                            switch choiceModel_P_chosen
                                case 1 % R high effort option (amount)
                                    choice_modVals(n_choiceMods,:) = raw_or_z(P_amount_chosen.(task_nm)(choice_trial_idx));
                                case 2 % R high effort option (level)
                                    choice_modVals(n_choiceMods,:) = raw_or_z(P_level_chosen.(task_nm)(choice_trial_idx));
                                case 3 % z(R high effort option) (amount)
                                    choice_modVals(n_choiceMods,:) = raw_or_z(z_P_amount_chosen.(task_nm)(choice_trial_idx));
                                case 4 % z(R high effort option) (level)
                                    choice_modVals(n_choiceMods,:) = raw_or_z(z_P_level_chosen.(task_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % money left
                        if choiceModel_moneyLeft > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'money left';
                            switch choiceModel_moneyLeft
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(money_amount_left.(task_nm)(choice_trial_idx));
                                case 2
                                    choice_modVals(n_choiceMods,:) = raw_or_z(abs(money_amount_left.(task_nm)(choice_trial_idx)));
                                case 3
                                    choice_modVals(n_choiceMods,:) = raw_or_z(money_level_left.(task_nm)(choice_trial_idx));
                                case 4
                                    choice_modVals(n_choiceMods,:) = raw_or_z(abs(money_level_left.(task_nm)(choice_trial_idx)));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % money right
                        if choiceModel_moneyRight > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'money right';
                            switch choiceModel_moneyRight
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(money_amount_right.(task_nm)(choice_trial_idx));
                                case 2
                                    choice_modVals(n_choiceMods,:) = raw_or_z(abs(money_amount_right.(task_nm)(choice_trial_idx)));
                                case 3
                                    choice_modVals(n_choiceMods,:) = raw_or_z(money_level_right.(task_nm)(choice_trial_idx));
                                case 4
                                    choice_modVals(n_choiceMods,:) = raw_or_z(abs(money_level_right.(task_nm)(choice_trial_idx)));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % money chosen
                        if choiceModel_moneyChosen > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'money chosen';
                            switch choiceModel_moneyChosen
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(money_amount_chosen.(task_nm)(choice_trial_idx));
                                case 2
                                    choice_modVals(n_choiceMods,:) = raw_or_z(abs_money_amount_chosen.(task_nm)(choice_trial_idx));
                                case 3
                                    choice_modVals(n_choiceMods,:) = raw_or_z(money_level_chosen.(task_nm)(choice_trial_idx));
                                case 4
                                    choice_modVals(n_choiceMods,:) = raw_or_z(abs_money_level_chosen.(task_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % money unchosen
                        if choiceModel_moneyUnchosen > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'money unchosen';
                            switch choiceModel_moneyUnchosen
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(money_amount_unchosen.(task_nm)(choice_trial_idx));
                                case 2
                                    choice_modVals(n_choiceMods,:) = raw_or_z(abs_money_amount_unchosen.(task_nm)(choice_trial_idx));
                                case 3
                                    choice_modVals(n_choiceMods,:) = raw_or_z(money_level_unchosen.(task_nm)(choice_trial_idx));
                                case 4
                                    choice_modVals(n_choiceMods,:) = raw_or_z(abs_money_level_unchosen.(task_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % money associated to the non-default option
                        if choiceModel_moneyNonDefault > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'money non-default';
                            switch choiceModel_moneyNonDefault
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(money_amount_varOption.(task_nm)(choice_trial_idx));
                                case 2
                                    choice_modVals(n_choiceMods,:) = raw_or_z(abs_money_amount_varOption.(task_nm)(choice_trial_idx));
                                case 3
                                    choice_modVals(n_choiceMods,:) = raw_or_z(money_level_varOption.(task_nm)(choice_trial_idx));
                                case 4
                                    choice_modVals(n_choiceMods,:) = raw_or_z(abs_money_level_varOption.(task_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % money chosen - money unchosen
                        if choiceModel_moneyChosen_min_moneyUnchosen > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'money chosen - unchosen';
                            switch choiceModel_moneyChosen_min_moneyUnchosen
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(moneyChosen_min_moneyUnchosen_amount.(task_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % money chosen - money fixed option (low R low E)
                        if choiceModel_moneyChosen_min_moneyDefault > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'money chosen - money fixed';
                            switch choiceModel_moneyChosen_min_moneyDefault
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(moneyChosen_min_moneyFixed_amount.(task_nm)(choice_trial_idx));
                                case 2
                                    choice_modVals(n_choiceMods,:) = raw_or_z(moneyChosen_min_moneyFixed_level.(task_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % sum of money
                        if choiceModel_moneySum > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'money sum';
                            switch choiceModel_moneySum
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(money_amount_sum.(task_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % effort left
                        if choiceModel_E_left > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'E left';
                            switch choiceModel_E_left
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(E_left.(task_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % effort right
                        if choiceModel_E_right > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'E right';
                            switch choiceModel_E_right
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(E_right.(task_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % effort chosen
                        if choiceModel_E_chosen > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'effort chosen';
                            switch choiceModel_E_chosen
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(E_chosen.(task_nm)(choice_trial_idx));
                                case 3
                                    choice_modVals(n_choiceMods,:) = raw_or_z(E_chosen_bis.(task_nm)(choice_trial_idx));
                                case 4
                                    choice_modVals(n_choiceMods,:) = zscore(E_chosen.(task_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % effort unchosen
                        if choiceModel_E_unchosen > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'effort unchosen';
                            switch choiceModel_E_unchosen
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(E_unchosen.(task_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % effort non-default option
                        if choiceModel_E_nonDefault > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'effort non-default option';
                            switch choiceModel_E_nonDefault
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(E_varOption.(task_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % effort chosen - unchosen
                        if choiceModel_Ech_min_Eunch > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'effort chosen - unchosen';
                            switch choiceModel_Ech_min_Eunch
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(E_chosen_min_E_unchosen.(task_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % effort chosen - fixed option
                        if choiceModel_Ech_min_Efixed > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'effort chosen - fixed';
                            switch choiceModel_Ech_min_Efixed
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(Ech_min_Efixed.(task_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % sum of efforts
                        if choiceModel_E_sum > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'effort sum';
                            switch choiceModel_E_sum
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(E_sum.(task_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % money level*effort level high E option
                        if choiceModel_money_level_x_E_varOption > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'money x effort non-default option';
                            switch choiceModel_money_level_x_E_varOption
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(money_level_x_E_varOption.(task_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % money level * effort chosen
                        if choiceModel_money_level_x_E_chosen > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'money x effort chosen';
                            switch choiceModel_money_level_x_E_chosen
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(money_level_x_E_chosen.(task_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % R level * effort high E option
                        if choiceModel_R_level_x_E_varOption > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'R x E non-default option';
                            switch choiceModel_R_level_x_E_varOption
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(R_level_x_E_varOption.(task_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % R level * effort chosen
                        if choiceModel_R_level_x_E_chosen > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'R x E chosen';
                            switch choiceModel_R_level_x_E_chosen
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(R_level_x_E_chosen.(task_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % P level * effort high E option
                        if choiceModel_P_level_x_E_varOption > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'P x E non-default option';
                            switch choiceModel_P_level_x_E_varOption
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(P_level_x_E_varOption.(task_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % P level * effort chosen
                        if choiceModel_P_level_x_E_chosen > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'P x E chosen';
                            switch choiceModel_P_level_x_E_chosen
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(P_level_x_E_chosen.(task_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % net value chosen
                        if choiceModel_NV_chosen > 0 && ~ismember(choiceModel_NV_chosen,[4,5,6])
                            n_choiceMods = n_choiceMods + 1;
                            switch choiceModel_NV_chosen
                                case 1
                                    choice_modNames{n_choiceMods} = 'NVch-NVunch';
                                    choice_modVals(n_choiceMods,:) = raw_or_z(NV_ch_min_unch.(task_nm)(choice_trial_idx));
                                case 2
                                    choice_modNames{n_choiceMods} = 'p(chosen)';
                                    choice_modVals(n_choiceMods,:) = raw_or_z(pChosen.(task_nm)(choice_trial_idx));
                                case 3
                                    choice_modNames{n_choiceMods} = 'NVch-NVunch';
                                    choice_modVals(n_choiceMods,:) = raw_or_z(NV_ch_min_unch_with_bias.(task_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % net value non-default option
                        if choiceModel_NV_varOption > 0
                            n_choiceMods = n_choiceMods + 1;
                            switch choiceModel_NV_varOption
                                case 1
                                    choice_modNames{n_choiceMods} = 'delta NV high E - low E';
                                    choice_modVals(n_choiceMods,:) = raw_or_z(NV_varOption.(task_nm)(choice_trial_idx));
                                case 2
                                    choice_modNames{n_choiceMods} = '|delta NV high E - low E|';
                                    choice_modVals(n_choiceMods,:) = raw_or_z(abs(NV_varOption.(task_nm)(choice_trial_idx)));
                                case 3
                                    choice_modNames{n_choiceMods} = 'p(choice=hE)';
                                    choice_modVals(n_choiceMods,:) = raw_or_z(pChoice_hE.(task_nm)(choice_trial_idx));
                                case 4
                                    choice_modNames{n_choiceMods} = 'delta NV high E - low E + bias';
                                    choice_modVals(n_choiceMods,:) = raw_or_z(NV_varOption_plus_bias.(task_nm)(choice_trial_idx));
                                case 5
                                    choice_modNames{n_choiceMods} = '|delta NV high E - low E + bias|';
                                    choice_modVals(n_choiceMods,:) = raw_or_z(abs(NV_varOption_plus_bias.(task_nm)(choice_trial_idx)));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % net value non-default option bis
                        if choiceModel_NV_varOption_bis > 0
                            n_choiceMods = n_choiceMods + 1;
                            switch choiceModel_NV_varOption_bis
                                case 1
                                    choice_modNames{n_choiceMods} = 'delta NV high E - low E';
                                    choice_modVals(n_choiceMods,:) = raw_or_z(NV_varOption.(task_nm)(choice_trial_idx));
                                case 2
                                    choice_modNames{n_choiceMods} = '|delta NV high E - low E|';
                                    choice_modVals(n_choiceMods,:) = raw_or_z(abs(NV_varOption.(task_nm)(choice_trial_idx)));
                                case 3
                                    choice_modNames{n_choiceMods} = 'p(choice=hE)';
                                    choice_modVals(n_choiceMods,:) = raw_or_z(pChoice_hE.(task_nm)(choice_trial_idx));
                                case 4
                                    choice_modNames{n_choiceMods} = 'delta NV high E - low E + bias';
                                    choice_modVals(n_choiceMods,:) = raw_or_z(NV_varOption_plus_bias.(task_nm)(choice_trial_idx));
                                case 5
                                    choice_modNames{n_choiceMods} = '|delta NV high E - low E + bias|';
                                    choice_modVals(n_choiceMods,:) = raw_or_z(abs(NV_varOption_plus_bias.(task_nm)(choice_trial_idx)));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        if strcmp(task_id,'Ep')
                            % force integral
                            if choiceModel_F_integral > 0
                                n_choiceMods = n_choiceMods + 1;
                                choice_modNames{n_choiceMods} = 'F integral';
                                switch choiceModel_F_integral
                                    case 1
                                        choice_modVals(n_choiceMods,:) = raw_or_z(AUC.(task_nm)(choice_trial_idx));
                                    case 2
                                        choice_modVals(n_choiceMods,:) = raw_or_z(AUC_overshoot.(task_nm)(choice_trial_idx));
                                    case 3
                                        choice_modVals(n_choiceMods,:) = raw_or_z(AUC_N.(task_nm)(choice_trial_idx));
                                    case 4
                                        choice_modVals(n_choiceMods,:) = raw_or_z(AUC_overshoot_N.(task_nm)(choice_trial_idx));
                                    case 5
                                        choice_modVals(n_choiceMods,:) = zscore(AUC.(task_nm)(choice_trial_idx));
                                    case 6
                                        choice_modVals(n_choiceMods,:) = zscore(AUC_overshoot.(task_nm)(choice_trial_idx));
                                    case 7
                                        choice_modVals(n_choiceMods,:) = zscore(AUC_N.(task_nm)(choice_trial_idx));
                                    case 8
                                        choice_modVals(n_choiceMods,:) = zscore(AUC_overshoot_N.(task_nm)(choice_trial_idx));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                            
                            % fatigue
                            if choiceModel_fatigue > 0
                                n_choiceMods = n_choiceMods + 1;
                                choice_modNames{n_choiceMods} = 'fatigue';
                                switch choiceModel_fatigue
                                    case 1
                                        choice_modVals(n_choiceMods,:) = raw_or_z(fatigue.(task_nm)(choice_trial_idx));
                                    case 2
                                        choice_modVals(n_choiceMods,:) = zscore(fatigue.(task_nm)(choice_trial_idx));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                            
                            % Ech*fatigue
                            if choiceModel_Ech_x_fatigue > 0
                                n_choiceMods = n_choiceMods + 1;
                                choice_modNames{n_choiceMods} = 'Ech_x_fatigue';
                                switch choiceModel_Ech_x_fatigue
                                    case 1
                                        choice_modVals(n_choiceMods,:) = raw_or_z(E_chosen.(task_nm)(choice_trial_idx).*fatigue.(task_nm)(choice_trial_idx));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                        end % physical effort filter
                        
                        if strcmp(task_id,'Em')
                            % efficacy
                            if choiceModel_efficacy > 0
                                n_choiceMods = n_choiceMods + 1;
                                choice_modNames{n_choiceMods} = 'efficacy';
                                switch choiceModel_efficacy
                                    case 1
                                        choice_modVals(n_choiceMods,:) = raw_or_z(efficacy_with2first.(task_nm)(choice_trial_idx));
                                    case 2
                                        choice_modVals(n_choiceMods,:) = raw_or_z(efficacy_pureNback.(task_nm)(choice_trial_idx));
                                    case 3
                                        choice_modVals(n_choiceMods,:) = raw_or_z(efficacy_bis_with2first.(task_nm)(choice_trial_idx));
                                    case 4
                                        choice_modVals(n_choiceMods,:) = raw_or_z(efficacy_bis_pureNback.(task_nm)(choice_trial_idx));
                                    case 5
                                        choice_modVals(n_choiceMods,:) = zscore(efficacy_with2first.(task_nm)(choice_trial_idx));
                                    case 6
                                        choice_modVals(n_choiceMods,:) = zscore(efficacy_pureNback.(task_nm)(choice_trial_idx));
                                    case 7
                                        choice_modVals(n_choiceMods,:) = zscore(efficacy_bis_with2first.(task_nm)(choice_trial_idx));
                                    case 8
                                        choice_modVals(n_choiceMods,:) = zscore(efficacy_bis_pureNback.(task_nm)(choice_trial_idx));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                            
                            % previous trial efficacy
                            if choiceModel_prevEfficacy > 0
                                n_choiceMods = n_choiceMods + 1;
                                choice_modNames{n_choiceMods} = 'previous trial efficacy';
                                switch choiceModel_prevEfficacy
                                    case 1
                                        choice_modVals(n_choiceMods,:) = raw_or_z(prevEfficacy_with2first.(task_nm)(choice_trial_idx));
                                    case 2
                                        choice_modVals(n_choiceMods,:) = raw_or_z(prevEfficacy_pureNback.(task_nm)(choice_trial_idx));
                                    case 3
                                        choice_modVals(n_choiceMods,:) = raw_or_z(prevEfficacy_bis_with2first.(task_nm)(choice_trial_idx));
                                    case 4
                                        choice_modVals(n_choiceMods,:) = raw_or_z(prevEfficacy_bis_pureNback.(task_nm)(choice_trial_idx));
                                    case 5
                                        choice_modVals(n_choiceMods,:) = zscore(prevEfficacy_with2first.(task_nm)(choice_trial_idx));
                                    case 6
                                        choice_modVals(n_choiceMods,:) = zscore(prevEfficacy_pureNback.(task_nm)(choice_trial_idx));
                                    case 7
                                        choice_modVals(n_choiceMods,:) = zscore(prevEfficacy_bis_with2first.(task_nm)(choice_trial_idx));
                                    case 8
                                        choice_modVals(n_choiceMods,:) = zscore(prevEfficacy_bis_pureNback.(task_nm)(choice_trial_idx));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                            
                            % Ech*(previous trial efficacy)
                            if choiceModel_Ech_x_prevEfficacy > 0
                                n_choiceMods = n_choiceMods + 1;
                                choice_modNames{n_choiceMods} = 'Ech_x_previous trial efficacy';
                                switch choiceModel_Ech_x_prevEfficacy
                                    case 1
                                        choice_modVals(n_choiceMods,:) = raw_or_z(E_chosen.(task_nm)(choice_trial_idx).*prevEfficacy_with2first.(task_nm)(choice_trial_idx));
                                    case 2
                                        choice_modVals(n_choiceMods,:) = raw_or_z(E_chosen.(task_nm)(choice_trial_idx).*prevEfficacy_pureNback.(task_nm)(choice_trial_idx));
                                    case 3
                                        choice_modVals(n_choiceMods,:) = raw_or_z(E_chosen.(task_nm)(choice_trial_idx).*prevEfficacy_bis_with2first.(task_nm)(choice_trial_idx));
                                    case 4
                                        choice_modVals(n_choiceMods,:) = raw_or_z(E_chosen.(task_nm)(choice_trial_idx).*prevEfficacy_bis_pureNback.(task_nm)(choice_trial_idx));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                        end % mental effort filter
                        
                        % trial number
                        if choiceModel_trialN > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'trial number';
                            switch choiceModel_trialN
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(trialN.(task_nm)(choice_trial_idx));
                                case 2
                                    choice_modVals(n_choiceMods,:) = raw_or_z(trialN_dEch.(task_nm)(choice_trial_idx));
                                case 3
                                    choice_modVals(n_choiceMods,:) = raw_or_z(trialN_dEnonDef_min_Edef.(task_nm)(choice_trial_idx));
                                case 4
                                    choice_modVals(n_choiceMods,:) = raw_or_z(trialN_dEnonDef.(task_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % choice confidence
                        if choiceModel_conf > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'confidence';
                            switch choiceModel_conf
                                case 1 % binary variable => no zscore
                                    choice_modVals(n_choiceMods,:) = confidence.(task_nm)(choice_trial_idx);
                                case {2,3,4} % confidence inferred by the model => ok to zscore
                                    choice_modVals(n_choiceMods,:) = raw_or_z(confidence.(task_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % RT (last regressor)
                        if choiceModel_RT > 0 && ~ismember(choiceModel_RT,[4,5,6])
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'choice RT';
                            switch choiceModel_RT
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(choice_RT.(task_nm)(choice_trial_idx));
                                case 2
                                    choice_modVals(n_choiceMods,:) = zscore(choice_RT.(task_nm)(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        [matlabbatch] = First_level_loadEachCondition_DCM(matlabbatch, sub_idx, iCond,...
                            ['choice_',RP_dispChoice_nm,'_',splitE_dispChoice_nm], modelChoiceOnset, modelChoiceDur,...
                            n_choiceMods, choice_modNames, choice_modVals,...
                            orth_vars, onsets_only_GLM);
                    end % task loop
                case {3,5} % all pooled across sessions
                    
                    % extract trial index for the current loop
                    switch RP_dispChoice_nm
                        case 'RP'
                            RPfilter_dispChoice = true(1,length(double(RP_var_binary.allTrials)));
                        case 'R'
                            RPfilter_dispChoice = (RP_var_binary.allTrials == 1);
                        case 'P'
                            RPfilter_dispChoice = (RP_var_binary.allTrials == 0);
                    end
                    switch splitE_dispChoice_nm
                        case 'E'
                            Efilter_dispChoice = true(1,length(double(RP_var_binary.allTrials)));
                        case 'E1'
                            Efilter_dispChoice = (E_varOption.allTrials == 1);
                        case 'E2'
                            Efilter_dispChoice = (E_varOption.allTrials == 2);
                        case 'E3'
                            Efilter_dispChoice = (E_varOption.allTrials == 3);
                        case 'Ech0'
                            Efilter_dispChoice = (E_chosen.allTrials == 0);
                        case 'Ech1'
                            Efilter_dispChoice = (E_chosen.allTrials == 1);
                        case 'Ech2'
                            Efilter_dispChoice = (E_chosen.allTrials == 2);
                        case 'Ech3'
                            Efilter_dispChoice = (E_chosen.allTrials == 3);
                        case 'lEch'
                            Efilter_dispChoice = (choice_hE.allTrials == 0);
                        case 'hEch'
                            Efilter_dispChoice = (choice_hE.allTrials == 1);
                    end
                    choice_trial_idx = (RPfilter_dispChoice.*Efilter_dispChoice) == 1; % NEED to transform it into logical or will just focus on the first trial
                    
                    %% choice onset
                    iCond = iCond + 1;
                    modelChoiceOnset = onsets.dispChoiceOptionOnsets.allTrials(choice_trial_idx);
                    % duration
                    switch choiceModel
                        case 'stick'
                            modelChoiceDur = 0;
                        case 'boxcar'
                            modelChoiceDur = durations.dispChoiceOptionsDur.allTrials(choice_trial_idx);
                        case 'boxcar_bis'
                            modelChoiceDur = durations.dispChoiceOptionsDur.allTrials(choice_trial_idx) + durations.dispChosenDur.allTrials(choice_trial_idx);
                    end
                    
                    %% choice modulators
                    n_choiceMods = 0;
                    choice_modNames = cell(1,1);
                    choice_modVals = [];
                    
                    % RT (first regressor)
                    if choiceModel_RT > 0 && ~ismember(choiceModel_RT,[1,2,3])
                        n_choiceMods = n_choiceMods + 1;
                        choice_modNames{n_choiceMods} = 'choice RT';
                        switch choiceModel_RT
                            case 4
                                choice_modVals(n_choiceMods,:) = raw_or_z(choice_RT.allTrials(choice_trial_idx));
                            case 5
                                choice_modVals(n_choiceMods,:) = zscore(choice_RT.allTrials(choice_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % net value chosen (first regressor)
                    if choiceModel_NV_chosen > 0 && ~ismember(choiceModel_NV_chosen,[1,2,3])
                        n_choiceMods = n_choiceMods + 1;
                        switch choiceModel_NV_chosen
                            case 4
                                choice_modNames{n_choiceMods} = 'NVch-NVunch';
                                choice_modVals(n_choiceMods,:) = raw_or_z(NV_ch_min_unch.allTrials(choice_trial_idx));
                            case 5
                                choice_modNames{n_choiceMods} = 'p(chosen)';
                                choice_modVals(n_choiceMods,:) = raw_or_z(pChosen.allTrials(choice_trial_idx));
                            case 6
                                choice_modNames{n_choiceMods} = 'NVch-NVunch';
                                choice_modVals(n_choiceMods,:) = raw_or_z(NV_ch_min_unch_with_bias.allTrials(choice_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % reward vs punishments
                    if choiceModel_RP == 1
                        if strcmp(RP_dispChoice_nm,'RP')
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'R vs P';
                            choice_modVals(n_choiceMods,:) = RP_var_binary.allTrials(choice_trial_idx); % binary variable => no zscore
                        else
                            error('cannot split R and P trials and add a variable representing R/P trials') ;
                        end
                    end
                    
                    % choice = high effort
                    switch choiceModel_choicehE
                        case 0
                        case 1
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'choice = high effort';
                            choice_modVals(n_choiceMods,:) = choice_hE.allTrials(choice_trial_idx); % binary variable => no zscore
                        case 2
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'choice = high effort';
                            choice_modVals(n_choiceMods,:) = choice_hE_bis.allTrials(choice_trial_idx); % binary variable => no zscore
                        otherwise
                            error('case not ready yet');
                    end
                    
                    % reward for the non-default option
                    if choiceModel_R_varOption > 0
                        n_choiceMods = n_choiceMods + 1;
                        choice_modNames{n_choiceMods} = 'R non-default';
                        switch choiceModel_R_varOption
                            case 1 % R high effort option (amount)
                                choice_modVals(n_choiceMods,:) = raw_or_z(R_amount_varOption.allTrials(choice_trial_idx));
                            case 2 % R high effort option (level)
                                choice_modVals(n_choiceMods,:) = raw_or_z(R_level_varOption.allTrials(choice_trial_idx));
                            case 3 % z(R high effort option) (amount)
                                choice_modVals(n_choiceMods,:) = raw_or_z(z_R_amount_varOption.allTrials(choice_trial_idx));
                            case 4 % z(R high effort option) (level)
                                choice_modVals(n_choiceMods,:) = raw_or_z(z_R_level_varOption.allTrials(choice_trial_idx));
                        end
                    end
                    
                    % reward for the chosen option
                    if choiceModel_R_chosen > 0
                        n_choiceMods = n_choiceMods + 1;
                        choice_modNames{n_choiceMods} = 'R chosen';
                        switch choiceModel_R_chosen
                            case 1 % R chosen option amount
                                choice_modVals(n_choiceMods,:) = raw_or_z(R_amount_chosen.allTrials(choice_trial_idx));
                            case 2 % R chosen option level (0/1/2/3)
                                choice_modVals(n_choiceMods,:) = raw_or_z(R_level_chosen.allTrials(choice_trial_idx));
                            case 3 % z(R chosen option) amount
                                choice_modVals(n_choiceMods,:) = raw_or_z(z_R_amount_chosen.allTrials(choice_trial_idx));
                            case 4 % z(R chosen option) level (0/1/2/3)
                                choice_modVals(n_choiceMods,:) = raw_or_z(z_R_level_chosen.allTrials(choice_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % punishment for the non-default option
                    if choiceModel_P_varOption > 0
                        n_choiceMods = n_choiceMods + 1;
                        choice_modNames{n_choiceMods} = 'P non-default';
                        switch choiceModel_P_varOption
                            case 1 % R high effort option (amount)
                                choice_modVals(n_choiceMods,:) = raw_or_z(P_amount_varOption.allTrials(choice_trial_idx));
                            case 2 % R high effort option (level)
                                choice_modVals(n_choiceMods,:) = raw_or_z(P_level_varOption.allTrials(choice_trial_idx));
                            case 3 % z(R high effort option) (amount)
                                choice_modVals(n_choiceMods,:) = raw_or_z(z_P_amount_varOption.allTrials(choice_trial_idx));
                            case 4 % z(R high effort option) (level)
                                choice_modVals(n_choiceMods,:) = raw_or_z(z_P_level_varOption.allTrials(choice_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % punishment for the chosen option
                    if choiceModel_P_chosen > 0
                        n_choiceMods = n_choiceMods + 1;
                        choice_modNames{n_choiceMods} = 'P chosen';
                        switch choiceModel_P_chosen
                            case 1 % R high effort option (amount)
                                choice_modVals(n_choiceMods,:) = raw_or_z(P_amount_chosen.allTrials(choice_trial_idx));
                            case 2 % R high effort option (level)
                                choice_modVals(n_choiceMods,:) = raw_or_z(P_level_chosen.allTrials(choice_trial_idx));
                            case 3 % z(R high effort option) (amount)
                                choice_modVals(n_choiceMods,:) = raw_or_z(z_P_amount_chosen.allTrials(choice_trial_idx));
                            case 4 % z(R high effort option) (level)
                                choice_modVals(n_choiceMods,:) = raw_or_z(z_P_level_chosen.allTrials(choice_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % money left
                    if choiceModel_moneyLeft > 0
                        n_choiceMods = n_choiceMods + 1;
                        choice_modNames{n_choiceMods} = 'money left';
                        switch choiceModel_moneyLeft
                            case 1
                                choice_modVals(n_choiceMods,:) = raw_or_z(money_amount_left.allTrials(choice_trial_idx));
                            case 2
                                choice_modVals(n_choiceMods,:) = raw_or_z(abs(money_amount_left.allTrials(choice_trial_idx)));
                            case 3
                                choice_modVals(n_choiceMods,:) = raw_or_z(money_level_left.allTrials(choice_trial_idx));
                            case 4
                                choice_modVals(n_choiceMods,:) = raw_or_z(abs(money_level_left.allTrials(choice_trial_idx)));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % money right
                    if choiceModel_moneyRight > 0
                        n_choiceMods = n_choiceMods + 1;
                        choice_modNames{n_choiceMods} = 'money right';
                        switch choiceModel_moneyRight
                            case 1
                                choice_modVals(n_choiceMods,:) = raw_or_z(money_amount_right.allTrials(choice_trial_idx));
                            case 2
                                choice_modVals(n_choiceMods,:) = raw_or_z(abs(money_amount_right.allTrials(choice_trial_idx)));
                            case 3
                                choice_modVals(n_choiceMods,:) = raw_or_z(money_level_right.allTrials(choice_trial_idx));
                            case 4
                                choice_modVals(n_choiceMods,:) = raw_or_z(abs(money_level_right.allTrials(choice_trial_idx)));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % money chosen
                    if choiceModel_moneyChosen > 0
                        n_choiceMods = n_choiceMods + 1;
                        choice_modNames{n_choiceMods} = 'money chosen';
                        switch choiceModel_moneyChosen
                            case 1
                                choice_modVals(n_choiceMods,:) = raw_or_z(money_amount_chosen.allTrials(choice_trial_idx));
                            case 2
                                choice_modVals(n_choiceMods,:) = raw_or_z(abs_money_amount_chosen.allTrials(choice_trial_idx));
                            case 3
                                choice_modVals(n_choiceMods,:) = raw_or_z(money_level_chosen.allTrials(choice_trial_idx));
                            case 4
                                choice_modVals(n_choiceMods,:) = raw_or_z(abs_money_level_chosen.allTrials(choice_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % money unchosen
                    if choiceModel_moneyUnchosen > 0
                        n_choiceMods = n_choiceMods + 1;
                        choice_modNames{n_choiceMods} = 'money unchosen';
                        switch choiceModel_moneyUnchosen
                            case 1
                                choice_modVals(n_choiceMods,:) = raw_or_z(money_amount_unchosen.allTrials(choice_trial_idx));
                            case 2
                                choice_modVals(n_choiceMods,:) = raw_or_z(abs_money_amount_unchosen.allTrials(choice_trial_idx));
                            case 3
                                choice_modVals(n_choiceMods,:) = raw_or_z(money_level_unchosen.allTrials(choice_trial_idx));
                            case 4
                                choice_modVals(n_choiceMods,:) = raw_or_z(abs_money_level_unchosen.allTrials(choice_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % money associated to the non-default option
                    if choiceModel_moneyNonDefault > 0
                        n_choiceMods = n_choiceMods + 1;
                        choice_modNames{n_choiceMods} = 'money non-default';
                        switch choiceModel_moneyNonDefault
                            case 1
                                choice_modVals(n_choiceMods,:) = raw_or_z(money_amount_varOption.allTrials(choice_trial_idx));
                            case 2
                                choice_modVals(n_choiceMods,:) = raw_or_z(abs_money_amount_varOption.allTrials(choice_trial_idx));
                            case 3
                                choice_modVals(n_choiceMods,:) = raw_or_z(money_level_varOption.allTrials(choice_trial_idx));
                            case 4
                                choice_modVals(n_choiceMods,:) = raw_or_z(abs_money_level_varOption.allTrials(choice_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % money chosen - money unchosen
                    if choiceModel_moneyChosen_min_moneyUnchosen > 0
                        n_choiceMods = n_choiceMods + 1;
                        choice_modNames{n_choiceMods} = 'money chosen - unchosen';
                        switch choiceModel_moneyChosen_min_moneyUnchosen
                            case 1
                                choice_modVals(n_choiceMods,:) = raw_or_z(moneyChosen_min_moneyUnchosen_amount.allTrials(choice_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % money chosen - money fixed option (low R low E)
                    if choiceModel_moneyChosen_min_moneyDefault > 0
                        n_choiceMods = n_choiceMods + 1;
                        choice_modNames{n_choiceMods} = 'money chosen - money fixed';
                        switch choiceModel_moneyChosen_min_moneyDefault
                            case 1
                                choice_modVals(n_choiceMods,:) = raw_or_z(moneyChosen_min_moneyFixed_amount.allTrials(choice_trial_idx));
                            case 2
                                choice_modVals(n_choiceMods,:) = raw_or_z(moneyChosen_min_moneyFixed_level.allTrials(choice_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % sum of money
                    if choiceModel_moneySum > 0
                        n_choiceMods = n_choiceMods + 1;
                        choice_modNames{n_choiceMods} = 'money sum';
                        switch choiceModel_moneySum
                            case 1
                                choice_modVals(n_choiceMods,:) = raw_or_z(money_amount_sum.allTrials(choice_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % effort left
                    if choiceModel_E_left > 0
                        n_choiceMods = n_choiceMods + 1;
                        choice_modNames{n_choiceMods} = 'E left';
                        switch choiceModel_E_left
                            case 1
                                choice_modVals(n_choiceMods,:) = raw_or_z(E_left.allTrials(choice_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % effort right
                    if choiceModel_E_right > 0
                        n_choiceMods = n_choiceMods + 1;
                        choice_modNames{n_choiceMods} = 'E right';
                        switch choiceModel_E_right
                            case 1
                                choice_modVals(n_choiceMods,:) = raw_or_z(E_right.allTrials(choice_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % effort chosen
                    if choiceModel_E_chosen > 0
                        n_choiceMods = n_choiceMods + 1;
                        choice_modNames{n_choiceMods} = 'effort chosen';
                        switch choiceModel_E_chosen
                            case 1
                                choice_modVals(n_choiceMods,:) = raw_or_z(E_chosen.allTrials(choice_trial_idx));
                            case 3
                                choice_modVals(n_choiceMods,:) = raw_or_z(E_chosen_bis.allTrials(choice_trial_idx));
                            case 4
                                choice_modVals(n_choiceMods,:) = zscore(E_chosen.allTrials(choice_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % effort unchosen
                    if choiceModel_E_unchosen > 0
                        n_choiceMods = n_choiceMods + 1;
                        choice_modNames{n_choiceMods} = 'effort unchosen';
                        switch choiceModel_E_unchosen
                            case 1
                                choice_modVals(n_choiceMods,:) = raw_or_z(E_unchosen.allTrials(choice_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % effort non-default option
                    if choiceModel_E_nonDefault > 0
                        n_choiceMods = n_choiceMods + 1;
                        choice_modNames{n_choiceMods} = 'effort non-default option';
                        switch choiceModel_E_nonDefault
                            case 1
                                choice_modVals(n_choiceMods,:) = raw_or_z(E_varOption.allTrials(choice_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % effort chosen - unchosen
                    if choiceModel_Ech_min_Eunch > 0
                        n_choiceMods = n_choiceMods + 1;
                        choice_modNames{n_choiceMods} = 'effort chosen - unchosen';
                        switch choiceModel_Ech_min_Eunch
                            case 1
                                choice_modVals(n_choiceMods,:) = raw_or_z(E_chosen_min_E_unchosen.allTrials(choice_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % effort chosen - fixed option
                    if choiceModel_Ech_min_Efixed > 0
                        n_choiceMods = n_choiceMods + 1;
                        choice_modNames{n_choiceMods} = 'effort chosen - fixed';
                        switch choiceModel_Ech_min_Efixed
                            case 1
                                choice_modVals(n_choiceMods,:) = raw_or_z(Ech_min_Efixed.allTrials(choice_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % sum of efforts
                    if choiceModel_E_sum > 0
                        n_choiceMods = n_choiceMods + 1;
                        choice_modNames{n_choiceMods} = 'effort sum';
                        switch choiceModel_E_sum
                            case 1
                                choice_modVals(n_choiceMods,:) = raw_or_z(E_sum.allTrials(choice_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % money level*effort level high E option
                    if choiceModel_money_level_x_E_varOption > 0
                        n_choiceMods = n_choiceMods + 1;
                        choice_modNames{n_choiceMods} = 'money x effort non-default option';
                        switch choiceModel_money_level_x_E_varOption
                            case 1
                                choice_modVals(n_choiceMods,:) = raw_or_z(money_level_x_E_varOption.allTrials(choice_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % money level * effort chosen
                    if choiceModel_money_level_x_E_chosen > 0
                        n_choiceMods = n_choiceMods + 1;
                        choice_modNames{n_choiceMods} = 'money x effort chosen';
                        switch choiceModel_money_level_x_E_chosen
                            case 1
                                choice_modVals(n_choiceMods,:) = raw_or_z(money_level_x_E_chosen.allTrials(choice_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % R level * effort high E option
                    if choiceModel_R_level_x_E_varOption > 0
                        n_choiceMods = n_choiceMods + 1;
                        choice_modNames{n_choiceMods} = 'R x E non-default option';
                        switch choiceModel_R_level_x_E_varOption
                            case 1
                                choice_modVals(n_choiceMods,:) = raw_or_z(R_level_x_E_varOption.allTrials(choice_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % R level * effort chosen
                    if choiceModel_R_level_x_E_chosen > 0
                        n_choiceMods = n_choiceMods + 1;
                        choice_modNames{n_choiceMods} = 'R x E chosen';
                        switch choiceModel_R_level_x_E_chosen
                            case 1
                                choice_modVals(n_choiceMods,:) = raw_or_z(R_level_x_E_chosen.allTrials(choice_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % P level * effort high E option
                    if choiceModel_P_level_x_E_varOption > 0
                        n_choiceMods = n_choiceMods + 1;
                        choice_modNames{n_choiceMods} = 'P x E non-default option';
                        switch choiceModel_P_level_x_E_varOption
                            case 1
                                choice_modVals(n_choiceMods,:) = raw_or_z(P_level_x_E_varOption.allTrials(choice_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % P level * effort chosen
                    if choiceModel_P_level_x_E_chosen > 0
                        n_choiceMods = n_choiceMods + 1;
                        choice_modNames{n_choiceMods} = 'P x E chosen';
                        switch choiceModel_P_level_x_E_chosen
                            case 1
                                choice_modVals(n_choiceMods,:) = raw_or_z(P_level_x_E_chosen.allTrials(choice_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % net value chosen
                    if choiceModel_NV_chosen > 0 && ~ismember(choiceModel_NV_chosen,[4,5,6])
                        n_choiceMods = n_choiceMods + 1;
                        switch choiceModel_NV_chosen
                            case 1
                                choice_modNames{n_choiceMods} = 'NVch-NVunch';
                                choice_modVals(n_choiceMods,:) = raw_or_z(NV_ch_min_unch.allTrials(choice_trial_idx));
                            case 2
                                choice_modNames{n_choiceMods} = 'p(chosen)';
                                choice_modVals(n_choiceMods,:) = raw_or_z(pChosen.allTrials(choice_trial_idx));
                            case 3
                                choice_modNames{n_choiceMods} = 'NVch-NVunch';
                                choice_modVals(n_choiceMods,:) = raw_or_z(NV_ch_min_unch_with_bias.allTrials(choice_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % net value non-default option
                    if choiceModel_NV_varOption > 0
                        n_choiceMods = n_choiceMods + 1;
                        switch choiceModel_NV_varOption
                            case 1
                                choice_modNames{n_choiceMods} = 'delta NV high E - low E';
                                choice_modVals(n_choiceMods,:) = raw_or_z(NV_varOption.allTrials(choice_trial_idx));
                            case 2
                                choice_modNames{n_choiceMods} = '|delta NV high E - low E|';
                                choice_modVals(n_choiceMods,:) = raw_or_z(abs(NV_varOption.allTrials(choice_trial_idx)));
                            case 3
                                choice_modNames{n_choiceMods} = 'p(choice=hE)';
                                choice_modVals(n_choiceMods,:) = raw_or_z(pChoice_hE.allTrials(choice_trial_idx));
                            case 4
                                choice_modNames{n_choiceMods} = 'delta NV high E - low E + bias';
                                choice_modVals(n_choiceMods,:) = raw_or_z(NV_varOption_plus_bias.allTrials(choice_trial_idx));
                            case 5
                                choice_modNames{n_choiceMods} = '|delta NV high E - low E + bias|';
                                choice_modVals(n_choiceMods,:) = raw_or_z(abs(NV_varOption_plus_bias.allTrials(choice_trial_idx)));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % net value non-default option bis
                    if choiceModel_NV_varOption_bis > 0
                        n_choiceMods = n_choiceMods + 1;
                        switch choiceModel_NV_varOption_bis
                            case 1
                                choice_modNames{n_choiceMods} = 'delta NV high E - low E';
                                choice_modVals(n_choiceMods,:) = raw_or_z(NV_varOption.allTrials(choice_trial_idx));
                            case 2
                                choice_modNames{n_choiceMods} = '|delta NV high E - low E|';
                                choice_modVals(n_choiceMods,:) = raw_or_z(abs(NV_varOption.allTrials(choice_trial_idx)));
                            case 3
                                choice_modNames{n_choiceMods} = 'p(choice=hE)';
                                choice_modVals(n_choiceMods,:) = raw_or_z(pChoice_hE.allTrials(choice_trial_idx));
                            case 4
                                choice_modNames{n_choiceMods} = 'delta NV high E - low E + bias';
                                choice_modVals(n_choiceMods,:) = raw_or_z(NV_varOption_plus_bias.allTrials(choice_trial_idx));
                            case 5
                                choice_modNames{n_choiceMods} = '|delta NV high E - low E + bias|';
                                choice_modVals(n_choiceMods,:) = raw_or_z(abs(NV_varOption_plus_bias.allTrials(choice_trial_idx)));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    if strcmp(task_id,'Ep')
                        % force integral
                        if choiceModel_F_integral > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'F integral';
                            switch choiceModel_F_integral
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(AUC.allTrials(choice_trial_idx));
                                case 2
                                    choice_modVals(n_choiceMods,:) = raw_or_z(AUC_overshoot.allTrials(choice_trial_idx));
                                case 3
                                    choice_modVals(n_choiceMods,:) = raw_or_z(AUC_N.allTrials(choice_trial_idx));
                                case 4
                                    choice_modVals(n_choiceMods,:) = raw_or_z(AUC_overshoot_N.allTrials(choice_trial_idx));
                                case 5
                                    choice_modVals(n_choiceMods,:) = zscore(AUC.allTrials(choice_trial_idx));
                                case 6
                                    choice_modVals(n_choiceMods,:) = zscore(AUC_overshoot.allTrials(choice_trial_idx));
                                case 7
                                    choice_modVals(n_choiceMods,:) = zscore(AUC_N.allTrials(choice_trial_idx));
                                case 8
                                    choice_modVals(n_choiceMods,:) = zscore(AUC_overshoot_N.allTrials(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % fatigue
                        if choiceModel_fatigue > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'fatigue';
                            switch choiceModel_fatigue
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(fatigue.allTrials(choice_trial_idx));
                                case 2
                                    choice_modVals(n_choiceMods,:) = zscore(fatigue.allTrials(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % Ech*fatigue
                        if choiceModel_Ech_x_fatigue > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'Ech_x_fatigue';
                            switch choiceModel_Ech_x_fatigue
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(E_chosen.allTrials(choice_trial_idx).*fatigue.allTrials(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                    end % physical effort filter
                    
                    if strcmp(task_id,'Em')
                        % efficacy
                        if choiceModel_efficacy > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'efficacy';
                            switch choiceModel_efficacy
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(efficacy_with2first.allTrials(choice_trial_idx));
                                case 2
                                    choice_modVals(n_choiceMods,:) = raw_or_z(efficacy_pureNback.allTrials(choice_trial_idx));
                                case 3
                                    choice_modVals(n_choiceMods,:) = raw_or_z(efficacy_bis_with2first.allTrials(choice_trial_idx));
                                case 4
                                    choice_modVals(n_choiceMods,:) = raw_or_z(efficacy_bis_pureNback.allTrials(choice_trial_idx));
                                case 5
                                    choice_modVals(n_choiceMods,:) = zscore(efficacy_with2first.allTrials(choice_trial_idx));
                                case 6
                                    choice_modVals(n_choiceMods,:) = zscore(efficacy_pureNback.allTrials(choice_trial_idx));
                                case 7
                                    choice_modVals(n_choiceMods,:) = zscore(efficacy_bis_with2first.allTrials(choice_trial_idx));
                                case 8
                                    choice_modVals(n_choiceMods,:) = zscore(efficacy_bis_pureNback.allTrials(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % previous trial efficacy
                        if choiceModel_prevEfficacy > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'previous trial efficacy';
                            switch choiceModel_prevEfficacy
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(prevEfficacy_with2first.allTrials(choice_trial_idx));
                                case 2
                                    choice_modVals(n_choiceMods,:) = raw_or_z(prevEfficacy_pureNback.allTrials(choice_trial_idx));
                                case 3
                                    choice_modVals(n_choiceMods,:) = raw_or_z(prevEfficacy_bis_with2first.allTrials(choice_trial_idx));
                                case 4
                                    choice_modVals(n_choiceMods,:) = raw_or_z(prevEfficacy_bis_pureNback.allTrials(choice_trial_idx));
                                case 5
                                    choice_modVals(n_choiceMods,:) = zscore(prevEfficacy_with2first.allTrials(choice_trial_idx));
                                case 6
                                    choice_modVals(n_choiceMods,:) = zscore(prevEfficacy_pureNback.allTrials(choice_trial_idx));
                                case 7
                                    choice_modVals(n_choiceMods,:) = zscore(prevEfficacy_bis_with2first.allTrials(choice_trial_idx));
                                case 8
                                    choice_modVals(n_choiceMods,:) = zscore(prevEfficacy_bis_pureNback.allTrials(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % Ech*(previous trial efficacy)
                        if choiceModel_Ech_x_prevEfficacy > 0
                            n_choiceMods = n_choiceMods + 1;
                            choice_modNames{n_choiceMods} = 'Ech_x_previous trial efficacy';
                            switch choiceModel_Ech_x_prevEfficacy
                                case 1
                                    choice_modVals(n_choiceMods,:) = raw_or_z(E_chosen.allTrials(choice_trial_idx).*prevEfficacy_with2first.allTrials(choice_trial_idx));
                                case 2
                                    choice_modVals(n_choiceMods,:) = raw_or_z(E_chosen.allTrials(choice_trial_idx).*prevEfficacy_pureNback.allTrials(choice_trial_idx));
                                case 3
                                    choice_modVals(n_choiceMods,:) = raw_or_z(E_chosen.allTrials(choice_trial_idx).*prevEfficacy_bis_with2first.allTrials(choice_trial_idx));
                                case 4
                                    choice_modVals(n_choiceMods,:) = raw_or_z(E_chosen.allTrials(choice_trial_idx).*prevEfficacy_bis_pureNback.allTrials(choice_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                    end % mental effort filter
                    
                    % trial number
                    if choiceModel_trialN > 0
                        n_choiceMods = n_choiceMods + 1;
                        choice_modNames{n_choiceMods} = 'trial number';
                        switch choiceModel_trialN
                            case 1
                                choice_modVals(n_choiceMods,:) = raw_or_z(trialN.allTrials(choice_trial_idx));
                            case 2
                                choice_modVals(n_choiceMods,:) = raw_or_z(trialN_dEch.allTrials(choice_trial_idx));
                            case 3
                                choice_modVals(n_choiceMods,:) = raw_or_z(trialN_dEnonDef_min_Edef.allTrials(choice_trial_idx));
                            case 4
                                choice_modVals(n_choiceMods,:) = raw_or_z(trialN_dEnonDef.allTrials(choice_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % choice confidence
                    if choiceModel_conf > 0
                        n_choiceMods = n_choiceMods + 1;
                        choice_modNames{n_choiceMods} = 'confidence';
                        switch choiceModel_conf
                            case 1 % binary variable => no zscore
                                choice_modVals(n_choiceMods,:) = confidence.allTrials(choice_trial_idx);
                            case {2,3,4} % confidence inferred by the model => ok to zscore
                                choice_modVals(n_choiceMods,:) = raw_or_z(confidence.allTrials(choice_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % RT (last regressor)
                    if choiceModel_RT > 0 && ~ismember(choiceModel_RT,[4,5,6])
                        n_choiceMods = n_choiceMods + 1;
                        choice_modNames{n_choiceMods} = 'choice RT';
                        switch choiceModel_RT
                            case 1
                                choice_modVals(n_choiceMods,:) = raw_or_z(choice_RT.allTrials(choice_trial_idx));
                            case 2
                                choice_modVals(n_choiceMods,:) = zscore(choice_RT.allTrials(choice_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    [matlabbatch] = First_level_loadEachCondition_DCM(matlabbatch, sub_idx, iCond,...
                        ['choice_',RP_dispChoice_nm,'_',splitE_dispChoice_nm], modelChoiceOnset, modelChoiceDur,...
                        n_choiceMods, choice_modNames, choice_modVals,...
                        orth_vars, onsets_only_GLM);
            end % DCM mode: info about how regressors concatenated or not across tasks or sessions
        end % E condition
    end % RP
end % model choice

%% chosen period
chosenModel = GLMprm.model_onset.(task_id).chosen;
if ismember(chosenModel,{'stick','boxcar','boxcar_bis','boxcar_ter'})
    
    for iRP_chosen = 1:length(RPchosenCond)
        RP_dispChosen_nm = RPchosenCond{iRP_chosen};
        for iEsplit_dispChosen = 1:length(EsplitchosenCond)
            splitE_dispChosen_nm = EsplitchosenCond{iEsplit_dispChosen};
            chosenModel_R_vs_P          = GLMprm.chosen.(task_id).(RP_dispChosen_nm).(splitE_dispChosen_nm).R_vs_P;
            chosenModel_choicehE        = GLMprm.chosen.(task_id).(RP_dispChosen_nm).(splitE_dispChosen_nm).choiceHighE;
            chosenModel_R_varOption     = GLMprm.chosen.(task_id).(RP_dispChosen_nm).(splitE_dispChosen_nm).R_varOption;
            chosenModel_P_varOption     = GLMprm.chosen.(task_id).(RP_dispChosen_nm).(splitE_dispChosen_nm).P_varOption;
            chosenModel_R_chosen        = GLMprm.chosen.(task_id).(RP_dispChosen_nm).(splitE_dispChosen_nm).R_chosen;
            chosenModel_P_chosen        = GLMprm.chosen.(task_id).(RP_dispChosen_nm).(splitE_dispChosen_nm).P_chosen;
            chosenModel_moneyLeft       = GLMprm.choice.(task_id).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_left;
            chosenModel_moneyRight      = GLMprm.choice.(task_id).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_right;
            chosenModel_moneyChosen     = GLMprm.chosen.(task_id).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_chosen;
            chosenModel_moneyUnchosen   = GLMprm.chosen.(task_id).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_unchosen;
            chosenModel_moneyNonDefault = GLMprm.chosen.(task_id).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_varOption;
            chosenModel_money_chosen_min_money_unchosen = GLMprm.chosen.(task_id).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_ch_min_unch;
            chosenModel_moneyChosen_min_moneyDefault    = GLMprm.chosen.(task_id).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_ch_min_fixOption;
            chosenModel_moneySum        = GLMprm.chosen.(task_id).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_sum;
            chosenModel_E_left          = GLMprm.chosen.(task_id).(RP_dispChosen_nm).(splitE_dispChosen_nm).E_left;
            chosenModel_E_right         = GLMprm.chosen.(task_id).(RP_dispChosen_nm).(splitE_dispChosen_nm).E_right;
            chosenModel_Echosen         = GLMprm.chosen.(task_id).(RP_dispChosen_nm).(splitE_dispChosen_nm).E_chosen;
            chosenModel_Eunchosen       = GLMprm.chosen.(task_id).(RP_dispChosen_nm).(splitE_dispChosen_nm).E_unchosen;
            chosenModel_EnonDefault     = GLMprm.chosen.(task_id).(RP_dispChosen_nm).(splitE_dispChosen_nm).E_varOption;
            chosenModel_Ech_min_Eunch   = GLMprm.chosen.(task_id).(RP_dispChosen_nm).(splitE_dispChosen_nm).E_ch_min_unch;
            chosenModel_Ech_min_Efixed  = GLMprm.chosen.(task_id).(RP_dispChosen_nm).(splitE_dispChosen_nm).E_ch_min_fixOption;
            chosenModel_E_sum           = GLMprm.chosen.(task_id).(RP_dispChosen_nm).(splitE_dispChosen_nm).E_sum;
            chosenModel_money_level_x_E_varOption = GLMprm.chosen.(task_id).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_level_x_E_varOption;
            chosenModel_money_level_x_E_chosen    = GLMprm.chosen.(task_id).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_level_x_E_chosen;
            chosenModel_R_level_x_E_varOption     = GLMprm.chosen.(task_id).(RP_dispChosen_nm).(splitE_dispChosen_nm).R_level_x_E_varOption;
            chosenModel_R_level_x_E_chosen        = GLMprm.chosen.(task_id).(RP_dispChosen_nm).(splitE_dispChosen_nm).R_level_x_E_chosen;
            chosenModel_P_level_x_E_varOption     = GLMprm.chosen.(task_id).(RP_dispChosen_nm).(splitE_dispChosen_nm).P_level_x_E_varOption;
            chosenModel_P_level_x_E_chosen        = GLMprm.chosen.(task_id).(RP_dispChosen_nm).(splitE_dispChosen_nm).P_level_x_E_chosen;
            chosenModel_NV_chosen           = GLMprm.chosen.(task_id).(RP_dispChosen_nm).(splitE_dispChosen_nm).NV_chosen;
            chosenModel_NV_varOption        = GLMprm.chosen.(task_id).(RP_dispChosen_nm).(splitE_dispChosen_nm).NV_varOption;
            chosenModel_NV_varOption_bis    = GLMprm.chosen.(task_id).(RP_dispChosen_nm).(splitE_dispChosen_nm).NV_varOption_bis;
            switch task_id
                case 'Ep'
                    chosenModel_F_integral = GLMprm.chosen.(task_id).(RP_dispChosen_nm).(splitE_dispChosen_nm).F_integral;
                    chosenModel_fatigue = GLMprm.chosen.(task_id).(RP_dispChosen_nm).(splitE_dispChosen_nm).fatigue;
                    chosenModel_Ech_x_fatigue = GLMprm.chosen.(task_id).(RP_dispChosen_nm).(splitE_dispChosen_nm).Ech_x_fatigue;
                case 'Em'
                    chosenModel_efficacy = GLMprm.chosen.(task_id).(RP_dispChosen_nm).(splitE_dispChosen_nm).efficacy;
                    chosenModel_prevEfficacy = GLMprm.chosen.(task_id).(RP_dispChosen_nm).(splitE_dispChosen_nm).prevEfficacy;
                    chosenModel_Ech_x_prevEfficacy = GLMprm.chosen.(task_id).(RP_dispChosen_nm).(splitE_dispChosen_nm).Ech_x_prevEfficacy;
            end
            chosenModel_confidence      = GLMprm.chosen.(task_id).(RP_dispChosen_nm).(splitE_dispChosen_nm).confidence;
            chosenModel_RT              = GLMprm.chosen.(task_id).(RP_dispChosen_nm).(splitE_dispChosen_nm).RT;
            chosenModel_trialN          = GLMprm.chosen.(task_id).(RP_dispChosen_nm).(splitE_dispChosen_nm).trialN;
            
            %% adapt depending on DCM_mode
            switch DCM_mode
                case 1 % all sessions independent
                    for iRun = 1:n_runs
                        run_full_nm = ['run',num2str(iRun)];
                        
                        % extract trial index for the current loop
                        switch RP_dispChosen_nm
                            case 'RP'
                                RPfilter_dispChosen = true(1,length(double(RP_var_binary.(run_full_nm))));
                            case 'R'
                                RPfilter_dispChosen = (RP_var_binary.(run_full_nm) == 1);
                            case 'P'
                                RPfilter_dispChosen = (RP_var_binary.(run_full_nm) == 0);
                        end
                        switch splitE_dispChosen_nm
                            case 'E'
                                Efilter_dispChosen = true(1,length(double(RP_var_binary.(run_full_nm))));
                            case 'E1'
                                Efilter_dispChosen = (E_varOption.(run_full_nm) == 1);
                            case 'E2'
                                Efilter_dispChosen = (E_varOption.(run_full_nm) == 2);
                            case 'E3'
                                Efilter_dispChosen = (E_varOption.(run_full_nm) == 3);
                            case 'Ech0'
                                Efilter_dispChosen = (E_chosen.(run_full_nm) == 0);
                            case 'Ech1'
                                Efilter_dispChosen = (E_chosen.(run_full_nm) == 1);
                            case 'Ech2'
                                Efilter_dispChosen = (E_chosen.(run_full_nm) == 2);
                            case 'Ech3'
                                Efilter_dispChosen = (E_chosen.(run_full_nm) == 3);
                            case 'lEch'
                                Efilter_dispChosen = (choice_hE.(run_full_nm) == 0);
                            case 'hEch'
                                Efilter_dispChosen = (choice_hE.(run_full_nm) == 1);
                        end
                        chosen_trial_idx = (RPfilter_dispChosen.*Efilter_dispChosen) == 1; % NEED to transform it into logical or will just focus on the first trial
                        
                        %% chosen onset
                        iCond = iCond + 1;
                        modelChosenOnset = onsets.dispChosenOnsets.(run_full_nm)(chosen_trial_idx);
                        % duration
                        switch chosenModel
                            case 'stick'
                                modelChosenDur = 0;
                            case 'boxcar' % duration displaying the chosen option
                                modelChosenDur = durations.dispChosenDur.(run_full_nm)(chosen_trial_idx);
                            case 'boxcar_bis' % duration going form display of chosen option
                                % until the end of the exertion of the effort
                                modelChosenDur = durations.dispChosenDur.(run_full_nm)(chosen_trial_idx) +...
                                    durations.preEffortCrossDur.(run_full_nm)(chosen_trial_idx) +...
                                    durations.EperfDur.(run_full_nm)(chosen_trial_idx);
                            case 'boxcar_ter' % duration going form display of chosen option
                                % until the end of the effort preparation cross
                                % (entailing the whole effort preparation period)
                                modelChosenDur = durations.dispChosenDur.(run_full_nm)(chosen_trial_idx) +...
                                    durations.preEffortCrossDur.(run_full_nm)(chosen_trial_idx);
                        end
                        
                        %% chosen modulators
                        n_chosenMods = 0;
                        chosen_modNames = cell(1,1);
                        chosen_modVals = [];
                        
                        % RT (first regressor)
                        if chosenModel_RT > 0 && ~ismember(chosenModel_RT,[1,2,3])
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'choice RT';
                            switch chosenModel_RT
                                case 4
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(choice_RT.(run_full_nm)(chosen_trial_idx));
                                case 5
                                    chosen_modVals(n_chosenMods,:) = zscore(choice_RT.(run_full_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % net value chosen (first regressor)
                        if chosenModel_NV_chosen > 0 && ~ismember(chosenModel_NV_chosen,[1,2,3])
                            n_chosenMods = n_chosenMods + 1;
                            switch chosenModel_NV_chosen
                                case 4
                                    chosen_modNames{n_chosenMods} = 'NVch-NVunch';
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(NV_ch_min_unch.(run_full_nm)(chosen_trial_idx));
                                case 5
                                    chosen_modNames{n_chosenMods} = 'p(chosen)';
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(pChosen.(run_full_nm)(chosen_trial_idx));
                                case 6
                                    chosen_modNames{n_chosenMods} = 'NVch-NVunch';
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(NV_ch_min_unch_with_bias.(run_full_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % reward vs punishment trials
                        if chosenModel_R_vs_P > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'R_vs_P';
                            switch chosenModel_R_vs_P
                                case 1
                                    chosen_modVals(n_chosenMods,:) = RP_var_binary.(run_full_nm)(chosen_trial_idx);
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % choice = high effort
                        switch chosenModel_choicehE
                            case 0
                            case 1
                                n_chosenMods = n_chosenMods + 1;
                                chosen_modNames{n_chosenMods} = 'choice = high effort';
                                chosen_modVals(n_chosenMods,:) = choice_hE.(run_full_nm)(chosen_trial_idx); % binary variable => no zscore
                            case 2
                                n_chosenMods = n_chosenMods + 1;
                                chosen_modNames{n_chosenMods} = 'choice = high effort';
                                chosen_modVals(n_chosenMods,:) = choice_hE_bis.(run_full_nm)(chosen_trial_idx); % binary variable => no zscore
                            otherwise
                                error('not ready yet');
                        end
                        
                        % reward for the non-default option
                        if chosenModel_R_varOption > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'R non-default';
                            switch chosenModel_R_varOption
                                case 1 % R high effort option (amount)
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(R_amount_varOption.(run_full_nm)(chosen_trial_idx));
                                case 2 % R high effort option (level)
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(R_level_varOption.(run_full_nm)(chosen_trial_idx));
                                case 3 % z(R high effort option) (amount)
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(z_R_amount_varOption.(run_full_nm)(chosen_trial_idx));
                                case 4 % z(R high effort option) (level)
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(z_R_level_varOption.(run_full_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % reward for the chosen option
                        if chosenModel_R_chosen > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'R chosen';
                            switch chosenModel_R_chosen
                                case 1 % R high effort option (amount)
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(R_amount_chosen.(run_full_nm)(chosen_trial_idx));
                                case 2 % R high effort option (level)
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(R_level_chosen.(run_full_nm)(chosen_trial_idx));
                                case 3 % z(R high effort option) (amount)
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(z_R_amount_chosen.(run_full_nm)(chosen_trial_idx));
                                case 4 % z(R high effort option) (level)
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(z_R_level_chosen.(run_full_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % punishment for the non-default option
                        if chosenModel_P_varOption > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'P non-default';
                            switch chosenModel_P_varOption
                                case 1 % R high effort option (amount)
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(P_amount_varOption.(run_full_nm)(chosen_trial_idx));
                                case 2 % R high effort option (level)
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(P_level_varOption.(run_full_nm)(chosen_trial_idx));
                                case 3 % z(R high effort option) (amount)
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(z_P_amount_varOption.(run_full_nm)(chosen_trial_idx));
                                case 4 % z(R high effort option) (level)
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(z_P_level_varOption.(run_full_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % punishment for the chosen option
                        if chosenModel_P_chosen > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'P chosen';
                            switch chosenModel_P_chosen
                                case 1 % R high effort option (amount)
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(P_amount_chosen.(run_full_nm)(chosen_trial_idx));
                                case 2 % R high effort option (level)
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(P_level_chosen.(run_full_nm)(chosen_trial_idx));
                                case 3 % z(R high effort option) (amount)
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(z_P_amount_chosen.(run_full_nm)(chosen_trial_idx));
                                case 4 % z(R high effort option) (level)
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(z_P_level_chosen.(run_full_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % money left
                        if chosenModel_moneyLeft > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'money left';
                            switch chosenModel_moneyLeft
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(money_amount_left.(run_full_nm)(chosen_trial_idx));
                                case 2
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(abs(money_amount_left.(run_full_nm)(chosen_trial_idx)));
                                case 3
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(money_level_left.(run_full_nm)(chosen_trial_idx));
                                case 4
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(abs(money_level_chosen.(run_full_nm)(chosen_trial_idx)));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % money right
                        if chosenModel_moneyRight > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'money right';
                            switch chosenModel_moneyRight
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(money_amount_right.(run_full_nm)(chosen_trial_idx));
                                case 2
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(abs(money_amount_right.(run_full_nm)(chosen_trial_idx)));
                                case 3
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(money_level_right.(run_full_nm)(chosen_trial_idx));
                                case 4
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(abs(money_level_right.(run_full_nm)(chosen_trial_idx)));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % money chosen
                        if chosenModel_moneyChosen > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'money chosen';
                            switch chosenModel_moneyChosen
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(money_amount_chosen.(run_full_nm)(chosen_trial_idx));
                                case 2
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(abs_money_amount_chosen.(run_full_nm)(chosen_trial_idx));
                                case 3
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(money_level_chosen.(run_full_nm)(chosen_trial_idx));
                                case 4
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(abs_money_level_chosen.(run_full_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % money unchosen
                        if chosenModel_moneyUnchosen > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'money unchosen';
                            switch chosenModel_moneyUnchosen
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(money_amount_unchosen.(run_full_nm)(chosen_trial_idx));
                                case 2
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(abs_money_amount_unchosen.(run_full_nm)(chosen_trial_idx));
                                case 3
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(money_level_unchosen.(run_full_nm)(chosen_trial_idx));
                                case 4
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(abs_money_level_unchosen.(run_full_nm)(chosen_trial_idx));
                                otherwise
                                    error('ready yet');
                            end
                        end
                        
                        % money non-default option
                        if chosenModel_moneyNonDefault > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'money non-default option';
                            switch chosenModel_moneyNonDefault
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(money_amount_varOption.(run_full_nm)(chosen_trial_idx));
                                case 2
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(abs_money_amount_varOption.(run_full_nm)(chosen_trial_idx));
                                case 3
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(money_level_varOption.(run_full_nm)(chosen_trial_idx));
                                case 4
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(abs_money_level_varOption.(run_full_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % money chosen - money unchosen
                        if chosenModel_money_chosen_min_money_unchosen > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'money chosen - unchosen';
                            switch chosenModel_money_chosen_min_money_unchosen
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(moneyChosen_min_moneyUnchosen_amount.(run_full_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % money chosen - money fixed option (low R low E)
                        if chosenModel_moneyChosen_min_moneyDefault > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'money chosen - money fixed';
                            switch chosenModel_moneyChosen_min_moneyDefault
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(moneyChosen_min_moneyFixed_amount.(run_full_nm)(chosen_trial_idx));
                                case 2
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(moneyChosen_min_moneyFixed_level.(run_full_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % sum of money
                        if chosenModel_moneySum > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'money sum';
                            switch chosenModel_moneySum
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(money_amount_sum.(run_full_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % effort left
                        if chosenModel_E_left > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'effort left';
                            switch chosenModel_E_left
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(E_left.(run_full_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % effort right
                        if chosenModel_E_right > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'effort right';
                            switch chosenModel_E_right
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(E_right.(run_full_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % effort chosen
                        if chosenModel_Echosen > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'effort chosen';
                            switch chosenModel_Echosen
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(E_chosen.(run_full_nm)(chosen_trial_idx));
                                case 3
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(E_chosen_bis.(run_full_nm)(chosen_trial_idx));
                                case 4
                                    chosen_modVals(n_chosenMods,:) = zscore(E_chosen.(run_full_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % effort unchosen
                        if chosenModel_Eunchosen > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'effort unchosen';
                            switch chosenModel_Eunchosen
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(E_unchosen.(run_full_nm)(chosen_trial_idx));
                                otherwise
                                    error('ready yet');
                            end
                        end
                        
                        % effort non-default option
                        if chosenModel_EnonDefault > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'effort non-default option';
                            switch chosenModel_EnonDefault
                                case 1
                                    chosen_modVals(n_chosenMods,:) = E_varOption.(run_full_nm)(chosen_trial_idx); % binary variable => no zscore
                                otherwise
                                    error('ready yet');
                            end
                        end
                        
                        % effort chosen - unchosen
                        if chosenModel_Ech_min_Eunch > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'effort chosen - unchosen';
                            switch chosenModel_Ech_min_Eunch
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(E_chosen_min_E_unchosen.(run_full_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % effort chosen - fixed option
                        if chosenModel_Ech_min_Efixed > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'effort chosen - fixed';
                            switch chosenModel_Ech_min_Efixed
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(Ech_min_Efixed.(run_full_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % sum of efforts
                        if chosenModel_E_sum > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'effort sum';
                            switch chosenModel_E_sum
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(E_sum.(run_full_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % money level*effort level high E option
                        if chosenModel_money_level_x_E_varOption > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'money x effort non-default option';
                            switch chosenModel_money_level_x_E_varOption
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(money_level_x_E_varOption.(run_full_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % money level * effort chosen
                        if chosenModel_money_level_x_E_chosen > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'money x effort chosen';
                            switch chosenModel_money_level_x_E_chosen
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(money_level_x_E_chosen.(run_full_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % R level * effort high E option
                        if chosenModel_R_level_x_E_varOption > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'R x E non-default option';
                            switch chosenModel_R_level_x_E_varOption
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(R_level_x_E_varOption.(run_full_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % R level * effort chosen
                        if chosenModel_R_level_x_E_chosen > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'R x E chosen';
                            switch chosenModel_R_level_x_E_chosen
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(R_level_x_E_chosen.(run_full_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % P level * effort high E option
                        if chosenModel_P_level_x_E_varOption > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'P x E non-default option';
                            switch chosenModel_P_level_x_E_varOption
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(P_level_x_E_varOption.(run_full_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % P level * effort chosen
                        if chosenModel_P_level_x_E_chosen > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'P x E chosen';
                            switch chosenModel_P_level_x_E_chosen
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(P_level_x_E_chosen.(run_full_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % net value chosen
                        if chosenModel_NV_chosen > 0 && ~ismember(chosenModel_NV_chosen,[4,5,6])
                            n_chosenMods = n_chosenMods + 1;
                            switch chosenModel_NV_chosen
                                case 1
                                    chosen_modNames{n_chosenMods} = 'NVch-NVunch';
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(NV_ch_min_unch.(run_full_nm)(chosen_trial_idx));
                                case 2
                                    chosen_modNames{n_chosenMods} = 'p(chosen)';
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(pChosen.(run_full_nm)(chosen_trial_idx));
                                case 3
                                    chosen_modNames{n_chosenMods} = 'NVch-NVunch';
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(NV_ch_min_unch_with_bias.(run_full_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % net value non-default option
                        if chosenModel_NV_varOption > 0
                            n_chosenMods = n_chosenMods + 1;
                            switch chosenModel_NV_varOption
                                case 1
                                    chosen_modNames{n_chosenMods} = 'delta NV high E - low E';
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(NV_varOption.(run_full_nm)(chosen_trial_idx));
                                case 2
                                    chosen_modNames{n_chosenMods} = '|delta NV high E - low E|';
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(abs(NV_varOption.(run_full_nm)(chosen_trial_idx)));
                                case 3
                                    chosen_modNames{n_chosenMods} = 'p(choice=hE)';
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(abs(pChoice_hE.(run_full_nm)(chosen_trial_idx)));
                                case 4
                                    chosen_modNames{n_chosenMods} = 'delta NV high E - low E + bias';
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(NV_varOption_plus_bias.(run_full_nm)(chosen_trial_idx));
                                case 5
                                    chosen_modNames{n_chosenMods} = '|delta NV high E - low E + bias|';
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(abs(NV_varOption_plus_bias.(run_full_nm)(chosen_trial_idx)));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % net value non-default option bis
                        if chosenModel_NV_varOption_bis > 0
                            n_chosenMods = n_chosenMods + 1;
                            switch chosenModel_NV_varOption_bis
                                case 1
                                    chosen_modNames{n_chosenMods} = 'delta NV high E - low E';
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(NV_varOption.(run_full_nm)(chosen_trial_idx));
                                case 2
                                    chosen_modNames{n_chosenMods} = '|delta NV high E - low E|';
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(abs(NV_varOption.(run_full_nm)(chosen_trial_idx)));
                                case 3
                                    chosen_modNames{n_chosenMods} = 'p(choice=hE)';
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(abs(pChoice_hE.(run_full_nm)(chosen_trial_idx)));
                                case 4
                                    chosen_modNames{n_chosenMods} = 'delta NV high E - low E + bias';
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(NV_varOption_plus_bias.(run_full_nm)(chosen_trial_idx));
                                case 5
                                    chosen_modNames{n_chosenMods} = '|delta NV high E - low E + bias|';
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(abs(NV_varOption_plus_bias.(run_full_nm)(chosen_trial_idx)));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        if strcmp(task_id,'Ep')
                            % force integral
                            if chosenModel_F_integral > 0
                                n_chosenMods = n_chosenMods + 1;
                                chosen_modNames{n_chosenMods} = 'F integral';
                                switch chosenModel_F_integral
                                    case 1
                                        chosen_modVals(n_chosenMods,:) = raw_or_z(AUC.(run_full_nm)(chosen_trial_idx));
                                    case 2
                                        chosen_modVals(n_chosenMods,:) = raw_or_z(AUC_overshoot.(run_full_nm)(chosen_trial_idx));
                                    case 3
                                        chosen_modVals(n_chosenMods,:) = raw_or_z(AUC_N.(run_full_nm)(chosen_trial_idx));
                                    case 4
                                        chosen_modVals(n_chosenMods,:) = raw_or_z(AUC_overshoot_N.(run_full_nm)(chosen_trial_idx));
                                    case 5
                                        chosen_modVals(n_chosenMods,:) = zscore(AUC.(run_full_nm)(chosen_trial_idx));
                                    case 6
                                        chosen_modVals(n_chosenMods,:) = zscore(AUC_overshoot.(run_full_nm)(chosen_trial_idx));
                                    case 7
                                        chosen_modVals(n_chosenMods,:) = zscore(AUC_N.(run_full_nm)(chosen_trial_idx));
                                    case 8
                                        chosen_modVals(n_chosenMods,:) = zscore(AUC_overshoot_N.(run_full_nm)(chosen_trial_idx));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                            
                            % fatigue
                            if chosenModel_fatigue > 0
                                n_chosenMods = n_chosenMods + 1;
                                chosen_modNames{n_chosenMods} = 'fatigue';
                                switch chosenModel_fatigue
                                    case 1
                                        chosen_modVals(n_chosenMods,:) = raw_or_z(fatigue.(run_full_nm)(chosen_trial_idx));
                                    case 2
                                        chosen_modVals(n_chosenMods,:) = zscore(fatigue.(run_full_nm)(chosen_trial_idx));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                            
                            % Ech*fatigue
                            if chosenModel_Ech_x_fatigue > 0
                                n_chosenMods = n_chosenMods + 1;
                                chosen_modNames{n_chosenMods} = 'Ech_x_fatigue';
                                switch chosenModel_Ech_x_fatigue
                                    case 1
                                        chosen_modVals(n_chosenMods,:) = raw_or_z(E_chosen.(run_full_nm)(chosen_trial_idx).*fatigue.(run_full_nm)(chosen_trial_idx));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                        end % physical effort filter
                        
                        if strcmp(task_id,'Em')
                            % efficacy
                            if chosenModel_efficacy > 0
                                n_chosenMods = n_chosenMods + 1;
                                chosen_modNames{n_chosenMods} = 'efficacy';
                                switch chosenModel_efficacy
                                    case 1
                                        chosen_modVals(n_chosenMods,:) = raw_or_z(efficacy_with2first.(run_full_nm)(chosen_trial_idx));
                                    case 2
                                        chosen_modVals(n_chosenMods,:) = raw_or_z(efficacy_pureNback.(run_full_nm)(chosen_trial_idx));
                                    case 3
                                        chosen_modVals(n_chosenMods,:) = raw_or_z(efficacy_bis_with2first.(run_full_nm)(chosen_trial_idx));
                                    case 4
                                        chosen_modVals(n_chosenMods,:) = raw_or_z(efficacy_bis_pureNback.(run_full_nm)(chosen_trial_idx));
                                    case 5
                                        chosen_modVals(n_chosenMods,:) = zscore(efficacy_with2first.(run_full_nm)(chosen_trial_idx));
                                    case 6
                                        chosen_modVals(n_chosenMods,:) = zscore(efficacy_pureNback.(run_full_nm)(chosen_trial_idx));
                                    case 7
                                        chosen_modVals(n_chosenMods,:) = zscore(efficacy_bis_with2first.(run_full_nm)(chosen_trial_idx));
                                    case 8
                                        chosen_modVals(n_chosenMods,:) = zscore(efficacy_bis_pureNback.(run_full_nm)(chosen_trial_idx));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                            
                            % previous trial efficacy
                            if chosenModel_prevEfficacy > 0
                                n_chosenMods = n_chosenMods + 1;
                                chosen_modNames{n_chosenMods} = 'previous trial efficacy';
                                switch chosenModel_prevEfficacy
                                    case 1
                                        chosen_modVals(n_chosenMods,:) = raw_or_z(prevEfficacy_with2first.(run_full_nm)(chosen_trial_idx));
                                    case 2
                                        chosen_modVals(n_chosenMods,:) = raw_or_z(prevEfficacy_pureNback.(run_full_nm)(chosen_trial_idx));
                                    case 3
                                        chosen_modVals(n_chosenMods,:) = raw_or_z(prevEfficacy_bis_with2first.(run_full_nm)(chosen_trial_idx));
                                    case 4
                                        chosen_modVals(n_chosenMods,:) = raw_or_z(prevEfficacy_bis_pureNback.(run_full_nm)(chosen_trial_idx));
                                    case 5
                                        chosen_modVals(n_chosenMods,:) = zscore(prevEfficacy_with2first.(run_full_nm)(chosen_trial_idx));
                                    case 6
                                        chosen_modVals(n_chosenMods,:) = zscore(prevEfficacy_pureNback.(run_full_nm)(chosen_trial_idx));
                                    case 7
                                        chosen_modVals(n_chosenMods,:) = zscore(prevEfficacy_bis_with2first.(run_full_nm)(chosen_trial_idx));
                                    case 8
                                        chosen_modVals(n_chosenMods,:) = zscore(prevEfficacy_bis_pureNback.(run_full_nm)(chosen_trial_idx));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                            
                            % Ech*(previous trial efficacy)
                            if chosenModel_Ech_x_prevEfficacy > 0
                                n_chosenMods = n_chosenMods + 1;
                                chosen_modNames{n_chosenMods} = 'Ech_x_previous trial efficacy';
                                switch chosenModel_Ech_x_prevEfficacy
                                    case 1
                                        chosen_modVals(n_chosenMods,:) = raw_or_z(Echosen.(run_full_nm)(chosen_trial_idx).*prevEfficacy_with2first.(run_full_nm)(chosen_trial_idx));
                                    case 2
                                        chosen_modVals(n_chosenMods,:) = raw_or_z(Echosen.(run_full_nm)(chosen_trial_idx).*prevEfficacy_pureNback.(run_full_nm)(chosen_trial_idx));
                                    case 3
                                        chosen_modVals(n_chosenMods,:) = raw_or_z(Echosen.(run_full_nm)(chosen_trial_idx).*prevEfficacy_bis_with2first.(run_full_nm)(chosen_trial_idx));
                                    case 4
                                        chosen_modVals(n_chosenMods,:) = raw_or_z(Echosen.(run_full_nm)(chosen_trial_idx).*prevEfficacy_bis_pureNback.(run_full_nm)(chosen_trial_idx));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                        end % mental effort filter
                        
                        % trial number
                        if chosenModel_trialN > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'trial number';
                            switch chosenModel_trialN
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(trialN.(run_full_nm)(chosen_trial_idx));
                                case 2
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(trialN_dEch.(run_full_nm)(chosen_trial_idx));
                                case 3
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(trialN_dEnonDef_min_Edef.(run_full_nm)(chosen_trial_idx));
                                case 4
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(trialN_dEnonDef.(run_full_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % confidence
                        if chosenModel_confidence > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'confidence';
                            switch chosenModel_confidence
                                case 1
                                    chosen_modVals(n_chosenMods,:) = confidence.(run_full_nm)(chosen_trial_idx); % binary variable => no zscore
                                case {2,3,4}
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(confidence.(run_full_nm)(chosen_trial_idx)); % confidence inferred by the model
                                otherwise
                                    error('ready yet');
                            end
                        end
                        
                        % RT (last regressor)
                        if chosenModel_RT > 0 && ~ismember(chosenModel_RT,[4,5,6])
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'choice RT';
                            switch chosenModel_RT
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(choice_RT.(run_full_nm)(chosen_trial_idx));
                                case 2
                                    chosen_modVals(n_chosenMods,:) = zscore(choice_RT.(run_full_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        [matlabbatch] = First_level_loadEachCondition_DCM(matlabbatch, sub_idx, iCond,...
                            ['dispChosen_',RP_dispChosen_nm,'_',splitE_dispChosen_nm], modelChosenOnset, modelChosenDur,...
                            n_chosenMods, chosen_modNames, chosen_modVals,...
                            orth_vars, onsets_only_GLM);
                    end % run loop
                case 2 % all tasks independent but sessions pooled
                    for iTask = 1:nTasks
                        task_nm = tasks{iTask};
                        
                        % extract trial index for the current loop
                        switch RP_dispChosen_nm
                            case 'RP'
                                RPfilter_dispChosen = true(1,length(double(RP_var_binary.(task_nm))));
                            case 'R'
                                RPfilter_dispChosen = (RP_var_binary.(task_nm) == 1);
                            case 'P'
                                RPfilter_dispChosen = (RP_var_binary.(task_nm) == 0);
                        end
                        switch splitE_dispChosen_nm
                            case 'E'
                                Efilter_dispChosen = true(1,length(double(RP_var_binary.(task_nm))));
                            case 'E1'
                                Efilter_dispChosen = (E_varOption.(task_nm) == 1);
                            case 'E2'
                                Efilter_dispChosen = (E_varOption.(task_nm) == 2);
                            case 'E3'
                                Efilter_dispChosen = (E_varOption.(task_nm) == 3);
                            case 'Ech0'
                                Efilter_dispChosen = (E_chosen.(task_nm) == 0);
                            case 'Ech1'
                                Efilter_dispChosen = (E_chosen.(task_nm) == 1);
                            case 'Ech2'
                                Efilter_dispChosen = (E_chosen.(task_nm) == 2);
                            case 'Ech3'
                                Efilter_dispChosen = (E_chosen.(task_nm) == 3);
                            case 'lEch'
                                Efilter_dispChosen = (choice_hE.(task_nm) == 0);
                            case 'hEch'
                                Efilter_dispChosen = (choice_hE.(task_nm) == 1);
                        end
                        chosen_trial_idx = (RPfilter_dispChosen.*Efilter_dispChosen) == 1; % NEED to transform it into logical or will just focus on the first trial
                        
                        %% chosen onset
                        iCond = iCond + 1;
                        modelChosenOnset = onsets.dispChosenOnsets.(task_nm)(chosen_trial_idx);
                        % duration
                        switch chosenModel
                            case 'stick'
                                modelChosenDur = 0;
                            case 'boxcar' % duration displaying the chosen option
                                modelChosenDur = durations.dispChosenDur.(task_nm)(chosen_trial_idx);
                            case 'boxcar_bis' % duration going form display of chosen option
                                % until the end of the exertion of the effort
                                modelChosenDur = durations.dispChosenDur.(task_nm)(chosen_trial_idx) +...
                                    durations.preEffortCrossDur.(task_nm)(chosen_trial_idx) +...
                                    durations.EperfDur.(task_nm)(chosen_trial_idx);
                            case 'boxcar_ter' % duration going form display of chosen option
                                % until the end of the effort preparation cross
                                % (entailing the whole effort preparation period)
                                modelChosenDur = durations.dispChosenDur.(task_nm)(chosen_trial_idx) +...
                                    durations.preEffortCrossDur.(task_nm)(chosen_trial_idx);
                        end
                        
                        %% chosen modulators
                        n_chosenMods = 0;
                        chosen_modNames = cell(1,1);
                        chosen_modVals = [];
                        
                        % RT (first regressor)
                        if chosenModel_RT > 0 && ~ismember(chosenModel_RT,[1,2,3])
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'choice RT';
                            switch chosenModel_RT
                                case 4
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(choice_RT.(task_nm)(chosen_trial_idx));
                                case 5
                                    chosen_modVals(n_chosenMods,:) = zscore(choice_RT.(task_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % net value chosen (first regressor)
                        if chosenModel_NV_chosen > 0 && ~ismember(chosenModel_NV_chosen,[1,2,3])
                            n_chosenMods = n_chosenMods + 1;
                            switch chosenModel_NV_chosen
                                case 4
                                    chosen_modNames{n_chosenMods} = 'NVch-NVunch';
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(NV_ch_min_unch.(task_nm)(chosen_trial_idx));
                                case 5
                                    chosen_modNames{n_chosenMods} = 'p(chosen)';
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(pChosen.(task_nm)(chosen_trial_idx));
                                case 6
                                    chosen_modNames{n_chosenMods} = 'NVch-NVunch';
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(NV_ch_min_unch_with_bias.(task_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % reward vs punishment trials
                        if chosenModel_R_vs_P > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'R_vs_P';
                            switch chosenModel_R_vs_P
                                case 1
                                    chosen_modVals(n_chosenMods,:) = RP_var_binary.(task_nm)(chosen_trial_idx);
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % choice = high effort
                        switch chosenModel_choicehE
                            case 0
                            case 1
                                n_chosenMods = n_chosenMods + 1;
                                chosen_modNames{n_chosenMods} = 'choice = high effort';
                                chosen_modVals(n_chosenMods,:) = choice_hE.(task_nm)(chosen_trial_idx); % binary variable => no zscore
                            case 2
                                n_chosenMods = n_chosenMods + 1;
                                chosen_modNames{n_chosenMods} = 'choice = high effort';
                                chosen_modVals(n_chosenMods,:) = choice_hE_bis.(task_nm)(chosen_trial_idx); % binary variable => no zscore
                            otherwise
                                error('not ready yet');
                        end
                        
                        % reward for the non-default option
                        if chosenModel_R_varOption > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'R non-default';
                            switch chosenModel_R_varOption
                                case 1 % R high effort option (amount)
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(R_amount_varOption.(task_nm)(chosen_trial_idx));
                                case 2 % R high effort option (level)
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(R_level_varOption.(task_nm)(chosen_trial_idx));
                                case 3 % z(R high effort option) (amount)
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(z_R_amount_varOption.(task_nm)(chosen_trial_idx));
                                case 4 % z(R high effort option) (level)
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(z_R_level_varOption.(task_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % reward for the chosen option
                        if chosenModel_R_chosen > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'R chosen';
                            switch chosenModel_R_chosen
                                case 1 % R high effort option (amount)
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(R_amount_chosen.(task_nm)(chosen_trial_idx));
                                case 2 % R high effort option (level)
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(R_level_chosen.(task_nm)(chosen_trial_idx));
                                case 3 % z(R high effort option) (amount)
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(z_R_amount_chosen.(task_nm)(chosen_trial_idx));
                                case 4 % z(R high effort option) (level)
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(z_R_level_chosen.(task_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % punishment for the non-default option
                        if chosenModel_P_varOption > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'P non-default';
                            switch chosenModel_P_varOption
                                case 1 % R high effort option (amount)
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(P_amount_varOption.(task_nm)(chosen_trial_idx));
                                case 2 % R high effort option (level)
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(P_level_varOption.(task_nm)(chosen_trial_idx));
                                case 3 % z(R high effort option) (amount)
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(z_P_amount_varOption.(task_nm)(chosen_trial_idx));
                                case 4 % z(R high effort option) (level)
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(z_P_level_varOption.(task_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % punishment for the chosen option
                        if chosenModel_P_chosen > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'P chosen';
                            switch chosenModel_P_chosen
                                case 1 % R high effort option (amount)
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(P_amount_chosen.(task_nm)(chosen_trial_idx));
                                case 2 % R high effort option (level)
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(P_level_chosen.(task_nm)(chosen_trial_idx));
                                case 3 % z(R high effort option) (amount)
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(z_P_amount_chosen.(task_nm)(chosen_trial_idx));
                                case 4 % z(R high effort option) (level)
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(z_P_level_chosen.(task_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % money left
                        if chosenModel_moneyLeft > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'money left';
                            switch chosenModel_moneyLeft
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(money_amount_left.(task_nm)(chosen_trial_idx));
                                case 2
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(abs(money_amount_left.(task_nm)(chosen_trial_idx)));
                                case 3
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(money_level_left.(task_nm)(chosen_trial_idx));
                                case 4
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(abs(money_level_chosen.(task_nm)(chosen_trial_idx)));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % money right
                        if chosenModel_moneyRight > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'money right';
                            switch chosenModel_moneyRight
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(money_amount_right.(task_nm)(chosen_trial_idx));
                                case 2
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(abs(money_amount_right.(task_nm)(chosen_trial_idx)));
                                case 3
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(money_level_right.(task_nm)(chosen_trial_idx));
                                case 4
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(abs(money_level_right.(task_nm)(chosen_trial_idx)));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % money chosen
                        if chosenModel_moneyChosen > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'money chosen';
                            switch chosenModel_moneyChosen
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(money_amount_chosen.(task_nm)(chosen_trial_idx));
                                case 2
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(abs_money_amount_chosen.(task_nm)(chosen_trial_idx));
                                case 3
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(money_level_chosen.(task_nm)(chosen_trial_idx));
                                case 4
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(abs_money_level_chosen.(task_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % money unchosen
                        if chosenModel_moneyUnchosen > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'money unchosen';
                            switch chosenModel_moneyUnchosen
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(money_amount_unchosen.(task_nm)(chosen_trial_idx));
                                case 2
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(abs_money_amount_unchosen.(task_nm)(chosen_trial_idx));
                                case 3
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(money_level_unchosen.(task_nm)(chosen_trial_idx));
                                case 4
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(abs_money_level_unchosen.(task_nm)(chosen_trial_idx));
                                otherwise
                                    error('ready yet');
                            end
                        end
                        
                        % money non-default option
                        if chosenModel_moneyNonDefault > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'money non-default option';
                            switch chosenModel_moneyNonDefault
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(money_amount_varOption.(task_nm)(chosen_trial_idx));
                                case 2
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(abs_money_amount_varOption.(task_nm)(chosen_trial_idx));
                                case 3
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(money_level_varOption.(task_nm)(chosen_trial_idx));
                                case 4
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(abs_money_level_varOption.(task_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % money chosen - money unchosen
                        if chosenModel_money_chosen_min_money_unchosen > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'money chosen - unchosen';
                            switch chosenModel_money_chosen_min_money_unchosen
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(moneyChosen_min_moneyUnchosen_amount.(task_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % money chosen - money fixed option (low R low E)
                        if chosenModel_moneyChosen_min_moneyDefault > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'money chosen - money fixed';
                            switch chosenModel_moneyChosen_min_moneyDefault
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(moneyChosen_min_moneyFixed_amount.(task_nm)(chosen_trial_idx));
                                case 2
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(moneyChosen_min_moneyFixed_level.(task_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % sum of money
                        if chosenModel_moneySum > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'money sum';
                            switch chosenModel_moneySum
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(money_amount_sum.(task_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % effort left
                        if chosenModel_E_left > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'effort left';
                            switch chosenModel_E_left
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(E_left.(task_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % effort right
                        if chosenModel_E_right > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'effort right';
                            switch chosenModel_E_right
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(E_right.(task_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % effort chosen
                        if chosenModel_Echosen > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'effort chosen';
                            switch chosenModel_Echosen
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(E_chosen.(task_nm)(chosen_trial_idx));
                                case 3
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(E_chosen_bis.(task_nm)(chosen_trial_idx));
                                case 4
                                    chosen_modVals(n_chosenMods,:) = zscore(E_chosen.(task_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % effort unchosen
                        if chosenModel_Eunchosen > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'effort unchosen';
                            switch chosenModel_Eunchosen
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(E_unchosen.(task_nm)(chosen_trial_idx));
                                otherwise
                                    error('ready yet');
                            end
                        end
                        
                        % effort non-default option
                        if chosenModel_EnonDefault > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'effort non-default option';
                            switch chosenModel_EnonDefault
                                case 1
                                    chosen_modVals(n_chosenMods,:) = E_varOption.(task_nm)(chosen_trial_idx); % binary variable => no zscore
                                otherwise
                                    error('ready yet');
                            end
                        end
                        
                        % effort chosen - unchosen
                        if chosenModel_Ech_min_Eunch > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'effort chosen - unchosen';
                            switch chosenModel_Ech_min_Eunch
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(E_chosen_min_E_unchosen.(task_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % effort chosen - fixed option
                        if chosenModel_Ech_min_Efixed > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'effort chosen - fixed';
                            switch chosenModel_Ech_min_Efixed
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(Ech_min_Efixed.(task_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % sum of efforts
                        if chosenModel_E_sum > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'effort sum';
                            switch chosenModel_E_sum
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(E_sum.(task_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % money level*effort level high E option
                        if chosenModel_money_level_x_E_varOption > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'money x effort non-default option';
                            switch chosenModel_money_level_x_E_varOption
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(money_level_x_E_varOption.(task_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % money level * effort chosen
                        if chosenModel_money_level_x_E_chosen > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'money x effort chosen';
                            switch chosenModel_money_level_x_E_chosen
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(money_level_x_E_chosen.(task_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % R level * effort high E option
                        if chosenModel_R_level_x_E_varOption > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'R x E non-default option';
                            switch chosenModel_R_level_x_E_varOption
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(R_level_x_E_varOption.(task_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % R level * effort chosen
                        if chosenModel_R_level_x_E_chosen > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'R x E chosen';
                            switch chosenModel_R_level_x_E_chosen
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(R_level_x_E_chosen.(task_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % P level * effort high E option
                        if chosenModel_P_level_x_E_varOption > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'P x E non-default option';
                            switch chosenModel_P_level_x_E_varOption
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(P_level_x_E_varOption.(task_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % P level * effort chosen
                        if chosenModel_P_level_x_E_chosen > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'P x E chosen';
                            switch chosenModel_P_level_x_E_chosen
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(P_level_x_E_chosen.(task_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % net value chosen
                        if chosenModel_NV_chosen > 0 && ~ismember(chosenModel_NV_chosen,[4,5,6])
                            n_chosenMods = n_chosenMods + 1;
                            switch chosenModel_NV_chosen
                                case 1
                                    chosen_modNames{n_chosenMods} = 'NVch-NVunch';
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(NV_ch_min_unch.(task_nm)(chosen_trial_idx));
                                case 2
                                    chosen_modNames{n_chosenMods} = 'p(chosen)';
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(pChosen.(task_nm)(chosen_trial_idx));
                                case 3
                                    chosen_modNames{n_chosenMods} = 'NVch-NVunch';
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(NV_ch_min_unch_with_bias.(task_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % net value non-default option
                        if chosenModel_NV_varOption > 0
                            n_chosenMods = n_chosenMods + 1;
                            switch chosenModel_NV_varOption
                                case 1
                                    chosen_modNames{n_chosenMods} = 'delta NV high E - low E';
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(NV_varOption.(task_nm)(chosen_trial_idx));
                                case 2
                                    chosen_modNames{n_chosenMods} = '|delta NV high E - low E|';
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(abs(NV_varOption.(task_nm)(chosen_trial_idx)));
                                case 3
                                    chosen_modNames{n_chosenMods} = 'p(choice=hE)';
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(abs(pChoice_hE.(task_nm)(chosen_trial_idx)));
                                case 4
                                    chosen_modNames{n_chosenMods} = 'delta NV high E - low E + bias';
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(NV_varOption_plus_bias.(task_nm)(chosen_trial_idx));
                                case 5
                                    chosen_modNames{n_chosenMods} = '|delta NV high E - low E + bias|';
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(abs(NV_varOption_plus_bias.(task_nm)(chosen_trial_idx)));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % net value non-default option bis
                        if chosenModel_NV_varOption_bis > 0
                            n_chosenMods = n_chosenMods + 1;
                            switch chosenModel_NV_varOption_bis
                                case 1
                                    chosen_modNames{n_chosenMods} = 'delta NV high E - low E';
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(NV_varOption.(task_nm)(chosen_trial_idx));
                                case 2
                                    chosen_modNames{n_chosenMods} = '|delta NV high E - low E|';
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(abs(NV_varOption.(task_nm)(chosen_trial_idx)));
                                case 3
                                    chosen_modNames{n_chosenMods} = 'p(choice=hE)';
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(abs(pChoice_hE.(task_nm)(chosen_trial_idx)));
                                case 4
                                    chosen_modNames{n_chosenMods} = 'delta NV high E - low E + bias';
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(NV_varOption_plus_bias.(task_nm)(chosen_trial_idx));
                                case 5
                                    chosen_modNames{n_chosenMods} = '|delta NV high E - low E + bias|';
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(abs(NV_varOption_plus_bias.(task_nm)(chosen_trial_idx)));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        if strcmp(task_id,'Ep')
                            % force integral
                            if chosenModel_F_integral > 0
                                n_chosenMods = n_chosenMods + 1;
                                chosen_modNames{n_chosenMods} = 'F integral';
                                switch chosenModel_F_integral
                                    case 1
                                        chosen_modVals(n_chosenMods,:) = raw_or_z(AUC.(task_nm)(chosen_trial_idx));
                                    case 2
                                        chosen_modVals(n_chosenMods,:) = raw_or_z(AUC_overshoot.(task_nm)(chosen_trial_idx));
                                    case 3
                                        chosen_modVals(n_chosenMods,:) = raw_or_z(AUC_N.(task_nm)(chosen_trial_idx));
                                    case 4
                                        chosen_modVals(n_chosenMods,:) = raw_or_z(AUC_overshoot_N.(task_nm)(chosen_trial_idx));
                                    case 5
                                        chosen_modVals(n_chosenMods,:) = zscore(AUC.(task_nm)(chosen_trial_idx));
                                    case 6
                                        chosen_modVals(n_chosenMods,:) = zscore(AUC_overshoot.(task_nm)(chosen_trial_idx));
                                    case 7
                                        chosen_modVals(n_chosenMods,:) = zscore(AUC_N.(task_nm)(chosen_trial_idx));
                                    case 8
                                        chosen_modVals(n_chosenMods,:) = zscore(AUC_overshoot_N.(task_nm)(chosen_trial_idx));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                            
                            % fatigue
                            if chosenModel_fatigue > 0
                                n_chosenMods = n_chosenMods + 1;
                                chosen_modNames{n_chosenMods} = 'fatigue';
                                switch chosenModel_fatigue
                                    case 1
                                        chosen_modVals(n_chosenMods,:) = raw_or_z(fatigue.(task_nm)(chosen_trial_idx));
                                    case 2
                                        chosen_modVals(n_chosenMods,:) = zscore(fatigue.(task_nm)(chosen_trial_idx));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                            
                            % Ech*fatigue
                            if chosenModel_Ech_x_fatigue > 0
                                n_chosenMods = n_chosenMods + 1;
                                chosen_modNames{n_chosenMods} = 'Ech_x_fatigue';
                                switch chosenModel_Ech_x_fatigue
                                    case 1
                                        chosen_modVals(n_chosenMods,:) = raw_or_z(E_chosen.(task_nm)(chosen_trial_idx).*fatigue.(task_nm)(chosen_trial_idx));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                        end % physical effort filter
                        
                        if strcmp(task_id,'Em')
                            % efficacy
                            if chosenModel_efficacy > 0
                                n_chosenMods = n_chosenMods + 1;
                                chosen_modNames{n_chosenMods} = 'efficacy';
                                switch chosenModel_efficacy
                                    case 1
                                        chosen_modVals(n_chosenMods,:) = raw_or_z(efficacy_with2first.(task_nm)(chosen_trial_idx));
                                    case 2
                                        chosen_modVals(n_chosenMods,:) = raw_or_z(efficacy_pureNback.(task_nm)(chosen_trial_idx));
                                    case 3
                                        chosen_modVals(n_chosenMods,:) = raw_or_z(efficacy_bis_with2first.(task_nm)(chosen_trial_idx));
                                    case 4
                                        chosen_modVals(n_chosenMods,:) = raw_or_z(efficacy_bis_pureNback.(task_nm)(chosen_trial_idx));
                                    case 5
                                        chosen_modVals(n_chosenMods,:) = zscore(efficacy_with2first.(task_nm)(chosen_trial_idx));
                                    case 6
                                        chosen_modVals(n_chosenMods,:) = zscore(efficacy_pureNback.(task_nm)(chosen_trial_idx));
                                    case 7
                                        chosen_modVals(n_chosenMods,:) = zscore(efficacy_bis_with2first.(task_nm)(chosen_trial_idx));
                                    case 8
                                        chosen_modVals(n_chosenMods,:) = zscore(efficacy_bis_pureNback.(task_nm)(chosen_trial_idx));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                            
                            % previous trial efficacy
                            if chosenModel_prevEfficacy > 0
                                n_chosenMods = n_chosenMods + 1;
                                chosen_modNames{n_chosenMods} = 'previous trial efficacy';
                                switch chosenModel_prevEfficacy
                                    case 1
                                        chosen_modVals(n_chosenMods,:) = raw_or_z(prevEfficacy_with2first.(task_nm)(chosen_trial_idx));
                                    case 2
                                        chosen_modVals(n_chosenMods,:) = raw_or_z(prevEfficacy_pureNback.(task_nm)(chosen_trial_idx));
                                    case 3
                                        chosen_modVals(n_chosenMods,:) = raw_or_z(prevEfficacy_bis_with2first.(task_nm)(chosen_trial_idx));
                                    case 4
                                        chosen_modVals(n_chosenMods,:) = raw_or_z(prevEfficacy_bis_pureNback.(task_nm)(chosen_trial_idx));
                                    case 5
                                        chosen_modVals(n_chosenMods,:) = zscore(prevEfficacy_with2first.(task_nm)(chosen_trial_idx));
                                    case 6
                                        chosen_modVals(n_chosenMods,:) = zscore(prevEfficacy_pureNback.(task_nm)(chosen_trial_idx));
                                    case 7
                                        chosen_modVals(n_chosenMods,:) = zscore(prevEfficacy_bis_with2first.(task_nm)(chosen_trial_idx));
                                    case 8
                                        chosen_modVals(n_chosenMods,:) = zscore(prevEfficacy_bis_pureNback.(task_nm)(chosen_trial_idx));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                            
                            % Ech*(previous trial efficacy)
                            if chosenModel_Ech_x_prevEfficacy > 0
                                n_chosenMods = n_chosenMods + 1;
                                chosen_modNames{n_chosenMods} = 'Ech_x_previous trial efficacy';
                                switch chosenModel_Ech_x_prevEfficacy
                                    case 1
                                        chosen_modVals(n_chosenMods,:) = raw_or_z(Echosen.(task_nm)(chosen_trial_idx).*prevEfficacy_with2first.(task_nm)(chosen_trial_idx));
                                    case 2
                                        chosen_modVals(n_chosenMods,:) = raw_or_z(Echosen.(task_nm)(chosen_trial_idx).*prevEfficacy_pureNback.(task_nm)(chosen_trial_idx));
                                    case 3
                                        chosen_modVals(n_chosenMods,:) = raw_or_z(Echosen.(task_nm)(chosen_trial_idx).*prevEfficacy_bis_with2first.(task_nm)(chosen_trial_idx));
                                    case 4
                                        chosen_modVals(n_chosenMods,:) = raw_or_z(Echosen.(task_nm)(chosen_trial_idx).*prevEfficacy_bis_pureNback.(task_nm)(chosen_trial_idx));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                        end % mental effort filter
                        
                        % trial number
                        if chosenModel_trialN > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'trial number';
                            switch chosenModel_trialN
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(trialN.(task_nm)(chosen_trial_idx));
                                case 2
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(trialN_dEch.(task_nm)(chosen_trial_idx));
                                case 3
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(trialN_dEnonDef_min_Edef.(task_nm)(chosen_trial_idx));
                                case 4
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(trialN_dEnonDef.(task_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % confidence
                        if chosenModel_confidence > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'confidence';
                            switch chosenModel_confidence
                                case 1
                                    chosen_modVals(n_chosenMods,:) = confidence.(task_nm)(chosen_trial_idx); % binary variable => no zscore
                                case {2,3,4}
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(confidence.(task_nm)(chosen_trial_idx)); % confidence inferred by the model
                                otherwise
                                    error('ready yet');
                            end
                        end
                        
                        % RT (last regressor)
                        if chosenModel_RT > 0 && ~ismember(chosenModel_RT,[4,5,6])
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'choice RT';
                            switch chosenModel_RT
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(choice_RT.(task_nm)(chosen_trial_idx));
                                case 2
                                    chosen_modVals(n_chosenMods,:) = zscore(choice_RT.(task_nm)(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        [matlabbatch] = First_level_loadEachCondition_DCM(matlabbatch, sub_idx, iCond,...
                            ['dispChosen_',RP_dispChosen_nm,'_',splitE_dispChosen_nm], modelChosenOnset, modelChosenDur,...
                            n_chosenMods, chosen_modNames, chosen_modVals,...
                            orth_vars, onsets_only_GLM);
                    end % task loop
                case {3,4,5} % all fixation crosses pooled across sessions
                    
                    % extract trial index for the current loop
                    switch RP_dispChosen_nm
                        case 'RP'
                            RPfilter_dispChosen = true(1,length(double(RP_var_binary.allTrials)));
                        case 'R'
                            RPfilter_dispChosen = (RP_var_binary.allTrials == 1);
                        case 'P'
                            RPfilter_dispChosen = (RP_var_binary.allTrials == 0);
                    end
                    switch splitE_dispChosen_nm
                        case 'E'
                            Efilter_dispChosen = true(1,length(double(RP_var_binary.allTrials)));
                        case 'E1'
                            Efilter_dispChosen = (E_varOption.allTrials == 1);
                        case 'E2'
                            Efilter_dispChosen = (E_varOption.allTrials == 2);
                        case 'E3'
                            Efilter_dispChosen = (E_varOption.allTrials == 3);
                        case 'Ech0'
                            Efilter_dispChosen = (E_chosen.allTrials == 0);
                        case 'Ech1'
                            Efilter_dispChosen = (E_chosen.allTrials == 1);
                        case 'Ech2'
                            Efilter_dispChosen = (E_chosen.allTrials == 2);
                        case 'Ech3'
                            Efilter_dispChosen = (E_chosen.allTrials == 3);
                        case 'lEch'
                            Efilter_dispChosen = (choice_hE.allTrials == 0);
                        case 'hEch'
                            Efilter_dispChosen = (choice_hE.allTrials == 1);
                    end
                    chosen_trial_idx = (RPfilter_dispChosen.*Efilter_dispChosen) == 1; % NEED to transform it into logical or will just focus on the first trial
                    
                    %% chosen onset
                    iCond = iCond + 1;
                    modelChosenOnset = onsets.dispChosenOnsets.allTrials(chosen_trial_idx);
                    % duration
                    switch chosenModel
                        case 'stick'
                            modelChosenDur = 0;
                        case 'boxcar' % duration displaying the chosen option
                            modelChosenDur = durations.dispChosenDur.allTrials(chosen_trial_idx);
                        case 'boxcar_bis' % duration going form display of chosen option
                            % until the end of the exertion of the effort
                            modelChosenDur = durations.dispChosenDur.allTrials(chosen_trial_idx) +...
                                durations.preEffortCrossDur.allTrials(chosen_trial_idx) +...
                                durations.EperfDur.allTrials(chosen_trial_idx);
                        case 'boxcar_ter' % duration going form display of chosen option
                            % until the end of the effort preparation cross
                            % (entailing the whole effort preparation period)
                            modelChosenDur = durations.dispChosenDur.allTrials(chosen_trial_idx) +...
                                durations.preEffortCrossDur.allTrials(chosen_trial_idx);
                    end
                    
                    %% chosen modulators
                    n_chosenMods = 0;
                    chosen_modNames = cell(1,1);
                    chosen_modVals = [];
                    
                    % RT (first regressor)
                    if chosenModel_RT > 0 && ~ismember(chosenModel_RT,[1,2,3])
                        n_chosenMods = n_chosenMods + 1;
                        chosen_modNames{n_chosenMods} = 'choice RT';
                        switch chosenModel_RT
                            case 4
                                chosen_modVals(n_chosenMods,:) = raw_or_z(choice_RT.allTrials(chosen_trial_idx));
                            case 5
                                chosen_modVals(n_chosenMods,:) = zscore(choice_RT.allTrials(chosen_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % net value chosen (first regressor)
                    if chosenModel_NV_chosen > 0 && ~ismember(chosenModel_NV_chosen,[1,2,3])
                        n_chosenMods = n_chosenMods + 1;
                        switch chosenModel_NV_chosen
                            case 4
                                chosen_modNames{n_chosenMods} = 'NVch-NVunch';
                                chosen_modVals(n_chosenMods,:) = raw_or_z(NV_ch_min_unch.allTrials(chosen_trial_idx));
                            case 5
                                chosen_modNames{n_chosenMods} = 'p(chosen)';
                                chosen_modVals(n_chosenMods,:) = raw_or_z(pChosen.allTrials(chosen_trial_idx));
                            case 6
                                chosen_modNames{n_chosenMods} = 'NVch-NVunch';
                                chosen_modVals(n_chosenMods,:) = raw_or_z(NV_ch_min_unch_with_bias.allTrials(chosen_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % reward vs punishment trials
                    if chosenModel_R_vs_P > 0
                        n_chosenMods = n_chosenMods + 1;
                        chosen_modNames{n_chosenMods} = 'R_vs_P';
                        switch chosenModel_R_vs_P
                            case 1
                                chosen_modVals(n_chosenMods,:) = RP_var_binary.allTrials(chosen_trial_idx);
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % choice = high effort
                    switch chosenModel_choicehE
                        case 0
                        case 1
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'choice = high effort';
                            chosen_modVals(n_chosenMods,:) = choice_hE.allTrials(chosen_trial_idx); % binary variable => no zscore
                        case 2
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'choice = high effort';
                            chosen_modVals(n_chosenMods,:) = choice_hE_bis.allTrials(chosen_trial_idx); % binary variable => no zscore
                        otherwise
                            error('not ready yet');
                    end
                    
                    % reward for the non-default option
                    if chosenModel_R_varOption > 0
                        n_chosenMods = n_chosenMods + 1;
                        chosen_modNames{n_chosenMods} = 'R non-default';
                        switch chosenModel_R_varOption
                            case 1 % R high effort option (amount)
                                chosen_modVals(n_chosenMods,:) = raw_or_z(R_amount_varOption.allTrials(chosen_trial_idx));
                            case 2 % R high effort option (level)
                                chosen_modVals(n_chosenMods,:) = raw_or_z(R_level_varOption.allTrials(chosen_trial_idx));
                            case 3 % z(R high effort option) (amount)
                                chosen_modVals(n_chosenMods,:) = raw_or_z(z_R_amount_varOption.allTrials(chosen_trial_idx));
                            case 4 % z(R high effort option) (level)
                                chosen_modVals(n_chosenMods,:) = raw_or_z(z_R_level_varOption.allTrials(chosen_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % reward for the chosen option
                    if chosenModel_R_chosen > 0
                        n_chosenMods = n_chosenMods + 1;
                        chosen_modNames{n_chosenMods} = 'R chosen';
                        switch chosenModel_R_chosen
                            case 1 % R high effort option (amount)
                                chosen_modVals(n_chosenMods,:) = raw_or_z(R_amount_chosen.allTrials(chosen_trial_idx));
                            case 2 % R high effort option (level)
                                chosen_modVals(n_chosenMods,:) = raw_or_z(R_level_chosen.allTrials(chosen_trial_idx));
                            case 3 % z(R high effort option) (amount)
                                chosen_modVals(n_chosenMods,:) = raw_or_z(z_R_amount_chosen.allTrials(chosen_trial_idx));
                            case 4 % z(R high effort option) (level)
                                chosen_modVals(n_chosenMods,:) = raw_or_z(z_R_level_chosen.allTrials(chosen_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % punishment for the non-default option
                    if chosenModel_P_varOption > 0
                        n_chosenMods = n_chosenMods + 1;
                        chosen_modNames{n_chosenMods} = 'P non-default';
                        switch chosenModel_P_varOption
                            case 1 % R high effort option (amount)
                                chosen_modVals(n_chosenMods,:) = raw_or_z(P_amount_varOption.allTrials(chosen_trial_idx));
                            case 2 % R high effort option (level)
                                chosen_modVals(n_chosenMods,:) = raw_or_z(P_level_varOption.allTrials(chosen_trial_idx));
                            case 3 % z(R high effort option) (amount)
                                chosen_modVals(n_chosenMods,:) = raw_or_z(z_P_amount_varOption.allTrials(chosen_trial_idx));
                            case 4 % z(R high effort option) (level)
                                chosen_modVals(n_chosenMods,:) = raw_or_z(z_P_level_varOption.allTrials(chosen_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % punishment for the chosen option
                    if chosenModel_P_chosen > 0
                        n_chosenMods = n_chosenMods + 1;
                        chosen_modNames{n_chosenMods} = 'P chosen';
                        switch chosenModel_P_chosen
                            case 1 % R high effort option (amount)
                                chosen_modVals(n_chosenMods,:) = raw_or_z(P_amount_chosen.allTrials(chosen_trial_idx));
                            case 2 % R high effort option (level)
                                chosen_modVals(n_chosenMods,:) = raw_or_z(P_level_chosen.allTrials(chosen_trial_idx));
                            case 3 % z(R high effort option) (amount)
                                chosen_modVals(n_chosenMods,:) = raw_or_z(z_P_amount_chosen.allTrials(chosen_trial_idx));
                            case 4 % z(R high effort option) (level)
                                chosen_modVals(n_chosenMods,:) = raw_or_z(z_P_level_chosen.allTrials(chosen_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % money left
                    if chosenModel_moneyLeft > 0
                        n_chosenMods = n_chosenMods + 1;
                        chosen_modNames{n_chosenMods} = 'money left';
                        switch chosenModel_moneyLeft
                            case 1
                                chosen_modVals(n_chosenMods,:) = raw_or_z(money_amount_left.allTrials(chosen_trial_idx));
                            case 2
                                chosen_modVals(n_chosenMods,:) = raw_or_z(abs(money_amount_left.allTrials(chosen_trial_idx)));
                            case 3
                                chosen_modVals(n_chosenMods,:) = raw_or_z(money_level_left.allTrials(chosen_trial_idx));
                            case 4
                                chosen_modVals(n_chosenMods,:) = raw_or_z(abs(money_level_chosen.allTrials(chosen_trial_idx)));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % money right
                    if chosenModel_moneyRight > 0
                        n_chosenMods = n_chosenMods + 1;
                        chosen_modNames{n_chosenMods} = 'money right';
                        switch chosenModel_moneyRight
                            case 1
                                chosen_modVals(n_chosenMods,:) = raw_or_z(money_amount_right.allTrials(chosen_trial_idx));
                            case 2
                                chosen_modVals(n_chosenMods,:) = raw_or_z(abs(money_amount_right.allTrials(chosen_trial_idx)));
                            case 3
                                chosen_modVals(n_chosenMods,:) = raw_or_z(money_level_right.allTrials(chosen_trial_idx));
                            case 4
                                chosen_modVals(n_chosenMods,:) = raw_or_z(abs(money_level_right.allTrials(chosen_trial_idx)));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % money chosen
                    if chosenModel_moneyChosen > 0
                        n_chosenMods = n_chosenMods + 1;
                        chosen_modNames{n_chosenMods} = 'money chosen';
                        switch chosenModel_moneyChosen
                            case 1
                                chosen_modVals(n_chosenMods,:) = raw_or_z(money_amount_chosen.allTrials(chosen_trial_idx));
                            case 2
                                chosen_modVals(n_chosenMods,:) = raw_or_z(abs_money_amount_chosen.allTrials(chosen_trial_idx));
                            case 3
                                chosen_modVals(n_chosenMods,:) = raw_or_z(money_level_chosen.allTrials(chosen_trial_idx));
                            case 4
                                chosen_modVals(n_chosenMods,:) = raw_or_z(abs_money_level_chosen.allTrials(chosen_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % money unchosen
                    if chosenModel_moneyUnchosen > 0
                        n_chosenMods = n_chosenMods + 1;
                        chosen_modNames{n_chosenMods} = 'money unchosen';
                        switch chosenModel_moneyUnchosen
                            case 1
                                chosen_modVals(n_chosenMods,:) = raw_or_z(money_amount_unchosen.allTrials(chosen_trial_idx));
                            case 2
                                chosen_modVals(n_chosenMods,:) = raw_or_z(abs_money_amount_unchosen.allTrials(chosen_trial_idx));
                            case 3
                                chosen_modVals(n_chosenMods,:) = raw_or_z(money_level_unchosen.allTrials(chosen_trial_idx));
                            case 4
                                chosen_modVals(n_chosenMods,:) = raw_or_z(abs_money_level_unchosen.allTrials(chosen_trial_idx));
                            otherwise
                                error('ready yet');
                        end
                    end
                    
                    % money non-default option
                    if chosenModel_moneyNonDefault > 0
                        n_chosenMods = n_chosenMods + 1;
                        chosen_modNames{n_chosenMods} = 'money non-default option';
                        switch chosenModel_moneyNonDefault
                            case 1
                                chosen_modVals(n_chosenMods,:) = raw_or_z(money_amount_varOption.allTrials(chosen_trial_idx));
                            case 2
                                chosen_modVals(n_chosenMods,:) = raw_or_z(abs_money_amount_varOption.allTrials(chosen_trial_idx));
                            case 3
                                chosen_modVals(n_chosenMods,:) = raw_or_z(money_level_varOption.allTrials(chosen_trial_idx));
                            case 4
                                chosen_modVals(n_chosenMods,:) = raw_or_z(abs_money_level_varOption.allTrials(chosen_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % money chosen - money unchosen
                    if chosenModel_money_chosen_min_money_unchosen > 0
                        n_chosenMods = n_chosenMods + 1;
                        chosen_modNames{n_chosenMods} = 'money chosen - unchosen';
                        switch chosenModel_money_chosen_min_money_unchosen
                            case 1
                                chosen_modVals(n_chosenMods,:) = raw_or_z(moneyChosen_min_moneyUnchosen_amount.allTrials(chosen_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % money chosen - money fixed option (low R low E)
                    if chosenModel_moneyChosen_min_moneyDefault > 0
                        n_chosenMods = n_chosenMods + 1;
                        chosen_modNames{n_chosenMods} = 'money chosen - money fixed';
                        switch chosenModel_moneyChosen_min_moneyDefault
                            case 1
                                chosen_modVals(n_chosenMods,:) = raw_or_z(moneyChosen_min_moneyFixed_amount.allTrials(chosen_trial_idx));
                            case 2
                                chosen_modVals(n_chosenMods,:) = raw_or_z(moneyChosen_min_moneyFixed_level.allTrials(chosen_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % sum of money
                    if chosenModel_moneySum > 0
                        n_chosenMods = n_chosenMods + 1;
                        chosen_modNames{n_chosenMods} = 'money sum';
                        switch chosenModel_moneySum
                            case 1
                                chosen_modVals(n_chosenMods,:) = raw_or_z(money_amount_sum.allTrials(chosen_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % effort left
                    if chosenModel_E_left > 0
                        n_chosenMods = n_chosenMods + 1;
                        chosen_modNames{n_chosenMods} = 'effort left';
                        switch chosenModel_E_left
                            case 1
                                chosen_modVals(n_chosenMods,:) = raw_or_z(E_left.allTrials(chosen_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % effort right
                    if chosenModel_E_right > 0
                        n_chosenMods = n_chosenMods + 1;
                        chosen_modNames{n_chosenMods} = 'effort right';
                        switch chosenModel_E_right
                            case 1
                                chosen_modVals(n_chosenMods,:) = raw_or_z(E_right.allTrials(chosen_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % effort chosen
                    if chosenModel_Echosen > 0
                        n_chosenMods = n_chosenMods + 1;
                        chosen_modNames{n_chosenMods} = 'effort chosen';
                        switch chosenModel_Echosen
                            case 1
                                chosen_modVals(n_chosenMods,:) = raw_or_z(E_chosen.allTrials(chosen_trial_idx));
                            case 3
                                chosen_modVals(n_chosenMods,:) = raw_or_z(E_chosen_bis.allTrials(chosen_trial_idx));
                            case 4
                                chosen_modVals(n_chosenMods,:) = zscore(E_chosen.allTrials(chosen_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % effort unchosen
                    if chosenModel_Eunchosen > 0
                        n_chosenMods = n_chosenMods + 1;
                        chosen_modNames{n_chosenMods} = 'effort unchosen';
                        switch chosenModel_Eunchosen
                            case 1
                                chosen_modVals(n_chosenMods,:) = raw_or_z(E_unchosen.allTrials(chosen_trial_idx));
                            otherwise
                                error('ready yet');
                        end
                    end
                    
                    % effort non-default option
                    if chosenModel_EnonDefault > 0
                        n_chosenMods = n_chosenMods + 1;
                        chosen_modNames{n_chosenMods} = 'effort non-default option';
                        switch chosenModel_EnonDefault
                            case 1
                                chosen_modVals(n_chosenMods,:) = E_varOption.allTrials(chosen_trial_idx); % binary variable => no zscore
                            otherwise
                                error('ready yet');
                        end
                    end
                    
                    % effort chosen - unchosen
                    if chosenModel_Ech_min_Eunch > 0
                        n_chosenMods = n_chosenMods + 1;
                        chosen_modNames{n_chosenMods} = 'effort chosen - unchosen';
                        switch chosenModel_Ech_min_Eunch
                            case 1
                                chosen_modVals(n_chosenMods,:) = raw_or_z(E_chosen_min_E_unchosen.allTrials(chosen_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % effort chosen - fixed option
                    if chosenModel_Ech_min_Efixed > 0
                        n_chosenMods = n_chosenMods + 1;
                        chosen_modNames{n_chosenMods} = 'effort chosen - fixed';
                        switch chosenModel_Ech_min_Efixed
                            case 1
                                chosen_modVals(n_chosenMods,:) = raw_or_z(Ech_min_Efixed.allTrials(chosen_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % sum of efforts
                    if chosenModel_E_sum > 0
                        n_chosenMods = n_chosenMods + 1;
                        chosen_modNames{n_chosenMods} = 'effort sum';
                        switch chosenModel_E_sum
                            case 1
                                chosen_modVals(n_chosenMods,:) = raw_or_z(E_sum.allTrials(chosen_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % money level*effort level high E option
                    if chosenModel_money_level_x_E_varOption > 0
                        n_chosenMods = n_chosenMods + 1;
                        chosen_modNames{n_chosenMods} = 'money x effort non-default option';
                        switch chosenModel_money_level_x_E_varOption
                            case 1
                                chosen_modVals(n_chosenMods,:) = raw_or_z(money_level_x_E_varOption.allTrials(chosen_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % money level * effort chosen
                    if chosenModel_money_level_x_E_chosen > 0
                        n_chosenMods = n_chosenMods + 1;
                        chosen_modNames{n_chosenMods} = 'money x effort chosen';
                        switch chosenModel_money_level_x_E_chosen
                            case 1
                                chosen_modVals(n_chosenMods,:) = raw_or_z(money_level_x_E_chosen.allTrials(chosen_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % R level * effort high E option
                    if chosenModel_R_level_x_E_varOption > 0
                        n_chosenMods = n_chosenMods + 1;
                        chosen_modNames{n_chosenMods} = 'R x E non-default option';
                        switch chosenModel_R_level_x_E_varOption
                            case 1
                                chosen_modVals(n_chosenMods,:) = raw_or_z(R_level_x_E_varOption.allTrials(chosen_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % R level * effort chosen
                    if chosenModel_R_level_x_E_chosen > 0
                        n_chosenMods = n_chosenMods + 1;
                        chosen_modNames{n_chosenMods} = 'R x E chosen';
                        switch chosenModel_R_level_x_E_chosen
                            case 1
                                chosen_modVals(n_chosenMods,:) = raw_or_z(R_level_x_E_chosen.allTrials(chosen_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % P level * effort high E option
                    if chosenModel_P_level_x_E_varOption > 0
                        n_chosenMods = n_chosenMods + 1;
                        chosen_modNames{n_chosenMods} = 'P x E non-default option';
                        switch chosenModel_P_level_x_E_varOption
                            case 1
                                chosen_modVals(n_chosenMods,:) = raw_or_z(P_level_x_E_varOption.allTrials(chosen_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % P level * effort chosen
                    if chosenModel_P_level_x_E_chosen > 0
                        n_chosenMods = n_chosenMods + 1;
                        chosen_modNames{n_chosenMods} = 'P x E chosen';
                        switch chosenModel_P_level_x_E_chosen
                            case 1
                                chosen_modVals(n_chosenMods,:) = raw_or_z(P_level_x_E_chosen.allTrials(chosen_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % net value chosen
                    if chosenModel_NV_chosen > 0 && ~ismember(chosenModel_NV_chosen,[4,5,6])
                        n_chosenMods = n_chosenMods + 1;
                        switch chosenModel_NV_chosen
                            case 1
                                chosen_modNames{n_chosenMods} = 'NVch-NVunch';
                                chosen_modVals(n_chosenMods,:) = raw_or_z(NV_ch_min_unch.allTrials(chosen_trial_idx));
                            case 2
                                chosen_modNames{n_chosenMods} = 'p(chosen)';
                                chosen_modVals(n_chosenMods,:) = raw_or_z(pChosen.allTrials(chosen_trial_idx));
                            case 3
                                chosen_modNames{n_chosenMods} = 'NVch-NVunch';
                                chosen_modVals(n_chosenMods,:) = raw_or_z(NV_ch_min_unch_with_bias.allTrials(chosen_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % net value non-default option
                    if chosenModel_NV_varOption > 0
                        n_chosenMods = n_chosenMods + 1;
                        switch chosenModel_NV_varOption
                            case 1
                                chosen_modNames{n_chosenMods} = 'delta NV high E - low E';
                                chosen_modVals(n_chosenMods,:) = raw_or_z(NV_varOption.allTrials(chosen_trial_idx));
                            case 2
                                chosen_modNames{n_chosenMods} = '|delta NV high E - low E|';
                                chosen_modVals(n_chosenMods,:) = raw_or_z(abs(NV_varOption.allTrials(chosen_trial_idx)));
                            case 3
                                chosen_modNames{n_chosenMods} = 'p(choice=hE)';
                                chosen_modVals(n_chosenMods,:) = raw_or_z(abs(pChoice_hE.allTrials(chosen_trial_idx)));
                            case 4
                                chosen_modNames{n_chosenMods} = 'delta NV high E - low E + bias';
                                chosen_modVals(n_chosenMods,:) = raw_or_z(NV_varOption_plus_bias.allTrials(chosen_trial_idx));
                            case 5
                                chosen_modNames{n_chosenMods} = '|delta NV high E - low E + bias|';
                                chosen_modVals(n_chosenMods,:) = raw_or_z(abs(NV_varOption_plus_bias.allTrials(chosen_trial_idx)));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % net value non-default option bis
                    if chosenModel_NV_varOption_bis > 0
                        n_chosenMods = n_chosenMods + 1;
                        switch chosenModel_NV_varOption_bis
                            case 1
                                chosen_modNames{n_chosenMods} = 'delta NV high E - low E';
                                chosen_modVals(n_chosenMods,:) = raw_or_z(NV_varOption.allTrials(chosen_trial_idx));
                            case 2
                                chosen_modNames{n_chosenMods} = '|delta NV high E - low E|';
                                chosen_modVals(n_chosenMods,:) = raw_or_z(abs(NV_varOption.allTrials(chosen_trial_idx)));
                            case 3
                                chosen_modNames{n_chosenMods} = 'p(choice=hE)';
                                chosen_modVals(n_chosenMods,:) = raw_or_z(abs(pChoice_hE.allTrials(chosen_trial_idx)));
                            case 4
                                chosen_modNames{n_chosenMods} = 'delta NV high E - low E + bias';
                                chosen_modVals(n_chosenMods,:) = raw_or_z(NV_varOption_plus_bias.allTrials(chosen_trial_idx));
                            case 5
                                chosen_modNames{n_chosenMods} = '|delta NV high E - low E + bias|';
                                chosen_modVals(n_chosenMods,:) = raw_or_z(abs(NV_varOption_plus_bias.allTrials(chosen_trial_idx)));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    if strcmp(task_id,'Ep')
                        % force integral
                        if chosenModel_F_integral > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'F integral';
                            switch chosenModel_F_integral
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(AUC.allTrials(chosen_trial_idx));
                                case 2
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(AUC_overshoot.allTrials(chosen_trial_idx));
                                case 3
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(AUC_N.allTrials(chosen_trial_idx));
                                case 4
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(AUC_overshoot_N.allTrials(chosen_trial_idx));
                                case 5
                                    chosen_modVals(n_chosenMods,:) = zscore(AUC.allTrials(chosen_trial_idx));
                                case 6
                                    chosen_modVals(n_chosenMods,:) = zscore(AUC_overshoot.allTrials(chosen_trial_idx));
                                case 7
                                    chosen_modVals(n_chosenMods,:) = zscore(AUC_N.allTrials(chosen_trial_idx));
                                case 8
                                    chosen_modVals(n_chosenMods,:) = zscore(AUC_overshoot_N.allTrials(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % fatigue
                        if chosenModel_fatigue > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'fatigue';
                            switch chosenModel_fatigue
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(fatigue.allTrials(chosen_trial_idx));
                                case 2
                                    chosen_modVals(n_chosenMods,:) = zscore(fatigue.allTrials(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % Ech*fatigue
                        if chosenModel_Ech_x_fatigue > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'Ech_x_fatigue';
                            switch chosenModel_Ech_x_fatigue
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(E_chosen.allTrials(chosen_trial_idx).*fatigue.allTrials(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                    end % physical effort filter
                    
                    if strcmp(task_id,'Em')
                        % efficacy
                        if chosenModel_efficacy > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'efficacy';
                            switch chosenModel_efficacy
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(efficacy_with2first.allTrials(chosen_trial_idx));
                                case 2
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(efficacy_pureNback.allTrials(chosen_trial_idx));
                                case 3
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(efficacy_bis_with2first.allTrials(chosen_trial_idx));
                                case 4
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(efficacy_bis_pureNback.allTrials(chosen_trial_idx));
                                case 5
                                    chosen_modVals(n_chosenMods,:) = zscore(efficacy_with2first.allTrials(chosen_trial_idx));
                                case 6
                                    chosen_modVals(n_chosenMods,:) = zscore(efficacy_pureNback.allTrials(chosen_trial_idx));
                                case 7
                                    chosen_modVals(n_chosenMods,:) = zscore(efficacy_bis_with2first.allTrials(chosen_trial_idx));
                                case 8
                                    chosen_modVals(n_chosenMods,:) = zscore(efficacy_bis_pureNback.allTrials(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % previous trial efficacy
                        if chosenModel_prevEfficacy > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'previous trial efficacy';
                            switch chosenModel_prevEfficacy
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(prevEfficacy_with2first.allTrials(chosen_trial_idx));
                                case 2
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(prevEfficacy_pureNback.allTrials(chosen_trial_idx));
                                case 3
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(prevEfficacy_bis_with2first.allTrials(chosen_trial_idx));
                                case 4
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(prevEfficacy_bis_pureNback.allTrials(chosen_trial_idx));
                                case 5
                                    chosen_modVals(n_chosenMods,:) = zscore(prevEfficacy_with2first.allTrials(chosen_trial_idx));
                                case 6
                                    chosen_modVals(n_chosenMods,:) = zscore(prevEfficacy_pureNback.allTrials(chosen_trial_idx));
                                case 7
                                    chosen_modVals(n_chosenMods,:) = zscore(prevEfficacy_bis_with2first.allTrials(chosen_trial_idx));
                                case 8
                                    chosen_modVals(n_chosenMods,:) = zscore(prevEfficacy_bis_pureNback.allTrials(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % Ech*(previous trial efficacy)
                        if chosenModel_Ech_x_prevEfficacy > 0
                            n_chosenMods = n_chosenMods + 1;
                            chosen_modNames{n_chosenMods} = 'Ech_x_previous trial efficacy';
                            switch chosenModel_Ech_x_prevEfficacy
                                case 1
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(Echosen.allTrials(chosen_trial_idx).*prevEfficacy_with2first.allTrials(chosen_trial_idx));
                                case 2
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(Echosen.allTrials(chosen_trial_idx).*prevEfficacy_pureNback.allTrials(chosen_trial_idx));
                                case 3
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(Echosen.allTrials(chosen_trial_idx).*prevEfficacy_bis_with2first.allTrials(chosen_trial_idx));
                                case 4
                                    chosen_modVals(n_chosenMods,:) = raw_or_z(Echosen.allTrials(chosen_trial_idx).*prevEfficacy_bis_pureNback.allTrials(chosen_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                    end % mental effort filter
                    
                    % trial number
                    if chosenModel_trialN > 0
                        n_chosenMods = n_chosenMods + 1;
                        chosen_modNames{n_chosenMods} = 'trial number';
                        switch chosenModel_trialN
                            case 1
                                chosen_modVals(n_chosenMods,:) = raw_or_z(trialN.allTrials(chosen_trial_idx));
                            case 2
                                chosen_modVals(n_chosenMods,:) = raw_or_z(trialN_dEch.allTrials(chosen_trial_idx));
                            case 3
                                chosen_modVals(n_chosenMods,:) = raw_or_z(trialN_dEnonDef_min_Edef.allTrials(chosen_trial_idx));
                            case 4
                                chosen_modVals(n_chosenMods,:) = raw_or_z(trialN_dEnonDef.allTrials(chosen_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % confidence
                    if chosenModel_confidence > 0
                        n_chosenMods = n_chosenMods + 1;
                        chosen_modNames{n_chosenMods} = 'confidence';
                        switch chosenModel_confidence
                            case 1
                                chosen_modVals(n_chosenMods,:) = confidence.allTrials(chosen_trial_idx); % binary variable => no zscore
                            case {2,3,4}
                                chosen_modVals(n_chosenMods,:) = raw_or_z(confidence.allTrials(chosen_trial_idx)); % confidence inferred by the model
                            otherwise
                                error('ready yet');
                        end
                    end
                    
                    % RT (last regressor)
                    if chosenModel_RT > 0 && ~ismember(chosenModel_RT,[4,5,6])
                        n_chosenMods = n_chosenMods + 1;
                        chosen_modNames{n_chosenMods} = 'choice RT';
                        switch chosenModel_RT
                            case 1
                                chosen_modVals(n_chosenMods,:) = raw_or_z(choice_RT.allTrials(chosen_trial_idx));
                            case 2
                                chosen_modVals(n_chosenMods,:) = zscore(choice_RT.allTrials(chosen_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    [matlabbatch] = First_level_loadEachCondition_DCM(matlabbatch, sub_idx, iCond,...
                        ['dispChosen_',RP_dispChosen_nm,'_',splitE_dispChosen_nm], modelChosenOnset, modelChosenDur,...
                        n_chosenMods, chosen_modNames, chosen_modVals,...
                        orth_vars, onsets_only_GLM);
            end % DCM mode: info about how regressors concatenated or not across tasks or sessions
        end % E condition
    end % RP
end % model chosen period

%% fixation cross before effort
preEffortCrossModel = GLMprm.model_onset.(task_id).preEffortCross;
if ismember(preEffortCrossModel,{'stick','boxcar','boxcar_bis'})
    
    for iRP_preEcross = 1:length(RPpreEcrossCond)
        RP_preEcross_nm = RPpreEcrossCond{iRP_preEcross};
        for iEsplit_preEcross = 1:length(EsplitpreEcrossCond)
            splitE_preEcross_nm = EsplitpreEcrossCond{iEsplit_preEcross};
            %
            preEcrossModel_choicehE         = GLMprm.preEffortCross.(task_id).(RP_preEcross_nm).(splitE_preEcross_nm).choiceHighE;
            preEcrossModel_money_chosen     = GLMprm.preEffortCross.(task_id).(RP_preEcross_nm).(splitE_preEcross_nm).money_chosen;
            preEcrossModel_effort_chosen    = GLMprm.preEffortCross.(task_id).(RP_preEcross_nm).(splitE_preEcross_nm).E_chosen;
            switch task_id
                case 'Ep'
                    preEcrossModel_F_peak           = GLMprm.preEffortCross.(task_id).(RP_preEcross_nm).(splitE_preEcross_nm).F_peak;
                    preEcrossModel_F_integral       = GLMprm.preEffortCross.(task_id).(RP_preEcross_nm).(splitE_preEcross_nm).F_integral;
                    preEcrossModel_RT_avg = 0;
                    preEcrossModel_n_correct = 0;
                    preEcrossModel_n_errors = 0;
                case 'Em'
                    preEcrossModel_F_peak = 0;
                    preEcrossModel_F_integral = 0;
                    preEcrossModel_RT_avg           = GLMprm.preEffortCross.(task_id).(RP_preEcross_nm).(splitE_preEcross_nm).RT_avg;
                    preEcrossModel_n_correct        = GLMprm.preEffortCross.(task_id).(RP_preEcross_nm).(splitE_preEcross_nm).n_correct;
                    preEcrossModel_n_errors         = GLMprm.preEffortCross.(task_id).(RP_preEcross_nm).(splitE_preEcross_nm).n_errors;
            end
            preEcrossModel_NV_chosen        = GLMprm.preEffortCross.(task_id).(RP_preEcross_nm).(splitE_preEcross_nm).NV_chosen;
            preEcrossModel_NV_varOption     = GLMprm.preEffortCross.(task_id).(RP_preEcross_nm).(splitE_preEcross_nm).NV_varOption;
            preEcrossModel_NV_varOption_bis = GLMprm.preEffortCross.(task_id).(RP_preEcross_nm).(splitE_preEcross_nm).NV_varOption_bis;
            preEcrossModel_RT1stAnswer      = GLMprm.preEffortCross.(task_id).(RP_preEcross_nm).(splitE_preEcross_nm).RT_1stAnswer;
            preEcrossModel_trialN           = GLMprm.preEffortCross.(task_id).(RP_preEcross_nm).(splitE_preEcross_nm).trialN;
            preEcrossModel_conf             = GLMprm.preEffortCross.(task_id).(RP_preEcross_nm).(splitE_preEcross_nm).confidence;
            
            %% adapt depending on DCM_mode
            switch DCM_mode
                case 1
                    for iRun = 1:n_runs
                        run_full_nm = ['run',num2str(iRun)];
                        
                        % extract trial index for the current loop
                        switch RP_preEcross_nm
                            case 'RP'
                                RPfilter_preEcross = true(1,length(double(RP_var_binary.(run_full_nm))));
                            case 'R'
                                RPfilter_preEcross = (RP_var_binary.(run_full_nm) == 1);
                            case 'P'
                                RPfilter_preEcross = (RP_var_binary.(run_full_nm) == 0);
                        end
                        switch splitE_preEcross_nm
                            case 'E'
                                Efilter_preEcross = true(1,length(double(RP_var_binary.(run_full_nm))));
                            case 'E1'
                                Efilter_preEcross = (E_varOption.(run_full_nm) == 1);
                            case 'E2'
                                Efilter_preEcross = (E_varOption.(run_full_nm) == 2);
                            case 'E3'
                                Efilter_preEcross = (E_varOption.(run_full_nm) == 3);
                            case 'Ech0'
                                Efilter_preEcross = (E_chosen.(run_full_nm) == 0);
                            case 'Ech1'
                                Efilter_preEcross = (E_chosen.(run_full_nm) == 1);
                            case 'Ech2'
                                Efilter_preEcross = (E_chosen.(run_full_nm) == 2);
                            case 'Ech3'
                                Efilter_preEcross = (E_chosen.(run_full_nm) == 3);
                            case 'lEch'
                                Efilter_preEcross = (choice_hE.(run_full_nm) == 0);
                            case 'hEch'
                                Efilter_preEcross = (choice_hE.(run_full_nm) == 1);
                        end
                        preEcross_trial_idx = (RPfilter_preEcross.*Efilter_preEcross) == 1; % NEED to transform it into logical or will just focus on the first trial
                        
                        %% pre-effort cross onset
                        iCond = iCond + 1;
                        modelpreEcrossOnset = onsets.preEffortCrossOnsets.(run_full_nm)(preEcross_trial_idx);
                        % duration
                        switch preEffortCrossModel
                            case 'stick'
                                modelPreEffortCrossdur = 0;
                            case 'boxcar'
                                modelPreEffortCrossdur = durations.preEffortCrossDur.(run_full_nm)(preEcross_trial_idx);
                            case 'boxcar_bis'
                                modelPreEffortCrossdur = durations.preEffortCrossDur.(run_full_nm)(preEcross_trial_idx) +...
                                    durations.EperfDur.(run_full_nm)(preEcross_trial_idx);
                        end
                        
                        %% pre-effort cross modulators
                        n_preEcrossMods = 0;
                        preEcross_modNames = cell(1,1);
                        preEcross_modVals = [];
                        
                        
                        % choice = high effort
                        switch preEcrossModel_choicehE
                            case 0
                            case 1
                                n_preEcrossMods = n_preEcrossMods + 1;
                                preEcross_modNames{n_preEcrossMods} = 'choice = high effort';
                                preEcross_modVals(n_preEcrossMods,:) = choice_hE.(run_full_nm)(preEcross_trial_idx); % binary variable => no zscore
                            case 2
                                n_preEcrossMods = n_preEcrossMods + 1;
                                preEcross_modNames{n_preEcrossMods} = 'choice = high effort';
                                preEcross_modVals(n_preEcrossMods,:) = choice_hE_bis.(run_full_nm)(preEcross_trial_idx); % binary variable => no zscore
                            otherwise
                                error('not ready yet');
                        end
                        
                        % money chosen
                        if preEcrossModel_money_chosen > 0
                            n_preEcrossMods = n_preEcrossMods + 1;
                            preEcross_modNames{n_preEcrossMods} = 'money chosen';
                            switch preEcrossModel_money_chosen
                                case 1
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(money_amount_chosen.(run_full_nm)(preEcross_trial_idx));
                                case 2
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(abs_money_amount_chosen.(run_full_nm)(preEcross_trial_idx));
                                case 3
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(money_level_chosen.(run_full_nm)(preEcross_trial_idx));
                                case 4
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(abs_money_level_chosen.(run_full_nm)(preEcross_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % effort chosen
                        if preEcrossModel_effort_chosen > 0
                            n_preEcrossMods = n_preEcrossMods + 1;
                            preEcross_modNames{n_preEcrossMods} = 'effort chosen';
                            switch preEcrossModel_effort_chosen
                                case 1
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(E_chosen.(run_full_nm)(preEcross_trial_idx));
                                case 3
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(E_chosen_bis.(run_full_nm)(preEcross_trial_idx));
                                case 4
                                    preEcross_modVals(n_preEcrossMods,:) = zscore(E_chosen.(run_full_nm)(preEcross_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % force peak
                        if strcmp(task_id,'Ep')
                            % force peak
                            if preEcrossModel_F_peak > 0
                                n_preEcrossMods = n_preEcrossMods + 1;
                                preEcross_modNames{n_preEcrossMods} = 'force peak';
                                switch preEcrossModel_F_peak
                                    case 1
                                        preEcross_modVals(n_preEcrossMods,:) = raw_or_z(forcePeak.(run_full_nm)(preEcross_trial_idx));
                                    case 2
                                        preEcross_modVals(n_preEcrossMods,:) = raw_or_z(forcePeak_N.(run_full_nm)(preEcross_trial_idx));
                                    case 3
                                        preEcross_modVals(n_preEcrossMods,:) = zscore(forcePeak.(run_full_nm)(preEcross_trial_idx));
                                    case 4
                                        preEcross_modVals(n_preEcrossMods,:) = zscore(forcePeak_N.(run_full_nm)(preEcross_trial_idx));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                        end
                        
                        % force integral
                        if preEcrossModel_F_integral > 0
                            error('case not ready yet.');
                        end
                        
                        % RT average
                        if preEcrossModel_RT_avg > 0
                            error('case not ready yet.');
                        end
                        
                        % number of correct answers
                        if preEcrossModel_n_correct > 0
                            error('case not ready yet.');
                        end
                        
                        % number of errors
                        if preEcrossModel_n_errors > 0
                            error('case not ready yet.');
                        end
                        
                        % net value chosen
                        if preEcrossModel_NV_chosen > 0
                            n_preEcrossMods = n_preEcrossMods + 1;
                            switch preEcrossModel_NV_chosen
                                case 1
                                    preEcross_modNames{n_preEcrossMods} = 'NVch-NVunch';
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(NV_ch_min_unch.(run_full_nm)(preEcross_trial_idx));
                                case 2
                                    preEcross_modNames{n_preEcrossMods} = 'p(chosen)';
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(pChosen.(run_full_nm)(preEcross_trial_idx));
                                case 3
                                    preEcross_modNames{n_preEcrossMods} = 'NVch-NVunch';
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(NV_ch_min_unch_with_bias.(run_full_nm)(preEcross_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % net value non-default option
                        if preEcrossModel_NV_varOption > 0
                            n_preEcrossMods = n_preEcrossMods + 1;
                            switch preEcrossModel_NV_varOption
                                case 1
                                    preEcross_modNames{n_preEcrossMods} = 'delta NV high E - low E';
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(NV_varOption.(run_full_nm)(preEcross_trial_idx));
                                case 2
                                    preEcross_modNames{n_preEcrossMods} = '|delta NV high E - low E|';
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(abs(NV_varOption.(run_full_nm)(preEcross_trial_idx)));
                                case 3
                                    preEcross_modNames{n_preEcrossMods} = 'p(choice=hE)';
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(abs(pChoice_hE.(run_full_nm)(preEcross_trial_idx)));
                                case 4
                                    preEcross_modNames{n_preEcrossMods} = 'delta NV high E - low E + bias';
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(NV_varOption_plus_bias.(run_full_nm)(preEcross_trial_idx));
                                case 5
                                    preEcross_modNames{n_preEcrossMods} = '|delta NV high E - low E + bias|';
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(abs(NV_varOption_plus_bias.(run_full_nm)(preEcross_trial_idx)));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % net value non-default option bis
                        if preEcrossModel_NV_varOption_bis > 0
                            n_preEcrossMods = n_preEcrossMods + 1;
                            switch preEcrossModel_NV_varOption_bis
                                case 1
                                    preEcross_modNames{n_preEcrossMods} = 'delta NV high E - low E';
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(NV_varOption.(run_full_nm)(preEcross_trial_idx));
                                case 2
                                    preEcross_modNames{n_preEcrossMods} = '|delta NV high E - low E|';
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(abs(NV_varOption.(run_full_nm)(preEcross_trial_idx)));
                                case 3
                                    preEcross_modNames{n_preEcrossMods} = 'p(choice=hE)';
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(abs(pChoice_hE.(run_full_nm)(preEcross_trial_idx)));
                                case 4
                                    preEcross_modNames{n_preEcrossMods} = 'delta NV high E - low E + bias';
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(NV_varOption_plus_bias.(run_full_nm)(preEcross_trial_idx));
                                case 5
                                    preEcross_modNames{n_preEcrossMods} = '|delta NV high E - low E + bias|';
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(abs(NV_varOption_plus_bias.(run_full_nm)(preEcross_trial_idx)));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % RT 1st answer
                        if preEcrossModel_RT1stAnswer > 0
                            n_preEcrossMods = n_preEcrossMods + 1;
                            preEcross_modNames{n_preEcrossMods} = 'Effort latency';
                            switch preEcrossModel_RT1stAnswer
                                case 1
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(latency.(run_full_nm)(preEcross_trial_idx));
                                    %             preEcross_modVals(n_preEcrossMods,:) = ;
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % trial number
                        if preEcrossModel_trialN > 0
                            n_preEcrossMods = n_preEcrossMods + 1;
                            preEcross_modNames{n_preEcrossMods} = 'trial number';
                            switch preEcrossModel_trialN
                                case 1
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(trialN.(run_full_nm)(preEcross_trial_idx));
                                case 2
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(trialN_dEch.(run_full_nm)(preEcross_trial_idx));
                                case 3
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(trialN_dEnonDef_min_Edef.(run_full_nm)(preEcross_trial_idx));
                                case 4
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(trialN_dEnonDef.(run_full_nm)(preEcross_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % choice confidence
                        if preEcrossModel_conf > 0
                            n_preEcrossMods = n_preEcrossMods + 1;
                            preEcross_modNames{n_preEcrossMods} = 'confidence';
                            switch preEcrossModel_conf
                                case 1 % binary variable => no zscore
                                    preEcross_modVals(n_preEcrossMods,:) = confidence.(run_full_nm)(preEcross_trial_idx);
                                case {2,3,4} % confidence inferred by the model => ok to zscore
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(confidence.(run_full_nm)(preEcross_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        [matlabbatch] = First_level_loadEachCondition_DCM(matlabbatch, sub_idx, iCond,...
                            ['preEffort fixation cross',RP_preEcross_nm,'_',splitE_preEcross_nm], modelpreEcrossOnset, modelPreEffortCrossdur,...
                            n_preEcrossMods, preEcross_modNames, preEcross_modVals,...
                            orth_vars, onsets_only_GLM);
                    end % run loop
                case 2 % all tasks independent but sessions pooled
                    for iTask = 1:nTasks
                        task_nm = tasks{iTask};
                        
                        % extract trial index for the current loop
                        switch RP_preEcross_nm
                            case 'RP'
                                RPfilter_preEcross = true(1,length(double(RP_var_binary.(task_nm))));
                            case 'R'
                                RPfilter_preEcross = (RP_var_binary.(task_nm) == 1);
                            case 'P'
                                RPfilter_preEcross = (RP_var_binary.(task_nm) == 0);
                        end
                        switch splitE_preEcross_nm
                            case 'E'
                                Efilter_preEcross = true(1,length(double(RP_var_binary.(task_nm))));
                            case 'E1'
                                Efilter_preEcross = (E_varOption.(task_nm) == 1);
                            case 'E2'
                                Efilter_preEcross = (E_varOption.(task_nm) == 2);
                            case 'E3'
                                Efilter_preEcross = (E_varOption.(task_nm) == 3);
                            case 'Ech0'
                                Efilter_preEcross = (E_chosen.(task_nm) == 0);
                            case 'Ech1'
                                Efilter_preEcross = (E_chosen.(task_nm) == 1);
                            case 'Ech2'
                                Efilter_preEcross = (E_chosen.(task_nm) == 2);
                            case 'Ech3'
                                Efilter_preEcross = (E_chosen.(task_nm) == 3);
                            case 'lEch'
                                Efilter_preEcross = (choice_hE.(task_nm) == 0);
                            case 'hEch'
                                Efilter_preEcross = (choice_hE.(task_nm) == 1);
                        end
                        preEcross_trial_idx = (RPfilter_preEcross.*Efilter_preEcross) == 1; % NEED to transform it into logical or will just focus on the first trial
                        
                        
                        %% pre-effort cross onset
                        iCond = iCond + 1;
                        modelpreEcrossOnset = onsets.preEffortCrossOnsets.(task_nm)(preEcross_trial_idx);
                        % duration
                        switch preEffortCrossModel
                            case 'stick'
                                modelPreEffortCrossdur = 0;
                            case 'boxcar'
                                modelPreEffortCrossdur = durations.preEffortCrossDur.(task_nm)(preEcross_trial_idx);
                            case 'boxcar_bis'
                                modelPreEffortCrossdur = durations.preEffortCrossDur.(task_nm)(preEcross_trial_idx) +...
                                    durations.EperfDur.(task_nm)(preEcross_trial_idx);
                        end
                        
                        %% pre-effort cross modulators
                        n_preEcrossMods = 0;
                        preEcross_modNames = cell(1,1);
                        preEcross_modVals = [];
                        
                        
                        % choice = high effort
                        switch preEcrossModel_choicehE
                            case 0
                            case 1
                                n_preEcrossMods = n_preEcrossMods + 1;
                                preEcross_modNames{n_preEcrossMods} = 'choice = high effort';
                                preEcross_modVals(n_preEcrossMods,:) = choice_hE.(task_nm)(preEcross_trial_idx); % binary variable => no zscore
                            case 2
                                n_preEcrossMods = n_preEcrossMods + 1;
                                preEcross_modNames{n_preEcrossMods} = 'choice = high effort';
                                preEcross_modVals(n_preEcrossMods,:) = choice_hE_bis.(task_nm)(preEcross_trial_idx); % binary variable => no zscore
                            otherwise
                                error('not ready yet');
                        end
                        
                        % money chosen
                        if preEcrossModel_money_chosen > 0
                            n_preEcrossMods = n_preEcrossMods + 1;
                            preEcross_modNames{n_preEcrossMods} = 'money chosen';
                            switch preEcrossModel_money_chosen
                                case 1
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(money_amount_chosen.(task_nm)(preEcross_trial_idx));
                                case 2
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(abs_money_amount_chosen.(task_nm)(preEcross_trial_idx));
                                case 3
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(money_level_chosen.(task_nm)(preEcross_trial_idx));
                                case 4
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(abs_money_level_chosen.(task_nm)(preEcross_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % effort chosen
                        if preEcrossModel_effort_chosen > 0
                            n_preEcrossMods = n_preEcrossMods + 1;
                            preEcross_modNames{n_preEcrossMods} = 'effort chosen';
                            switch preEcrossModel_effort_chosen
                                case 1
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(E_chosen.(task_nm)(preEcross_trial_idx));
                                case 3
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(E_chosen_bis.(task_nm)(preEcross_trial_idx));
                                case 4
                                    preEcross_modVals(n_preEcrossMods,:) = zscore(E_chosen.(task_nm)(preEcross_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % force peak
                        if strcmp(task_id,'Ep')
                            % force peak
                            if preEcrossModel_F_peak > 0
                                n_preEcrossMods = n_preEcrossMods + 1;
                                preEcross_modNames{n_preEcrossMods} = 'force peak';
                                switch preEcrossModel_F_peak
                                    case 1
                                        preEcross_modVals(n_preEcrossMods,:) = raw_or_z(forcePeak.(task_nm)(preEcross_trial_idx));
                                    case 2
                                        preEcross_modVals(n_preEcrossMods,:) = raw_or_z(forcePeak_N.(task_nm)(preEcross_trial_idx));
                                    case 3
                                        preEcross_modVals(n_preEcrossMods,:) = zscore(forcePeak.(task_nm)(preEcross_trial_idx));
                                    case 4
                                        preEcross_modVals(n_preEcrossMods,:) = zscore(forcePeak_N.(task_nm)(preEcross_trial_idx));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                        end
                        
                        % force integral
                        if preEcrossModel_F_integral > 0
                            error('case not ready yet.');
                        end
                        
                        % RT average
                        if preEcrossModel_RT_avg > 0
                            error('case not ready yet.');
                        end
                        
                        % number of correct answers
                        if preEcrossModel_n_correct > 0
                            error('case not ready yet.');
                        end
                        
                        % number of errors
                        if preEcrossModel_n_errors > 0
                            error('case not ready yet.');
                        end
                        
                        % net value chosen
                        if preEcrossModel_NV_chosen > 0
                            n_preEcrossMods = n_preEcrossMods + 1;
                            switch preEcrossModel_NV_chosen
                                case 1
                                    preEcross_modNames{n_preEcrossMods} = 'NVch-NVunch';
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(NV_ch_min_unch.(task_nm)(preEcross_trial_idx));
                                case 2
                                    preEcross_modNames{n_preEcrossMods} = 'p(chosen)';
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(pChosen.(task_nm)(preEcross_trial_idx));
                                case 3
                                    preEcross_modNames{n_preEcrossMods} = 'NVch-NVunch';
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(NV_ch_min_unch_with_bias.(task_nm)(preEcross_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % net value non-default option
                        if preEcrossModel_NV_varOption > 0
                            n_preEcrossMods = n_preEcrossMods + 1;
                            switch preEcrossModel_NV_varOption
                                case 1
                                    preEcross_modNames{n_preEcrossMods} = 'delta NV high E - low E';
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(NV_varOption.(task_nm)(preEcross_trial_idx));
                                case 2
                                    preEcross_modNames{n_preEcrossMods} = '|delta NV high E - low E|';
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(abs(NV_varOption.(task_nm)(preEcross_trial_idx)));
                                case 3
                                    preEcross_modNames{n_preEcrossMods} = 'p(choice=hE)';
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(abs(pChoice_hE.(task_nm)(preEcross_trial_idx)));
                                case 4
                                    preEcross_modNames{n_preEcrossMods} = 'delta NV high E - low E + bias';
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(NV_varOption_plus_bias.(task_nm)(preEcross_trial_idx));
                                case 5
                                    preEcross_modNames{n_preEcrossMods} = '|delta NV high E - low E + bias|';
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(abs(NV_varOption_plus_bias.(task_nm)(preEcross_trial_idx)));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % net value non-default option bis
                        if preEcrossModel_NV_varOption_bis > 0
                            n_preEcrossMods = n_preEcrossMods + 1;
                            switch preEcrossModel_NV_varOption_bis
                                case 1
                                    preEcross_modNames{n_preEcrossMods} = 'delta NV high E - low E';
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(NV_varOption.(task_nm)(preEcross_trial_idx));
                                case 2
                                    preEcross_modNames{n_preEcrossMods} = '|delta NV high E - low E|';
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(abs(NV_varOption.(task_nm)(preEcross_trial_idx)));
                                case 3
                                    preEcross_modNames{n_preEcrossMods} = 'p(choice=hE)';
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(abs(pChoice_hE.(task_nm)(preEcross_trial_idx)));
                                case 4
                                    preEcross_modNames{n_preEcrossMods} = 'delta NV high E - low E + bias';
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(NV_varOption_plus_bias.(task_nm)(preEcross_trial_idx));
                                case 5
                                    preEcross_modNames{n_preEcrossMods} = '|delta NV high E - low E + bias|';
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(abs(NV_varOption_plus_bias.(task_nm)(preEcross_trial_idx)));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % RT 1st answer
                        if preEcrossModel_RT1stAnswer > 0
                            n_preEcrossMods = n_preEcrossMods + 1;
                            preEcross_modNames{n_preEcrossMods} = 'Effort latency';
                            switch preEcrossModel_RT1stAnswer
                                case 1
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(latency.(task_nm)(preEcross_trial_idx));
                                    %             preEcross_modVals(n_preEcrossMods,:) = ;
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % trial number
                        if preEcrossModel_trialN > 0
                            n_preEcrossMods = n_preEcrossMods + 1;
                            preEcross_modNames{n_preEcrossMods} = 'trial number';
                            switch preEcrossModel_trialN
                                case 1
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(trialN.(task_nm)(preEcross_trial_idx));
                                case 2
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(trialN_dEch.(task_nm)(preEcross_trial_idx));
                                case 3
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(trialN_dEnonDef_min_Edef.(task_nm)(preEcross_trial_idx));
                                case 4
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(trialN_dEnonDef.(task_nm)(preEcross_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % choice confidence
                        if preEcrossModel_conf > 0
                            n_preEcrossMods = n_preEcrossMods + 1;
                            preEcross_modNames{n_preEcrossMods} = 'confidence';
                            switch preEcrossModel_conf
                                case 1 % binary variable => no zscore
                                    preEcross_modVals(n_preEcrossMods,:) = confidence.(task_nm)(preEcross_trial_idx);
                                case {2,3,4} % confidence inferred by the model => ok to zscore
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(confidence.(task_nm)(preEcross_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        [matlabbatch] = First_level_loadEachCondition_DCM(matlabbatch, sub_idx, iCond,...
                            ['preEffort fixation cross',RP_preEcross_nm,'_',splitE_preEcross_nm], modelpreEcrossOnset, modelPreEffortCrossdur,...
                            n_preEcrossMods, preEcross_modNames, preEcross_modVals,...
                            orth_vars, onsets_only_GLM);
                    end % task loop
                case {3,4,5} % all fixation crosses pooled across sessions
                    % extract trial index for the current loop
                    switch RP_preEcross_nm
                        case 'RP'
                            RPfilter_preEcross = true(1,length(double(RP_var_binary.allTrials)));
                        case 'R'
                            RPfilter_preEcross = (RP_var_binary.allTrials == 1);
                        case 'P'
                            RPfilter_preEcross = (RP_var_binary.allTrials == 0);
                    end
                    switch splitE_preEcross_nm
                        case 'E'
                            Efilter_preEcross = true(1,length(double(RP_var_binary.allTrials)));
                        case 'E1'
                            Efilter_preEcross = (E_varOption.allTrials == 1);
                        case 'E2'
                            Efilter_preEcross = (E_varOption.allTrials == 2);
                        case 'E3'
                            Efilter_preEcross = (E_varOption.allTrials == 3);
                        case 'Ech0'
                            Efilter_preEcross = (E_chosen.allTrials == 0);
                        case 'Ech1'
                            Efilter_preEcross = (E_chosen.allTrials == 1);
                        case 'Ech2'
                            Efilter_preEcross = (E_chosen.allTrials == 2);
                        case 'Ech3'
                            Efilter_preEcross = (E_chosen.allTrials == 3);
                        case 'lEch'
                            Efilter_preEcross = (choice_hE.allTrials == 0);
                        case 'hEch'
                            Efilter_preEcross = (choice_hE.allTrials == 1);
                    end
                    preEcross_trial_idx = (RPfilter_preEcross.*Efilter_preEcross) == 1; % NEED to transform it into logical or will just focus on the first trial
                    
                    
                    %% pre-effort cross onset
                    iCond = iCond + 1;
                    modelpreEcrossOnset = onsets.preEffortCrossOnsets.allTrials(preEcross_trial_idx);
                    % duration
                    switch preEffortCrossModel
                        case 'stick'
                            modelPreEffortCrossdur = 0;
                        case 'boxcar'
                            modelPreEffortCrossdur = durations.preEffortCrossDur.allTrials(preEcross_trial_idx);
                        case 'boxcar_bis'
                            modelPreEffortCrossdur = durations.preEffortCrossDur.allTrials(preEcross_trial_idx) +...
                                durations.EperfDur.allTrials(preEcross_trial_idx);
                    end
                    
                    %% pre-effort cross modulators
                    n_preEcrossMods = 0;
                    preEcross_modNames = cell(1,1);
                    preEcross_modVals = [];
                    
                    
                    % choice = high effort
                    switch preEcrossModel_choicehE
                        case 0
                        case 1
                            n_preEcrossMods = n_preEcrossMods + 1;
                            preEcross_modNames{n_preEcrossMods} = 'choice = high effort';
                            preEcross_modVals(n_preEcrossMods,:) = choice_hE.allTrials(preEcross_trial_idx); % binary variable => no zscore
                        case 2
                            n_preEcrossMods = n_preEcrossMods + 1;
                            preEcross_modNames{n_preEcrossMods} = 'choice = high effort';
                            preEcross_modVals(n_preEcrossMods,:) = choice_hE_bis.allTrials(preEcross_trial_idx); % binary variable => no zscore
                        otherwise
                            error('not ready yet');
                    end
                    
                    % money chosen
                    if preEcrossModel_money_chosen > 0
                        n_preEcrossMods = n_preEcrossMods + 1;
                        preEcross_modNames{n_preEcrossMods} = 'money chosen';
                        switch preEcrossModel_money_chosen
                            case 1
                                preEcross_modVals(n_preEcrossMods,:) = raw_or_z(money_amount_chosen.allTrials(preEcross_trial_idx));
                            case 2
                                preEcross_modVals(n_preEcrossMods,:) = raw_or_z(abs_money_amount_chosen.allTrials(preEcross_trial_idx));
                            case 3
                                preEcross_modVals(n_preEcrossMods,:) = raw_or_z(money_level_chosen.allTrials(preEcross_trial_idx));
                            case 4
                                preEcross_modVals(n_preEcrossMods,:) = raw_or_z(abs_money_level_chosen.allTrials(preEcross_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % effort chosen
                    if preEcrossModel_effort_chosen > 0
                        n_preEcrossMods = n_preEcrossMods + 1;
                        preEcross_modNames{n_preEcrossMods} = 'effort chosen';
                        switch preEcrossModel_effort_chosen
                            case 1
                                preEcross_modVals(n_preEcrossMods,:) = raw_or_z(E_chosen.allTrials(preEcross_trial_idx));
                            case 3
                                preEcross_modVals(n_preEcrossMods,:) = raw_or_z(E_chosen_bis.allTrials(preEcross_trial_idx));
                            case 4
                                preEcross_modVals(n_preEcrossMods,:) = zscore(E_chosen.allTrials(preEcross_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % force peak
                    if strcmp(task_id,'Ep')
                        % force peak
                        if preEcrossModel_F_peak > 0
                            n_preEcrossMods = n_preEcrossMods + 1;
                            preEcross_modNames{n_preEcrossMods} = 'force peak';
                            switch preEcrossModel_F_peak
                                case 1
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(forcePeak.allTrials(preEcross_trial_idx));
                                case 2
                                    preEcross_modVals(n_preEcrossMods,:) = raw_or_z(forcePeak_N.allTrials(preEcross_trial_idx));
                                case 3
                                    preEcross_modVals(n_preEcrossMods,:) = zscore(forcePeak.allTrials(preEcross_trial_idx));
                                case 4
                                    preEcross_modVals(n_preEcrossMods,:) = zscore(forcePeak_N.allTrials(preEcross_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                    end
                    
                    % force integral
                    if preEcrossModel_F_integral > 0
                        error('case not ready yet.');
                    end
                    
                    % RT average
                    if preEcrossModel_RT_avg > 0
                        error('case not ready yet.');
                    end
                    
                    % number of correct answers
                    if preEcrossModel_n_correct > 0
                        error('case not ready yet.');
                    end
                    
                    % number of errors
                    if preEcrossModel_n_errors > 0
                        error('case not ready yet.');
                    end
                    
                    % net value chosen
                    if preEcrossModel_NV_chosen > 0
                        n_preEcrossMods = n_preEcrossMods + 1;
                        switch preEcrossModel_NV_chosen
                            case 1
                                preEcross_modNames{n_preEcrossMods} = 'NVch-NVunch';
                                preEcross_modVals(n_preEcrossMods,:) = raw_or_z(NV_ch_min_unch.allTrials(preEcross_trial_idx));
                            case 2
                                preEcross_modNames{n_preEcrossMods} = 'p(chosen)';
                                preEcross_modVals(n_preEcrossMods,:) = raw_or_z(pChosen.allTrials(preEcross_trial_idx));
                            case 3
                                preEcross_modNames{n_preEcrossMods} = 'NVch-NVunch';
                                preEcross_modVals(n_preEcrossMods,:) = raw_or_z(NV_ch_min_unch_with_bias.allTrials(preEcross_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % net value non-default option
                    if preEcrossModel_NV_varOption > 0
                        n_preEcrossMods = n_preEcrossMods + 1;
                        switch preEcrossModel_NV_varOption
                            case 1
                                preEcross_modNames{n_preEcrossMods} = 'delta NV high E - low E';
                                preEcross_modVals(n_preEcrossMods,:) = raw_or_z(NV_varOption.allTrials(preEcross_trial_idx));
                            case 2
                                preEcross_modNames{n_preEcrossMods} = '|delta NV high E - low E|';
                                preEcross_modVals(n_preEcrossMods,:) = raw_or_z(abs(NV_varOption.allTrials(preEcross_trial_idx)));
                            case 3
                                preEcross_modNames{n_preEcrossMods} = 'p(choice=hE)';
                                preEcross_modVals(n_preEcrossMods,:) = raw_or_z(abs(pChoice_hE.allTrials(preEcross_trial_idx)));
                            case 4
                                preEcross_modNames{n_preEcrossMods} = 'delta NV high E - low E + bias';
                                preEcross_modVals(n_preEcrossMods,:) = raw_or_z(NV_varOption_plus_bias.allTrials(preEcross_trial_idx));
                            case 5
                                preEcross_modNames{n_preEcrossMods} = '|delta NV high E - low E + bias|';
                                preEcross_modVals(n_preEcrossMods,:) = raw_or_z(abs(NV_varOption_plus_bias.allTrials(preEcross_trial_idx)));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % net value non-default option bis
                    if preEcrossModel_NV_varOption_bis > 0
                        n_preEcrossMods = n_preEcrossMods + 1;
                        switch preEcrossModel_NV_varOption_bis
                            case 1
                                preEcross_modNames{n_preEcrossMods} = 'delta NV high E - low E';
                                preEcross_modVals(n_preEcrossMods,:) = raw_or_z(NV_varOption.allTrials(preEcross_trial_idx));
                            case 2
                                preEcross_modNames{n_preEcrossMods} = '|delta NV high E - low E|';
                                preEcross_modVals(n_preEcrossMods,:) = raw_or_z(abs(NV_varOption.allTrials(preEcross_trial_idx)));
                            case 3
                                preEcross_modNames{n_preEcrossMods} = 'p(choice=hE)';
                                preEcross_modVals(n_preEcrossMods,:) = raw_or_z(abs(pChoice_hE.allTrials(preEcross_trial_idx)));
                            case 4
                                preEcross_modNames{n_preEcrossMods} = 'delta NV high E - low E + bias';
                                preEcross_modVals(n_preEcrossMods,:) = raw_or_z(NV_varOption_plus_bias.allTrials(preEcross_trial_idx));
                            case 5
                                preEcross_modNames{n_preEcrossMods} = '|delta NV high E - low E + bias|';
                                preEcross_modVals(n_preEcrossMods,:) = raw_or_z(abs(NV_varOption_plus_bias.allTrials(preEcross_trial_idx)));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % RT 1st answer
                    if preEcrossModel_RT1stAnswer > 0
                        n_preEcrossMods = n_preEcrossMods + 1;
                        preEcross_modNames{n_preEcrossMods} = 'Effort latency';
                        switch preEcrossModel_RT1stAnswer
                            case 1
                                preEcross_modVals(n_preEcrossMods,:) = raw_or_z(latency.allTrials(preEcross_trial_idx));
                                %             preEcross_modVals(n_preEcrossMods,:) = ;
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % trial number
                    if preEcrossModel_trialN > 0
                        n_preEcrossMods = n_preEcrossMods + 1;
                        preEcross_modNames{n_preEcrossMods} = 'trial number';
                        switch preEcrossModel_trialN
                            case 1
                                preEcross_modVals(n_preEcrossMods,:) = raw_or_z(trialN.allTrials(preEcross_trial_idx));
                            case 2
                                preEcross_modVals(n_preEcrossMods,:) = raw_or_z(trialN_dEch.allTrials(preEcross_trial_idx));
                            case 3
                                preEcross_modVals(n_preEcrossMods,:) = raw_or_z(trialN_dEnonDef_min_Edef.allTrials(preEcross_trial_idx));
                            case 4
                                preEcross_modVals(n_preEcrossMods,:) = raw_or_z(trialN_dEnonDef.allTrials(preEcross_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % choice confidence
                    if preEcrossModel_conf > 0
                        n_preEcrossMods = n_preEcrossMods + 1;
                        preEcross_modNames{n_preEcrossMods} = 'confidence';
                        switch preEcrossModel_conf
                            case 1 % binary variable => no zscore
                                preEcross_modVals(n_preEcrossMods,:) = confidence.allTrials(preEcross_trial_idx);
                            case {2,3,4} % confidence inferred by the model => ok to zscore
                                preEcross_modVals(n_preEcrossMods,:) = raw_or_z(confidence.allTrials(preEcross_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    [matlabbatch] = First_level_loadEachCondition_DCM(matlabbatch, sub_idx, iCond,...
                        ['preEffort fixation cross',RP_preEcross_nm,'_',splitE_preEcross_nm], modelpreEcrossOnset, modelPreEffortCrossdur,...
                        n_preEcrossMods, preEcross_modNames, preEcross_modVals,...
                        orth_vars, onsets_only_GLM);
            end % DCM mode: info about how regressors concatenated or not across tasks or sessions
        end % E condition
    end % RP
end % model pre-effort cross

%% effort performance
EperfModel = GLMprm.model_onset.(task_id).Eperf;
if ismember(EperfModel,{'stick','boxcar'})
    
    for iRP_Eperf = 1:length(RPperfCond)
        RP_Eperf_nm = RPperfCond{iRP_Eperf};
        for iEsplit_Eperf = 1:length(EsplitEperfCond)
            splitE_Eperf_nm = EsplitEperfCond{iEsplit_Eperf};
            %
            EperfModel_choicehE         = GLMprm.Eperf.(task_id).(RP_Eperf_nm).(splitE_Eperf_nm).choiceHighE;
            EperfModel_money_chosen     = GLMprm.Eperf.(task_id).(RP_Eperf_nm).(splitE_Eperf_nm).money_chosen;
            EperfModel_effort_chosen    = GLMprm.Eperf.(task_id).(RP_Eperf_nm).(splitE_Eperf_nm).E_chosen;
            switch task_id
                case 'Ep'
                    EperfModel_F_peak           = GLMprm.Eperf.(task_id).(RP_Eperf_nm).(splitE_Eperf_nm).F_peak;
                    EperfModel_F_integral       = GLMprm.Eperf.(task_id).(RP_Eperf_nm).(splitE_Eperf_nm).F_integral;
                    EperfModel_fatigue          = GLMprm.Eperf.(task_id).(RP_Eperf_nm).(splitE_Eperf_nm).fatigue;
                    EperfModel_Ech_x_fatigue    = GLMprm.Eperf.(task_id).(RP_Eperf_nm).(splitE_Eperf_nm).Ech_x_fatigue;
                case 'Em'
                    EperfModel_efficacy             = GLMprm.Eperf.(task_id).(RP_Eperf_nm).(splitE_Eperf_nm).efficacy;
                    EperfModel_prevEfficacy         = GLMprm.Eperf.(task_id).(RP_Eperf_nm).(splitE_Eperf_nm).prevEfficacy;
                    EperfModel_Ech_x_prevEfficacy   = GLMprm.Eperf.(task_id).(RP_Eperf_nm).(splitE_Eperf_nm).Ech_x_prevEfficacy;
                    EperfModel_RT_avg               = GLMprm.Eperf.(task_id).(RP_Eperf_nm).(splitE_Eperf_nm).RT_avg;
                    EperfModel_n_correct            = GLMprm.Eperf.(task_id).(RP_Eperf_nm).(splitE_Eperf_nm).n_correct;
                    EperfModel_n_errors             = GLMprm.Eperf.(task_id).(RP_Eperf_nm).(splitE_Eperf_nm).n_errors;
            end
            EperfModel_NV_chosen        = GLMprm.Eperf.(task_id).(RP_Eperf_nm).(splitE_Eperf_nm).NV_chosen;
            EperfModel_NV_varOption     = GLMprm.Eperf.(task_id).(RP_Eperf_nm).(splitE_Eperf_nm).NV_varOption;
            EperfModel_NV_varOption_bis = GLMprm.Eperf.(task_id).(RP_Eperf_nm).(splitE_Eperf_nm).NV_varOption_bis;
            EperfModel_RT1stAnswer      = GLMprm.Eperf.(task_id).(RP_Eperf_nm).(splitE_Eperf_nm).RT_1stAnswer;
            EperfModel_trialN           = GLMprm.Eperf.(task_id).(RP_Eperf_nm).(splitE_Eperf_nm).trialN;
            EperfModel_conf             = GLMprm.Eperf.(task_id).(RP_Eperf_nm).(splitE_Eperf_nm).confidence;
            
            %% adapt depending on how sessions are pooled
            switch DCM_mode
                case 1
                    for iRun = 1:n_runs
                        run_full_nm = ['run',num2str(iRun)];
                        
                        % extract trial index for the current loop
                        switch RP_Eperf_nm
                            case 'RP'
                                RPfilter_Eperf = true(1,length(double(RP_var_binary.(run_full_nm))));
                            case 'R'
                                RPfilter_Eperf = (RP_var_binary.(run_full_nm) == 1);
                            case 'P'
                                RPfilter_Eperf = (RP_var_binary.(run_full_nm) == 0);
                        end
                        switch splitE_Eperf_nm
                            case 'E'
                                Efilter_Eperf = true(1,length(double(RP_var_binary.(run_full_nm))));
                            case 'E1'
                                Efilter_Eperf = (E_varOption.(run_full_nm) == 1);
                            case 'E2'
                                Efilter_Eperf = (E_varOption.(run_full_nm) == 2);
                            case 'E3'
                                Efilter_Eperf = (E_varOption.(run_full_nm) == 3);
                            case 'Ech0'
                                Efilter_Eperf = (E_chosen.(run_full_nm) == 0);
                            case 'Ech1'
                                Efilter_Eperf = (E_chosen.(run_full_nm) == 1);
                            case 'Ech2'
                                Efilter_Eperf = (E_chosen.(run_full_nm) == 2);
                            case 'Ech3'
                                Efilter_Eperf = (E_chosen.(run_full_nm) == 3);
                            case 'lEch'
                                Efilter_Eperf = (choice_hE.(run_full_nm) == 0);
                            case 'hEch'
                                Efilter_Eperf = (choice_hE.(run_full_nm) == 1);
                        end
                        switch onsets_only_GLM
                            case 0
                                Eperf_trial_idx = (RPfilter_Eperf.*Efilter_Eperf.*perfOkTrials) == 1; % NEED to transform it into logical or will just focus on the first trial
                            case 1 % keep same number of trials as in choices and feedback to avoid problems in the analysis
                                Eperf_trial_idx = (RPfilter_Eperf.*Efilter_Eperf) == 1;
                        end
                        
                        %% Effort performance onset
                        iCond = iCond + 1;
                        modelEperfOnset = onsets.EperfOnsets.(run_full_nm)(Eperf_trial_idx);
                        % duration
                        switch EperfModel
                            case 'stick'
                                modelEperfDur = 0;
                            case 'boxcar'
                                modelEperfDur = durations.EperfDur.(run_full_nm)(Eperf_trial_idx);
                        end
                        
                        %% Effort performance modulators
                        n_EperfMods = 0;
                        Eperf_modNames = cell(1,1);
                        Eperf_modVals = [];
                        
                        % choice = high effort
                        switch EperfModel_choicehE
                            case 0
                            case 1
                                n_EperfMods = n_EperfMods + 1;
                                Eperf_modNames{n_EperfMods} = 'choice = high effort';
                                Eperf_modVals(n_EperfMods,:) = choice_hE.(run_full_nm)(Eperf_trial_idx); % binary variable => no zscore
                            case 2
                                n_EperfMods = n_EperfMods + 1;
                                Eperf_modNames{n_EperfMods} = 'choice = high effort';
                                Eperf_modVals(n_EperfMods,:) = choice_hE_bis.(run_full_nm)(Eperf_trial_idx); % binary variable => no zscore
                            otherwise
                                error('not ready yet');
                        end
                        
                        % money chosen
                        if EperfModel_money_chosen > 0
                            n_EperfMods = n_EperfMods + 1;
                            Eperf_modNames{n_EperfMods} = 'money chosen';
                            switch EperfModel_money_chosen
                                case 1
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(money_amount_chosen.(run_full_nm)(Eperf_trial_idx));
                                case 2
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(abs_money_amount_chosen.(run_full_nm)(Eperf_trial_idx));
                                case 3
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(money_level_chosen.(run_full_nm)(Eperf_trial_idx));
                                case 4
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(abs_money_level_chosen.(run_full_nm)(Eperf_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % effort chosen
                        if EperfModel_effort_chosen > 0
                            n_EperfMods = n_EperfMods + 1;
                            Eperf_modNames{n_EperfMods} = 'effort chosen';
                            switch EperfModel_effort_chosen
                                case 1
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(E_chosen.(run_full_nm)(Eperf_trial_idx));
                                case 3
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(E_chosen_bis.(run_full_nm)(Eperf_trial_idx));
                                case 4
                                    Eperf_modVals(n_EperfMods,:) = zscore(E_chosen.(run_full_nm)(Eperf_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        if strcmp(task_id,'Ep')
                            % force peak
                            if EperfModel_F_peak> 0
                                n_EperfMods = n_EperfMods + 1;
                                Eperf_modNames{n_EperfMods} = 'force peak';
                                switch EperfModel_F_peak
                                    case 1
                                        Eperf_modVals(n_EperfMods,:) = raw_or_z(forcePeak.(run_full_nm)(Eperf_trial_idx));
                                    case 2
                                        Eperf_modVals(n_EperfMods,:) = raw_or_z(forcePeak_N.(run_full_nm)(Eperf_trial_idx));
                                    case 3
                                        Eperf_modVals(n_EperfMods,:) = zscore(forcePeak.(run_full_nm)(Eperf_trial_idx));
                                    case 4
                                        Eperf_modVals(n_EperfMods,:) = zscore(forcePeak_N.(run_full_nm)(Eperf_trial_idx));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                            
                            % force integral
                            if EperfModel_F_integral > 0
                                n_EperfMods = n_EperfMods + 1;
                                Eperf_modNames{n_EperfMods} = 'force integral';
                                switch EperfModel_F_integral
                                    case 1
                                        Eperf_modVals(n_EperfMods,:) = raw_or_z(AUC.(run_full_nm)(Eperf_trial_idx));
                                    case 2
                                        Eperf_modVals(n_EperfMods,:) = raw_or_z(AUC_overshoot.(run_full_nm)(Eperf_trial_idx));
                                    case 3
                                        Eperf_modVals(n_EperfMods,:) = raw_or_z(AUC_N.(run_full_nm)(Eperf_trial_idx));
                                    case 4
                                        Eperf_modVals(n_EperfMods,:) = raw_or_z(AUC_overshoot_N.(run_full_nm)(Eperf_trial_idx));
                                    case 5
                                        Eperf_modVals(n_EperfMods,:) = zscore(AUC.(run_full_nm)(Eperf_trial_idx));
                                    case 6
                                        Eperf_modVals(n_EperfMods,:) = zscore(AUC_overshoot.(run_full_nm)(Eperf_trial_idx));
                                    case 7
                                        Eperf_modVals(n_EperfMods,:) = zscore(AUC_N.(run_full_nm)(Eperf_trial_idx));
                                    case 8
                                        Eperf_modVals(n_EperfMods,:) = zscore(AUC_overshoot_N.(run_full_nm)(Eperf_trial_idx));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                        end % physical effort filter
                        
                        if strcmp(task_id,'Em')
                            % efficacy
                            if EperfModel_efficacy > 0
                                n_EperfMods = n_EperfMods + 1;
                                Eperf_modNames{n_EperfMods} = 'Em efficacy';
                                switch EperfModel_efficacy
                                    case 1
                                        Eperf_modVals(n_EperfMods,:) = raw_or_z(efficacy_with2first.(run_full_nm)(Eperf_trial_idx));
                                    case 2
                                        Eperf_modVals(n_EperfMods,:) = raw_or_z(efficacy_pureNback.(run_full_nm)(Eperf_trial_idx));
                                    case 3
                                        Eperf_modVals(n_EperfMods,:) = raw_or_z(efficacy_bis_with2first.(run_full_nm)(Eperf_trial_idx));
                                    case 4
                                        Eperf_modVals(n_EperfMods,:) = raw_or_z(efficacy_bis_pureNback.(run_full_nm)(Eperf_trial_idx));
                                    case 5
                                        Eperf_modVals(n_EperfMods,:) = zscore(efficacy_with2first.(run_full_nm)(Eperf_trial_idx));
                                    case 6
                                        Eperf_modVals(n_EperfMods,:) = zscore(efficacy_pureNback.(run_full_nm)(Eperf_trial_idx));
                                    case 7
                                        Eperf_modVals(n_EperfMods,:) = zscore(efficacy_bis_with2first.(run_full_nm)(Eperf_trial_idx));
                                    case 8
                                        Eperf_modVals(n_EperfMods,:) = zscore(efficacy_bis_pureNback.(run_full_nm)(Eperf_trial_idx));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                            
                            % RT average
                            if EperfModel_RT_avg > 0
                                n_EperfMods = n_EperfMods + 1;
                                Eperf_modNames{n_EperfMods} = 'avg RT N-back perf';
                                switch EperfModel_RT_avg
                                    case 1
                                        Eperf_modVals(n_EperfMods,:) = raw_or_z(RT_avg.(run_full_nm)(Eperf_trial_idx));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                            
                            % number of correct answers
                            if EperfModel_n_correct > 0
                                n_EperfMods = n_EperfMods + 1;
                                Eperf_modNames{n_EperfMods} = 'nb correct';
                                switch EperfModel_n_correct
                                    case 1
                                        Eperf_modVals(n_EperfMods,:) = raw_or_z(n_correct.(run_full_nm)(Eperf_trial_idx));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                            
                            % number of errors
                            if EperfModel_n_errors > 0
                                n_EperfMods = n_EperfMods + 1;
                                Eperf_modNames{n_EperfMods} = 'nb errors';
                                switch EperfModel_n_errors
                                    case 1
                                        Eperf_modVals(n_EperfMods,:) = raw_or_z(n_errors.(run_full_nm)(Eperf_trial_idx));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                        end % mental effort filter
                        
                        % net value chosen
                        if EperfModel_NV_chosen > 0
                            n_EperfMods = n_EperfMods + 1;
                            switch EperfModel_NV_chosen
                                case 1
                                    Eperf_modNames{n_EperfMods} = 'NV chosen';
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(NV_ch_min_unch.(run_full_nm)(Eperf_trial_idx));
                                case 2
                                    Eperf_modNames{n_EperfMods} = 'p(chosen)';
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(pChosen.(run_full_nm)(Eperf_trial_idx));
                                case 3
                                    Eperf_modNames{n_EperfMods} = 'NVch-NVunch';
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(NV_ch_min_unch_with_bias.(run_full_nm)(Eperf_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % net value non-default option
                        if EperfModel_NV_varOption > 0
                            n_EperfMods = n_EperfMods + 1;
                            switch EperfModel_NV_varOption
                                case 1
                                    Eperf_modNames{n_EperfMods} = 'delta NV high E - low E';
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(NV_varOption.(run_full_nm)(Eperf_trial_idx));
                                case 2
                                    Eperf_modNames{n_EperfMods} = '|delta NV high E - low E|';
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(abs(NV_varOption.(run_full_nm)(Eperf_trial_idx)));
                                case 3
                                    Eperf_modNames{n_EperfMods} = 'p(choice=hE)';
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(abs(pChoice_hE.(run_full_nm)(Eperf_trial_idx)));
                                case 4
                                    Eperf_modNames{n_EperfMods} = 'delta NV high E - low E + bias';
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(NV_varOption_plus_bias.(run_full_nm)(Eperf_trial_idx));
                                case 5
                                    Eperf_modNames{n_EperfMods} = '|delta NV high E - low E + bias|';
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(abs(NV_varOption_plus_bias.(run_full_nm)(Eperf_trial_idx)));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % net value non-default option bis
                        if EperfModel_NV_varOption_bis > 0
                            n_EperfMods = n_EperfMods + 1;
                            switch EperfModel_NV_varOption_bis
                                case 1
                                    Eperf_modNames{n_EperfMods} = 'delta NV high E - low E';
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(NV_varOption.(run_full_nm)(Eperf_trial_idx));
                                case 2
                                    Eperf_modNames{n_EperfMods} = '|delta NV high E - low E|';
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(abs(NV_varOption.(run_full_nm)(Eperf_trial_idx)));
                                case 3
                                    Eperf_modNames{n_EperfMods} = 'p(choice=hE)';
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(abs(pChoice_hE.(run_full_nm)(Eperf_trial_idx)));
                                case 4
                                    Eperf_modNames{n_EperfMods} = 'delta NV high E - low E + bias';
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(NV_varOption_plus_bias.(run_full_nm)(Eperf_trial_idx));
                                case 5
                                    Eperf_modNames{n_EperfMods} = '|delta NV high E - low E + bias|';
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(abs(NV_varOption_plus_bias.(run_full_nm)(Eperf_trial_idx)));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % RT 1st answer
                        if EperfModel_RT1stAnswer > 0
                            n_EperfMods = n_EperfMods + 1;
                            Eperf_modNames{n_EperfMods} = 'Effort latency';
                            switch EperfModel_RT1stAnswer
                                case 1
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(latency.(run_full_nm)(Eperf_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        switch task_id
                            case 'Ep'
                                % physical fatigue
                                if EperfModel_fatigue > 0
                                    n_EperfMods = n_EperfMods + 1;
                                    Eperf_modNames{n_EperfMods} = 'fatigue';
                                    switch EperfModel_fatigue
                                        case 1
                                            Eperf_modVals(n_EperfMods,:) = raw_or_z(fatigue.(run_full_nm)(Eperf_trial_idx));
                                        case 2
                                            Eperf_modVals(n_EperfMods,:) = zscore(fatigue.(run_full_nm)(Eperf_trial_idx));
                                        otherwise
                                            error('not ready yet');
                                    end
                                end
                                
                                % Ech*(physical fatigue)
                                if EperfModel_Ech_x_fatigue > 0
                                    n_EperfMods = n_EperfMods + 1;
                                    Eperf_modNames{n_EperfMods} = 'Ech_x_fatigue';
                                    switch EperfModel_Ech_x_fatigue
                                        case 1
                                            Eperf_modVals(n_EperfMods,:) = raw_or_z(E_chosen.(run_full_nm)(Eperf_trial_idx).*fatigue.(run_full_nm)(Eperf_trial_idx));
                                        otherwise
                                            error('not ready yet');
                                    end
                                end
                            case 'Em'
                                % mental facilitation based on previous trial
                                % performance
                                if EperfModel_prevEfficacy > 0
                                    n_EperfMods = n_EperfMods + 1;
                                    Eperf_modNames{n_EperfMods} = 'previous trial efficacy';
                                    switch EperfModel_prevEfficacy
                                        case 1
                                            Eperf_modVals(n_EperfMods,:) = raw_or_z(prevEfficacy_with2first.(run_full_nm)(Eperf_trial_idx));
                                        case 2
                                            Eperf_modVals(n_EperfMods,:) = raw_or_z(prevEfficacy_pureNback.(run_full_nm)(Eperf_trial_idx));
                                        case 3
                                            Eperf_modVals(n_EperfMods,:) = raw_or_z(prevEfficacy_bis_with2first.(run_full_nm)(Eperf_trial_idx));
                                        case 4
                                            Eperf_modVals(n_EperfMods,:) = raw_or_z(prevEfficacy_bis_pureNback.(run_full_nm)(Eperf_trial_idx));
                                        case 5
                                            Eperf_modVals(n_EperfMods,:) = zscore(prevEfficacy_with2first.(run_full_nm)(Eperf_trial_idx));
                                        case 6
                                            Eperf_modVals(n_EperfMods,:) = zscore(prevEfficacy_pureNback.(run_full_nm)(Eperf_trial_idx));
                                        case 7
                                            Eperf_modVals(n_EperfMods,:) = zscore(prevEfficacy_bis_with2first.(run_full_nm)(Eperf_trial_idx));
                                        case 8
                                            Eperf_modVals(n_EperfMods,:) = zscore(prevEfficacy_bis_pureNback.(run_full_nm)(Eperf_trial_idx));
                                        otherwise
                                            error('not ready yet');
                                    end
                                end
                                
                                % Echosen*(previous trial efficacy)
                                if EperfModel_Ech_x_prevEfficacy > 0
                                    n_EperfMods = n_EperfMods + 1;
                                    Eperf_modNames{n_EperfMods} = 'Ech_x_previous trial efficacy';
                                    switch EperfModel_Ech_x_prevEfficacy
                                        case 1
                                            Eperf_modVals(n_EperfMods,:) = raw_or_z(E_chosen.(run_full_nm)(Eperf_trial_idx).*prevEfficacy_with2first.(run_full_nm)(Eperf_trial_idx));
                                        case 2
                                            Eperf_modVals(n_EperfMods,:) = raw_or_z(E_chosen.(run_full_nm)(Eperf_trial_idx).*prevEfficacy_pureNback.(run_full_nm)(Eperf_trial_idx));
                                        case 3
                                            Eperf_modVals(n_EperfMods,:) = raw_or_z(E_chosen.(run_full_nm)(Eperf_trial_idx).*prevEfficacy_bis_with2first.(run_full_nm)(Eperf_trial_idx));
                                        case 4
                                            Eperf_modVals(n_EperfMods,:) = raw_or_z(E_chosen.(run_full_nm)(Eperf_trial_idx).*prevEfficacy_bis_pureNback.(run_full_nm)(Eperf_trial_idx));
                                        otherwise
                                            error('not ready yet');
                                    end
                                end
                        end % physical/mental effort filter
                        
                        % trial number
                        if EperfModel_trialN > 0
                            n_EperfMods = n_EperfMods + 1;
                            Eperf_modNames{n_EperfMods} = 'trial number';
                            switch EperfModel_trialN
                                case 1
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(trialN.(run_full_nm)(Eperf_trial_idx));
                                case 2
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(trialN_dEch.(run_full_nm)(Eperf_trial_idx));
                                case 3
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(trialN_dEnonDef_min_Edef.(run_full_nm)(Eperf_trial_idx));
                                case 4
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(trialN_dEnonDef.(run_full_nm)(Eperf_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % choice confidence
                        if EperfModel_conf > 0
                            n_EperfMods = n_EperfMods + 1;
                            Eperf_modNames{n_EperfMods} = 'confidence';
                            switch EperfModel_conf
                                case 1 % binary variable => no zscore
                                    Eperf_modVals(n_EperfMods,:) = confidence.(run_full_nm)(Eperf_trial_idx);
                                case {2,3,4} % confidence inferred by the model => ok to zscore
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(confidence.(run_full_nm)(Eperf_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        [matlabbatch] = First_level_loadEachCondition_DCM(matlabbatch, sub_idx, iCond,...
                            ['Eperf_',RP_Eperf_nm,'_',splitE_Eperf_nm], modelEperfOnset, modelEperfDur,...
                            n_EperfMods, Eperf_modNames, Eperf_modVals,...
                            orth_vars, onsets_only_GLM);
                    end % run loop
                case {2,4,5} % all tasks independent but sessions pooled
                    for iTask = 1:nTasks
                        task_nm = tasks{iTask};
                        
                        % extract trial index for the current loop
                        switch RP_Eperf_nm
                            case 'RP'
                                RPfilter_Eperf = true(1,length(double(RP_var_binary.(task_nm))));
                            case 'R'
                                RPfilter_Eperf = (RP_var_binary.(task_nm) == 1);
                            case 'P'
                                RPfilter_Eperf = (RP_var_binary.(task_nm) == 0);
                        end
                        switch splitE_Eperf_nm
                            case 'E'
                                Efilter_Eperf = true(1,length(double(RP_var_binary.(task_nm))));
                            case 'E1'
                                Efilter_Eperf = (E_varOption.(task_nm) == 1);
                            case 'E2'
                                Efilter_Eperf = (E_varOption.(task_nm) == 2);
                            case 'E3'
                                Efilter_Eperf = (E_varOption.(task_nm) == 3);
                            case 'Ech0'
                                Efilter_Eperf = (E_chosen.(task_nm) == 0);
                            case 'Ech1'
                                Efilter_Eperf = (E_chosen.(task_nm) == 1);
                            case 'Ech2'
                                Efilter_Eperf = (E_chosen.(task_nm) == 2);
                            case 'Ech3'
                                Efilter_Eperf = (E_chosen.(task_nm) == 3);
                            case 'lEch'
                                Efilter_Eperf = (choice_hE.(task_nm) == 0);
                            case 'hEch'
                                Efilter_Eperf = (choice_hE.(task_nm) == 1);
                        end
                        switch onsets_only_GLM
                            case 0
                                Eperf_trial_idx = (RPfilter_Eperf.*Efilter_Eperf.*perfOkTrials) == 1; % NEED to transform it into logical or will just focus on the first trial
                            case 1 % keep same number of trials as in choices and feedback to avoid problems in the analysis
                                Eperf_trial_idx = (RPfilter_Eperf.*Efilter_Eperf) == 1;
                        end
                        
                        %% Effort performance onset
                        iCond = iCond + 1;
                        modelEperfOnset = onsets.EperfOnsets.(run_full_nm)(Eperf_trial_idx);
                        % duration
                        switch EperfModel
                            case 'stick'
                                modelEperfDur = 0;
                            case 'boxcar'
                                modelEperfDur = durations.EperfDur.(run_full_nm)(Eperf_trial_idx);
                        end
                        
                        %% Effort performance modulators
                        n_EperfMods = 0;
                        Eperf_modNames = cell(1,1);
                        Eperf_modVals = [];
                        
                        % choice = high effort
                        switch EperfModel_choicehE
                            case 0
                            case 1
                                n_EperfMods = n_EperfMods + 1;
                                Eperf_modNames{n_EperfMods} = 'choice = high effort';
                                Eperf_modVals(n_EperfMods,:) = choice_hE.(run_full_nm)(Eperf_trial_idx); % binary variable => no zscore
                            case 2
                                n_EperfMods = n_EperfMods + 1;
                                Eperf_modNames{n_EperfMods} = 'choice = high effort';
                                Eperf_modVals(n_EperfMods,:) = choice_hE_bis.(run_full_nm)(Eperf_trial_idx); % binary variable => no zscore
                            otherwise
                                error('not ready yet');
                        end
                        
                        % money chosen
                        if EperfModel_money_chosen > 0
                            n_EperfMods = n_EperfMods + 1;
                            Eperf_modNames{n_EperfMods} = 'money chosen';
                            switch EperfModel_money_chosen
                                case 1
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(money_amount_chosen.(run_full_nm)(Eperf_trial_idx));
                                case 2
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(abs_money_amount_chosen.(run_full_nm)(Eperf_trial_idx));
                                case 3
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(money_level_chosen.(run_full_nm)(Eperf_trial_idx));
                                case 4
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(abs_money_level_chosen.(run_full_nm)(Eperf_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % effort chosen
                        if EperfModel_effort_chosen > 0
                            n_EperfMods = n_EperfMods + 1;
                            Eperf_modNames{n_EperfMods} = 'effort chosen';
                            switch EperfModel_effort_chosen
                                case 1
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(E_chosen.(run_full_nm)(Eperf_trial_idx));
                                case 3
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(E_chosen_bis.(run_full_nm)(Eperf_trial_idx));
                                case 4
                                    Eperf_modVals(n_EperfMods,:) = zscore(E_chosen.(run_full_nm)(Eperf_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        if strcmp(task_id,'Ep')
                            % force peak
                            if EperfModel_F_peak> 0
                                n_EperfMods = n_EperfMods + 1;
                                Eperf_modNames{n_EperfMods} = 'force peak';
                                switch EperfModel_F_peak
                                    case 1
                                        Eperf_modVals(n_EperfMods,:) = raw_or_z(forcePeak.(run_full_nm)(Eperf_trial_idx));
                                    case 2
                                        Eperf_modVals(n_EperfMods,:) = raw_or_z(forcePeak_N.(run_full_nm)(Eperf_trial_idx));
                                    case 3
                                        Eperf_modVals(n_EperfMods,:) = zscore(forcePeak.(run_full_nm)(Eperf_trial_idx));
                                    case 4
                                        Eperf_modVals(n_EperfMods,:) = zscore(forcePeak_N.(run_full_nm)(Eperf_trial_idx));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                            
                            % force integral
                            if EperfModel_F_integral > 0
                                n_EperfMods = n_EperfMods + 1;
                                Eperf_modNames{n_EperfMods} = 'force integral';
                                switch EperfModel_F_integral
                                    case 1
                                        Eperf_modVals(n_EperfMods,:) = raw_or_z(AUC.(run_full_nm)(Eperf_trial_idx));
                                    case 2
                                        Eperf_modVals(n_EperfMods,:) = raw_or_z(AUC_overshoot.(run_full_nm)(Eperf_trial_idx));
                                    case 3
                                        Eperf_modVals(n_EperfMods,:) = raw_or_z(AUC_N.(run_full_nm)(Eperf_trial_idx));
                                    case 4
                                        Eperf_modVals(n_EperfMods,:) = raw_or_z(AUC_overshoot_N.(run_full_nm)(Eperf_trial_idx));
                                    case 5
                                        Eperf_modVals(n_EperfMods,:) = zscore(AUC.(run_full_nm)(Eperf_trial_idx));
                                    case 6
                                        Eperf_modVals(n_EperfMods,:) = zscore(AUC_overshoot.(run_full_nm)(Eperf_trial_idx));
                                    case 7
                                        Eperf_modVals(n_EperfMods,:) = zscore(AUC_N.(run_full_nm)(Eperf_trial_idx));
                                    case 8
                                        Eperf_modVals(n_EperfMods,:) = zscore(AUC_overshoot_N.(run_full_nm)(Eperf_trial_idx));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                        end % physical effort filter
                        
                        if strcmp(task_id,'Em')
                            % efficacy
                            if EperfModel_efficacy > 0
                                n_EperfMods = n_EperfMods + 1;
                                Eperf_modNames{n_EperfMods} = 'Em efficacy';
                                switch EperfModel_efficacy
                                    case 1
                                        Eperf_modVals(n_EperfMods,:) = raw_or_z(efficacy_with2first.(run_full_nm)(Eperf_trial_idx));
                                    case 2
                                        Eperf_modVals(n_EperfMods,:) = raw_or_z(efficacy_pureNback.(run_full_nm)(Eperf_trial_idx));
                                    case 3
                                        Eperf_modVals(n_EperfMods,:) = raw_or_z(efficacy_bis_with2first.(run_full_nm)(Eperf_trial_idx));
                                    case 4
                                        Eperf_modVals(n_EperfMods,:) = raw_or_z(efficacy_bis_pureNback.(run_full_nm)(Eperf_trial_idx));
                                    case 5
                                        Eperf_modVals(n_EperfMods,:) = zscore(efficacy_with2first.(run_full_nm)(Eperf_trial_idx));
                                    case 6
                                        Eperf_modVals(n_EperfMods,:) = zscore(efficacy_pureNback.(run_full_nm)(Eperf_trial_idx));
                                    case 7
                                        Eperf_modVals(n_EperfMods,:) = zscore(efficacy_bis_with2first.(run_full_nm)(Eperf_trial_idx));
                                    case 8
                                        Eperf_modVals(n_EperfMods,:) = zscore(efficacy_bis_pureNback.(run_full_nm)(Eperf_trial_idx));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                            
                            % RT average
                            if EperfModel_RT_avg > 0
                                n_EperfMods = n_EperfMods + 1;
                                Eperf_modNames{n_EperfMods} = 'avg RT N-back perf';
                                switch EperfModel_RT_avg
                                    case 1
                                        Eperf_modVals(n_EperfMods,:) = raw_or_z(RT_avg.(run_full_nm)(Eperf_trial_idx));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                            
                            % number of correct answers
                            if EperfModel_n_correct > 0
                                n_EperfMods = n_EperfMods + 1;
                                Eperf_modNames{n_EperfMods} = 'nb correct';
                                switch EperfModel_n_correct
                                    case 1
                                        Eperf_modVals(n_EperfMods,:) = raw_or_z(n_correct.(run_full_nm)(Eperf_trial_idx));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                            
                            % number of errors
                            if EperfModel_n_errors > 0
                                n_EperfMods = n_EperfMods + 1;
                                Eperf_modNames{n_EperfMods} = 'nb errors';
                                switch EperfModel_n_errors
                                    case 1
                                        Eperf_modVals(n_EperfMods,:) = raw_or_z(n_errors.(run_full_nm)(Eperf_trial_idx));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                        end % mental effort filter
                        
                        % net value chosen
                        if EperfModel_NV_chosen > 0
                            n_EperfMods = n_EperfMods + 1;
                            switch EperfModel_NV_chosen
                                case 1
                                    Eperf_modNames{n_EperfMods} = 'NV chosen';
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(NV_ch_min_unch.(run_full_nm)(Eperf_trial_idx));
                                case 2
                                    Eperf_modNames{n_EperfMods} = 'p(chosen)';
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(pChosen.(run_full_nm)(Eperf_trial_idx));
                                case 3
                                    Eperf_modNames{n_EperfMods} = 'NVch-NVunch';
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(NV_ch_min_unch_with_bias.(run_full_nm)(Eperf_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % net value non-default option
                        if EperfModel_NV_varOption > 0
                            n_EperfMods = n_EperfMods + 1;
                            switch EperfModel_NV_varOption
                                case 1
                                    Eperf_modNames{n_EperfMods} = 'delta NV high E - low E';
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(NV_varOption.(run_full_nm)(Eperf_trial_idx));
                                case 2
                                    Eperf_modNames{n_EperfMods} = '|delta NV high E - low E|';
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(abs(NV_varOption.(run_full_nm)(Eperf_trial_idx)));
                                case 3
                                    Eperf_modNames{n_EperfMods} = 'p(choice=hE)';
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(abs(pChoice_hE.(run_full_nm)(Eperf_trial_idx)));
                                case 4
                                    Eperf_modNames{n_EperfMods} = 'delta NV high E - low E + bias';
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(NV_varOption_plus_bias.(run_full_nm)(Eperf_trial_idx));
                                case 5
                                    Eperf_modNames{n_EperfMods} = '|delta NV high E - low E + bias|';
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(abs(NV_varOption_plus_bias.(run_full_nm)(Eperf_trial_idx)));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % net value non-default option bis
                        if EperfModel_NV_varOption_bis > 0
                            n_EperfMods = n_EperfMods + 1;
                            switch EperfModel_NV_varOption_bis
                                case 1
                                    Eperf_modNames{n_EperfMods} = 'delta NV high E - low E';
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(NV_varOption.(run_full_nm)(Eperf_trial_idx));
                                case 2
                                    Eperf_modNames{n_EperfMods} = '|delta NV high E - low E|';
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(abs(NV_varOption.(run_full_nm)(Eperf_trial_idx)));
                                case 3
                                    Eperf_modNames{n_EperfMods} = 'p(choice=hE)';
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(abs(pChoice_hE.(run_full_nm)(Eperf_trial_idx)));
                                case 4
                                    Eperf_modNames{n_EperfMods} = 'delta NV high E - low E + bias';
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(NV_varOption_plus_bias.(run_full_nm)(Eperf_trial_idx));
                                case 5
                                    Eperf_modNames{n_EperfMods} = '|delta NV high E - low E + bias|';
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(abs(NV_varOption_plus_bias.(run_full_nm)(Eperf_trial_idx)));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % RT 1st answer
                        if EperfModel_RT1stAnswer > 0
                            n_EperfMods = n_EperfMods + 1;
                            Eperf_modNames{n_EperfMods} = 'Effort latency';
                            switch EperfModel_RT1stAnswer
                                case 1
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(latency.(run_full_nm)(Eperf_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        switch task_id
                            case 'Ep'
                                % physical fatigue
                                if EperfModel_fatigue > 0
                                    n_EperfMods = n_EperfMods + 1;
                                    Eperf_modNames{n_EperfMods} = 'fatigue';
                                    switch EperfModel_fatigue
                                        case 1
                                            Eperf_modVals(n_EperfMods,:) = raw_or_z(fatigue.(run_full_nm)(Eperf_trial_idx));
                                        case 2
                                            Eperf_modVals(n_EperfMods,:) = zscore(fatigue.(run_full_nm)(Eperf_trial_idx));
                                        otherwise
                                            error('not ready yet');
                                    end
                                end
                                
                                % Ech*(physical fatigue)
                                if EperfModel_Ech_x_fatigue > 0
                                    n_EperfMods = n_EperfMods + 1;
                                    Eperf_modNames{n_EperfMods} = 'Ech_x_fatigue';
                                    switch EperfModel_Ech_x_fatigue
                                        case 1
                                            Eperf_modVals(n_EperfMods,:) = raw_or_z(E_chosen.(run_full_nm)(Eperf_trial_idx).*fatigue.(run_full_nm)(Eperf_trial_idx));
                                        otherwise
                                            error('not ready yet');
                                    end
                                end
                            case 'Em'
                                % mental facilitation based on previous trial
                                % performance
                                if EperfModel_prevEfficacy > 0
                                    n_EperfMods = n_EperfMods + 1;
                                    Eperf_modNames{n_EperfMods} = 'previous trial efficacy';
                                    switch EperfModel_prevEfficacy
                                        case 1
                                            Eperf_modVals(n_EperfMods,:) = raw_or_z(prevEfficacy_with2first.(run_full_nm)(Eperf_trial_idx));
                                        case 2
                                            Eperf_modVals(n_EperfMods,:) = raw_or_z(prevEfficacy_pureNback.(run_full_nm)(Eperf_trial_idx));
                                        case 3
                                            Eperf_modVals(n_EperfMods,:) = raw_or_z(prevEfficacy_bis_with2first.(run_full_nm)(Eperf_trial_idx));
                                        case 4
                                            Eperf_modVals(n_EperfMods,:) = raw_or_z(prevEfficacy_bis_pureNback.(run_full_nm)(Eperf_trial_idx));
                                        case 5
                                            Eperf_modVals(n_EperfMods,:) = zscore(prevEfficacy_with2first.(run_full_nm)(Eperf_trial_idx));
                                        case 6
                                            Eperf_modVals(n_EperfMods,:) = zscore(prevEfficacy_pureNback.(run_full_nm)(Eperf_trial_idx));
                                        case 7
                                            Eperf_modVals(n_EperfMods,:) = zscore(prevEfficacy_bis_with2first.(run_full_nm)(Eperf_trial_idx));
                                        case 8
                                            Eperf_modVals(n_EperfMods,:) = zscore(prevEfficacy_bis_pureNback.(run_full_nm)(Eperf_trial_idx));
                                        otherwise
                                            error('not ready yet');
                                    end
                                end
                                
                                % Echosen*(previous trial efficacy)
                                if EperfModel_Ech_x_prevEfficacy > 0
                                    n_EperfMods = n_EperfMods + 1;
                                    Eperf_modNames{n_EperfMods} = 'Ech_x_previous trial efficacy';
                                    switch EperfModel_Ech_x_prevEfficacy
                                        case 1
                                            Eperf_modVals(n_EperfMods,:) = raw_or_z(E_chosen.(run_full_nm)(Eperf_trial_idx).*prevEfficacy_with2first.(run_full_nm)(Eperf_trial_idx));
                                        case 2
                                            Eperf_modVals(n_EperfMods,:) = raw_or_z(E_chosen.(run_full_nm)(Eperf_trial_idx).*prevEfficacy_pureNback.(run_full_nm)(Eperf_trial_idx));
                                        case 3
                                            Eperf_modVals(n_EperfMods,:) = raw_or_z(E_chosen.(run_full_nm)(Eperf_trial_idx).*prevEfficacy_bis_with2first.(run_full_nm)(Eperf_trial_idx));
                                        case 4
                                            Eperf_modVals(n_EperfMods,:) = raw_or_z(E_chosen.(run_full_nm)(Eperf_trial_idx).*prevEfficacy_bis_pureNback.(run_full_nm)(Eperf_trial_idx));
                                        otherwise
                                            error('not ready yet');
                                    end
                                end
                        end % physical/mental effort filter
                        
                        % trial number
                        if EperfModel_trialN > 0
                            n_EperfMods = n_EperfMods + 1;
                            Eperf_modNames{n_EperfMods} = 'trial number';
                            switch EperfModel_trialN
                                case 1
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(trialN.(run_full_nm)(Eperf_trial_idx));
                                case 2
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(trialN_dEch.(run_full_nm)(Eperf_trial_idx));
                                case 3
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(trialN_dEnonDef_min_Edef.(run_full_nm)(Eperf_trial_idx));
                                case 4
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(trialN_dEnonDef.(run_full_nm)(Eperf_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % choice confidence
                        if EperfModel_conf > 0
                            n_EperfMods = n_EperfMods + 1;
                            Eperf_modNames{n_EperfMods} = 'confidence';
                            switch EperfModel_conf
                                case 1 % binary variable => no zscore
                                    Eperf_modVals(n_EperfMods,:) = confidence.(run_full_nm)(Eperf_trial_idx);
                                case {2,3,4} % confidence inferred by the model => ok to zscore
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(confidence.(run_full_nm)(Eperf_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        [matlabbatch] = First_level_loadEachCondition_DCM(matlabbatch, sub_idx, iCond,...
                            ['Eperf_',RP_Eperf_nm,'_',splitE_Eperf_nm], modelEperfOnset, modelEperfDur,...
                            n_EperfMods, Eperf_modNames, Eperf_modVals,...
                            orth_vars, onsets_only_GLM);
                    end % task loop
                case 3 % all crosses pooled across sessions
                    % extract trial index for the current loop
                    switch RP_Eperf_nm
                        case 'RP'
                            RPfilter_Eperf = true(1,length(double(RP_var_binary.allTrials)));
                        case 'R'
                            RPfilter_Eperf = (RP_var_binary.allTrials == 1);
                        case 'P'
                            RPfilter_Eperf = (RP_var_binary.allTrials == 0);
                    end
                    switch splitE_Eperf_nm
                        case 'E'
                            Efilter_Eperf = true(1,length(double(RP_var_binary.allTrials)));
                        case 'E1'
                            Efilter_Eperf = (E_varOption.allTrials == 1);
                        case 'E2'
                            Efilter_Eperf = (E_varOption.allTrials == 2);
                        case 'E3'
                            Efilter_Eperf = (E_varOption.allTrials == 3);
                        case 'Ech0'
                            Efilter_Eperf = (E_chosen.allTrials == 0);
                        case 'Ech1'
                            Efilter_Eperf = (E_chosen.allTrials == 1);
                        case 'Ech2'
                            Efilter_Eperf = (E_chosen.allTrials == 2);
                        case 'Ech3'
                            Efilter_Eperf = (E_chosen.allTrials == 3);
                        case 'lEch'
                            Efilter_Eperf = (choice_hE.allTrials == 0);
                        case 'hEch'
                            Efilter_Eperf = (choice_hE.allTrials == 1);
                    end
                    switch onsets_only_GLM
                        case 0
                            Eperf_trial_idx = (RPfilter_Eperf.*Efilter_Eperf.*perfOkTrials) == 1; % NEED to transform it into logical or will just focus on the first trial
                        case 1 % keep same number of trials as in choices and feedback to avoid problems in the analysis
                            Eperf_trial_idx = (RPfilter_Eperf.*Efilter_Eperf) == 1;
                    end
                    
                    %% Effort performance onset
                    iCond = iCond + 1;
                    modelEperfOnset = onsets.EperfOnsets.allTrials(Eperf_trial_idx);
                    % duration
                    switch EperfModel
                        case 'stick'
                            modelEperfDur = 0;
                        case 'boxcar'
                            modelEperfDur = durations.EperfDur.allTrials(Eperf_trial_idx);
                    end
                    
                    %% Effort performance modulators
                    n_EperfMods = 0;
                    Eperf_modNames = cell(1,1);
                    Eperf_modVals = [];
                    
                    % choice = high effort
                    switch EperfModel_choicehE
                        case 0
                        case 1
                            n_EperfMods = n_EperfMods + 1;
                            Eperf_modNames{n_EperfMods} = 'choice = high effort';
                            Eperf_modVals(n_EperfMods,:) = choice_hE.allTrials(Eperf_trial_idx); % binary variable => no zscore
                        case 2
                            n_EperfMods = n_EperfMods + 1;
                            Eperf_modNames{n_EperfMods} = 'choice = high effort';
                            Eperf_modVals(n_EperfMods,:) = choice_hE_bis.allTrials(Eperf_trial_idx); % binary variable => no zscore
                        otherwise
                            error('not ready yet');
                    end
                    
                    % money chosen
                    if EperfModel_money_chosen > 0
                        n_EperfMods = n_EperfMods + 1;
                        Eperf_modNames{n_EperfMods} = 'money chosen';
                        switch EperfModel_money_chosen
                            case 1
                                Eperf_modVals(n_EperfMods,:) = raw_or_z(money_amount_chosen.allTrials(Eperf_trial_idx));
                            case 2
                                Eperf_modVals(n_EperfMods,:) = raw_or_z(abs_money_amount_chosen.allTrials(Eperf_trial_idx));
                            case 3
                                Eperf_modVals(n_EperfMods,:) = raw_or_z(money_level_chosen.allTrials(Eperf_trial_idx));
                            case 4
                                Eperf_modVals(n_EperfMods,:) = raw_or_z(abs_money_level_chosen.allTrials(Eperf_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % effort chosen
                    if EperfModel_effort_chosen > 0
                        n_EperfMods = n_EperfMods + 1;
                        Eperf_modNames{n_EperfMods} = 'effort chosen';
                        switch EperfModel_effort_chosen
                            case 1
                                Eperf_modVals(n_EperfMods,:) = raw_or_z(E_chosen.allTrials(Eperf_trial_idx));
                            case 3
                                Eperf_modVals(n_EperfMods,:) = raw_or_z(E_chosen_bis.allTrials(Eperf_trial_idx));
                            case 4
                                Eperf_modVals(n_EperfMods,:) = zscore(E_chosen.allTrials(Eperf_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    if strcmp(task_id,'Ep')
                        % force peak
                        if EperfModel_F_peak> 0
                            n_EperfMods = n_EperfMods + 1;
                            Eperf_modNames{n_EperfMods} = 'force peak';
                            switch EperfModel_F_peak
                                case 1
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(forcePeak.allTrials(Eperf_trial_idx));
                                case 2
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(forcePeak_N.allTrials(Eperf_trial_idx));
                                case 3
                                    Eperf_modVals(n_EperfMods,:) = zscore(forcePeak.allTrials(Eperf_trial_idx));
                                case 4
                                    Eperf_modVals(n_EperfMods,:) = zscore(forcePeak_N.allTrials(Eperf_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % force integral
                        if EperfModel_F_integral > 0
                            n_EperfMods = n_EperfMods + 1;
                            Eperf_modNames{n_EperfMods} = 'force integral';
                            switch EperfModel_F_integral
                                case 1
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(AUC.allTrials(Eperf_trial_idx));
                                case 2
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(AUC_overshoot.allTrials(Eperf_trial_idx));
                                case 3
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(AUC_N.allTrials(Eperf_trial_idx));
                                case 4
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(AUC_overshoot_N.allTrials(Eperf_trial_idx));
                                case 5
                                    Eperf_modVals(n_EperfMods,:) = zscore(AUC.allTrials(Eperf_trial_idx));
                                case 6
                                    Eperf_modVals(n_EperfMods,:) = zscore(AUC_overshoot.allTrials(Eperf_trial_idx));
                                case 7
                                    Eperf_modVals(n_EperfMods,:) = zscore(AUC_N.allTrials(Eperf_trial_idx));
                                case 8
                                    Eperf_modVals(n_EperfMods,:) = zscore(AUC_overshoot_N.allTrials(Eperf_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                    end % physical effort filter
                    
                    if strcmp(task_id,'Em')
                        % efficacy
                        if EperfModel_efficacy > 0
                            n_EperfMods = n_EperfMods + 1;
                            Eperf_modNames{n_EperfMods} = 'Em efficacy';
                            switch EperfModel_efficacy
                                case 1
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(efficacy_with2first.allTrials(Eperf_trial_idx));
                                case 2
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(efficacy_pureNback.allTrials(Eperf_trial_idx));
                                case 3
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(efficacy_bis_with2first.allTrials(Eperf_trial_idx));
                                case 4
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(efficacy_bis_pureNback.allTrials(Eperf_trial_idx));
                                case 5
                                    Eperf_modVals(n_EperfMods,:) = zscore(efficacy_with2first.allTrials(Eperf_trial_idx));
                                case 6
                                    Eperf_modVals(n_EperfMods,:) = zscore(efficacy_pureNback.allTrials(Eperf_trial_idx));
                                case 7
                                    Eperf_modVals(n_EperfMods,:) = zscore(efficacy_bis_with2first.allTrials(Eperf_trial_idx));
                                case 8
                                    Eperf_modVals(n_EperfMods,:) = zscore(efficacy_bis_pureNback.allTrials(Eperf_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % RT average
                        if EperfModel_RT_avg > 0
                            n_EperfMods = n_EperfMods + 1;
                            Eperf_modNames{n_EperfMods} = 'avg RT N-back perf';
                            switch EperfModel_RT_avg
                                case 1
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(RT_avg.allTrials(Eperf_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % number of correct answers
                        if EperfModel_n_correct > 0
                            n_EperfMods = n_EperfMods + 1;
                            Eperf_modNames{n_EperfMods} = 'nb correct';
                            switch EperfModel_n_correct
                                case 1
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(n_correct.allTrials(Eperf_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % number of errors
                        if EperfModel_n_errors > 0
                            n_EperfMods = n_EperfMods + 1;
                            Eperf_modNames{n_EperfMods} = 'nb errors';
                            switch EperfModel_n_errors
                                case 1
                                    Eperf_modVals(n_EperfMods,:) = raw_or_z(n_errors.allTrials(Eperf_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                    end % mental effort filter
                    
                    % net value chosen
                    if EperfModel_NV_chosen > 0
                        n_EperfMods = n_EperfMods + 1;
                        switch EperfModel_NV_chosen
                            case 1
                                Eperf_modNames{n_EperfMods} = 'NV chosen';
                                Eperf_modVals(n_EperfMods,:) = raw_or_z(NV_ch_min_unch.allTrials(Eperf_trial_idx));
                            case 2
                                Eperf_modNames{n_EperfMods} = 'p(chosen)';
                                Eperf_modVals(n_EperfMods,:) = raw_or_z(pChosen.allTrials(Eperf_trial_idx));
                            case 3
                                Eperf_modNames{n_EperfMods} = 'NVch-NVunch';
                                Eperf_modVals(n_EperfMods,:) = raw_or_z(NV_ch_min_unch_with_bias.allTrials(Eperf_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % net value non-default option
                    if EperfModel_NV_varOption > 0
                        n_EperfMods = n_EperfMods + 1;
                        switch EperfModel_NV_varOption
                            case 1
                                Eperf_modNames{n_EperfMods} = 'delta NV high E - low E';
                                Eperf_modVals(n_EperfMods,:) = raw_or_z(NV_varOption.allTrials(Eperf_trial_idx));
                            case 2
                                Eperf_modNames{n_EperfMods} = '|delta NV high E - low E|';
                                Eperf_modVals(n_EperfMods,:) = raw_or_z(abs(NV_varOption.allTrials(Eperf_trial_idx)));
                            case 3
                                Eperf_modNames{n_EperfMods} = 'p(choice=hE)';
                                Eperf_modVals(n_EperfMods,:) = raw_or_z(abs(pChoice_hE.allTrials(Eperf_trial_idx)));
                            case 4
                                Eperf_modNames{n_EperfMods} = 'delta NV high E - low E + bias';
                                Eperf_modVals(n_EperfMods,:) = raw_or_z(NV_varOption_plus_bias.allTrials(Eperf_trial_idx));
                            case 5
                                Eperf_modNames{n_EperfMods} = '|delta NV high E - low E + bias|';
                                Eperf_modVals(n_EperfMods,:) = raw_or_z(abs(NV_varOption_plus_bias.allTrials(Eperf_trial_idx)));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % net value non-default option bis
                    if EperfModel_NV_varOption_bis > 0
                        n_EperfMods = n_EperfMods + 1;
                        switch EperfModel_NV_varOption_bis
                            case 1
                                Eperf_modNames{n_EperfMods} = 'delta NV high E - low E';
                                Eperf_modVals(n_EperfMods,:) = raw_or_z(NV_varOption.allTrials(Eperf_trial_idx));
                            case 2
                                Eperf_modNames{n_EperfMods} = '|delta NV high E - low E|';
                                Eperf_modVals(n_EperfMods,:) = raw_or_z(abs(NV_varOption.allTrials(Eperf_trial_idx)));
                            case 3
                                Eperf_modNames{n_EperfMods} = 'p(choice=hE)';
                                Eperf_modVals(n_EperfMods,:) = raw_or_z(abs(pChoice_hE.allTrials(Eperf_trial_idx)));
                            case 4
                                Eperf_modNames{n_EperfMods} = 'delta NV high E - low E + bias';
                                Eperf_modVals(n_EperfMods,:) = raw_or_z(NV_varOption_plus_bias.allTrials(Eperf_trial_idx));
                            case 5
                                Eperf_modNames{n_EperfMods} = '|delta NV high E - low E + bias|';
                                Eperf_modVals(n_EperfMods,:) = raw_or_z(abs(NV_varOption_plus_bias.allTrials(Eperf_trial_idx)));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % RT 1st answer
                    if EperfModel_RT1stAnswer > 0
                        n_EperfMods = n_EperfMods + 1;
                        Eperf_modNames{n_EperfMods} = 'Effort latency';
                        switch EperfModel_RT1stAnswer
                            case 1
                                Eperf_modVals(n_EperfMods,:) = raw_or_z(latency.allTrials(Eperf_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    switch task_id
                        case 'Ep'
                            % physical fatigue
                            if EperfModel_fatigue > 0
                                n_EperfMods = n_EperfMods + 1;
                                Eperf_modNames{n_EperfMods} = 'fatigue';
                                switch EperfModel_fatigue
                                    case 1
                                        Eperf_modVals(n_EperfMods,:) = raw_or_z(fatigue.allTrials(Eperf_trial_idx));
                                    case 2
                                        Eperf_modVals(n_EperfMods,:) = zscore(fatigue.allTrials(Eperf_trial_idx));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                            
                            % Ech*(physical fatigue)
                            if EperfModel_Ech_x_fatigue > 0
                                n_EperfMods = n_EperfMods + 1;
                                Eperf_modNames{n_EperfMods} = 'Ech_x_fatigue';
                                switch EperfModel_Ech_x_fatigue
                                    case 1
                                        Eperf_modVals(n_EperfMods,:) = raw_or_z(E_chosen.allTrials(Eperf_trial_idx).*fatigue.allTrials(Eperf_trial_idx));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                        case 'Em'
                            % mental facilitation based on previous trial
                            % performance
                            if EperfModel_prevEfficacy > 0
                                n_EperfMods = n_EperfMods + 1;
                                Eperf_modNames{n_EperfMods} = 'previous trial efficacy';
                                switch EperfModel_prevEfficacy
                                    case 1
                                        Eperf_modVals(n_EperfMods,:) = raw_or_z(prevEfficacy_with2first.allTrials(Eperf_trial_idx));
                                    case 2
                                        Eperf_modVals(n_EperfMods,:) = raw_or_z(prevEfficacy_pureNback.allTrials(Eperf_trial_idx));
                                    case 3
                                        Eperf_modVals(n_EperfMods,:) = raw_or_z(prevEfficacy_bis_with2first.allTrials(Eperf_trial_idx));
                                    case 4
                                        Eperf_modVals(n_EperfMods,:) = raw_or_z(prevEfficacy_bis_pureNback.allTrials(Eperf_trial_idx));
                                    case 5
                                        Eperf_modVals(n_EperfMods,:) = zscore(prevEfficacy_with2first.allTrials(Eperf_trial_idx));
                                    case 6
                                        Eperf_modVals(n_EperfMods,:) = zscore(prevEfficacy_pureNback.allTrials(Eperf_trial_idx));
                                    case 7
                                        Eperf_modVals(n_EperfMods,:) = zscore(prevEfficacy_bis_with2first.allTrials(Eperf_trial_idx));
                                    case 8
                                        Eperf_modVals(n_EperfMods,:) = zscore(prevEfficacy_bis_pureNback.allTrials(Eperf_trial_idx));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                            
                            % Echosen*(previous trial efficacy)
                            if EperfModel_Ech_x_prevEfficacy > 0
                                n_EperfMods = n_EperfMods + 1;
                                Eperf_modNames{n_EperfMods} = 'Ech_x_previous trial efficacy';
                                switch EperfModel_Ech_x_prevEfficacy
                                    case 1
                                        Eperf_modVals(n_EperfMods,:) = raw_or_z(E_chosen.allTrials(Eperf_trial_idx).*prevEfficacy_with2first.allTrials(Eperf_trial_idx));
                                    case 2
                                        Eperf_modVals(n_EperfMods,:) = raw_or_z(E_chosen.allTrials(Eperf_trial_idx).*prevEfficacy_pureNback.allTrials(Eperf_trial_idx));
                                    case 3
                                        Eperf_modVals(n_EperfMods,:) = raw_or_z(E_chosen.allTrials(Eperf_trial_idx).*prevEfficacy_bis_with2first.allTrials(Eperf_trial_idx));
                                    case 4
                                        Eperf_modVals(n_EperfMods,:) = raw_or_z(E_chosen.allTrials(Eperf_trial_idx).*prevEfficacy_bis_pureNback.allTrials(Eperf_trial_idx));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                    end % physical/mental effort filter
                    
                    % trial number
                    if EperfModel_trialN > 0
                        n_EperfMods = n_EperfMods + 1;
                        Eperf_modNames{n_EperfMods} = 'trial number';
                        switch EperfModel_trialN
                            case 1
                                Eperf_modVals(n_EperfMods,:) = raw_or_z(trialN.allTrials(Eperf_trial_idx));
                            case 2
                                Eperf_modVals(n_EperfMods,:) = raw_or_z(trialN_dEch.allTrials(Eperf_trial_idx));
                            case 3
                                Eperf_modVals(n_EperfMods,:) = raw_or_z(trialN_dEnonDef_min_Edef.allTrials(Eperf_trial_idx));
                            case 4
                                Eperf_modVals(n_EperfMods,:) = raw_or_z(trialN_dEnonDef.allTrials(Eperf_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % choice confidence
                    if EperfModel_conf > 0
                        n_EperfMods = n_EperfMods + 1;
                        Eperf_modNames{n_EperfMods} = 'confidence';
                        switch EperfModel_conf
                            case 1 % binary variable => no zscore
                                Eperf_modVals(n_EperfMods,:) = confidence.allTrials(Eperf_trial_idx);
                            case {2,3,4} % confidence inferred by the model => ok to zscore
                                Eperf_modVals(n_EperfMods,:) = raw_or_z(confidence.allTrials(Eperf_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    [matlabbatch] = First_level_loadEachCondition_DCM(matlabbatch, sub_idx, iCond,...
                        ['Eperf_',RP_Eperf_nm,'_',splitE_Eperf_nm], modelEperfOnset, modelEperfDur,...
                        n_EperfMods, Eperf_modNames, Eperf_modVals,...
                        orth_vars, onsets_only_GLM);
            end % DCM mode: info about how regressors concatenated or not across tasks or sessions
        end % E condition
    end % RP
end % model effort performance period

%% feedback
fbkModel = GLMprm.model_onset.(task_id).fbk;
if ismember(fbkModel,{'stick','boxcar'})
    
    for iRP_fbk = 1:length(RPfbkCond)
        RP_fbk_nm = RPfbkCond{iRP_fbk};
        for iEsplit_fbk = 1:length(EsplitFbkCond)
            splitE_fbk_nm = EsplitFbkCond{iEsplit_fbk};
            %
            fbkModel_moneyObtained  = GLMprm.fbk.(task_id).(RP_fbk_nm).(splitE_fbk_nm).money_obtained;
            fbkModel_winVSloss      = GLMprm.fbk.(task_id).(RP_fbk_nm).(splitE_fbk_nm).win_vs_loss;
            fbkModel_choicehE       = GLMprm.fbk.(task_id).(RP_fbk_nm).(splitE_fbk_nm).choiceHighE;
            fbkModel_Emade          = GLMprm.fbk.(task_id).(RP_fbk_nm).(splitE_fbk_nm).E_made;
            fbkModel_confidence     = GLMprm.fbk.(task_id).(RP_fbk_nm).(splitE_fbk_nm).confidence;
            fbkModel_trialN         = GLMprm.fbk.(task_id).(RP_fbk_nm).(splitE_fbk_nm).trialN;
            
            %% adapt depending on how sessions are pooled
            switch DCM_mode
                case 1 % all sessions independent
                    for iRun = 1:n_runs
                        run_full_nm = ['run',num2str(iRun)];
                        
                        % extract trial index for the current loop
                        switch RP_fbk_nm
                            case 'RP'
                                RPfilter_fbk = true(1,length(double(RP_var_binary.(run_full_nm))));
                            case 'R'
                                RPfilter_fbk = (RP_var_binary.(run_full_nm) == 1);
                            case 'P'
                                RPfilter_fbk = (RP_var_binary.(run_full_nm) == 0);
                        end
                        switch splitE_fbk_nm
                            case 'E'
                                Efilter_fbk = true(1,length(double(RP_var_binary.(run_full_nm))));
                            case 'E1'
                                Efilter_fbk = (E_varOption.(run_full_nm) == 1);
                            case 'E2'
                                Efilter_fbk = (E_varOption.(run_full_nm) == 2);
                            case 'E3'
                                Efilter_fbk = (E_varOption.(run_full_nm) == 3);
                            case 'Ech0'
                                Efilter_fbk = (E_chosen.(run_full_nm) == 0);
                            case 'Ech1'
                                Efilter_fbk = (E_chosen.(run_full_nm) == 1);
                            case 'Ech2'
                                Efilter_fbk = (E_chosen.(run_full_nm) == 2);
                            case 'Ech3'
                                Efilter_fbk = (E_chosen.(run_full_nm) == 3);
                            case 'lEch'
                                Efilter_fbk = (choice_hE.(run_full_nm) == 0);
                            case 'hEch'
                                Efilter_fbk = (choice_hE.(run_full_nm) == 1);
                        end
                        fbk_trial_idx = (RPfilter_fbk.*Efilter_fbk) == 1; % NEED to transform it into logical or will just focus on the first trial
                        
                        %% feedback onset
                        iCond = iCond + 1;
                        modelFbkOnset = onsets.fbkOnsets.(run_full_nm)(fbk_trial_idx);
                        % duration
                        switch fbkModel
                            case 'stick'
                                modelFbkDur = 0;
                            case 'boxcar'
                                modelFbkDur = durations.fbkDur.(run_full_nm)(fbk_trial_idx);
                        end
                        
                        %% feedback modulators
                        n_fbkMods = 0;
                        fbk_modNames = cell(1,1);
                        fbk_modVals = [];
                        
                        % win vs loss
                        if fbkModel_winVSloss > 0
                            if strcmp(RP_fbk_nm,'RP')
                                n_fbkMods = n_fbkMods + 1;
                                fbk_modNames{n_fbkMods} = 'win vs loss';
                                switch fbkModel_winVSloss % binary variable => no zscore
                                    case 1
                                        fbk_modVals(n_fbkMods,:) = win_vs_loss_fbk.(run_full_nm);
                                    otherwise
                                        error('not ready yet');
                                end
                            else
                                error('cannot split R and P trials and add a variable representing R/P trials') ;
                            end
                        end
                        
                        % choice = high effort
                        switch fbkModel_choicehE
                            case 0
                            case 1
                                n_fbkMods = n_fbkMods + 1;
                                fbk_modNames{n_fbkMods} = 'choice = high effort';
                                fbk_modVals(n_fbkMods,:) = choice_hE.(run_full_nm)(fbk_trial_idx); % binary variable => no zscore
                            case 2
                                n_fbkMods = n_fbkMods + 1;
                                fbk_modNames{n_fbkMods} = 'choice = high effort';
                                fbk_modVals(n_fbkMods,:) = choice_hE_bis.(run_full_nm)(fbk_trial_idx); % binary variable => no zscore
                            otherwise
                                error('not ready yet');
                        end
                        
                        % money chosen
                        if fbkModel_moneyObtained > 0
                            n_fbkMods = n_fbkMods + 1;
                            fbk_modNames{n_fbkMods} = 'money obtained';
                            switch fbkModel_moneyObtained
                                case 1 % money amount
                                    fbk_modVals(n_fbkMods,:) = raw_or_z(money_amount_obtained.(run_full_nm)(fbk_trial_idx));
                                case 2 % |money amount|
                                    fbk_modVals(n_fbkMods,:) = raw_or_z(abs(money_amount_obtained.(run_full_nm)(fbk_trial_idx)));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % effort performed
                        if fbkModel_Emade > 0
                            n_fbkMods = n_fbkMods + 1;
                            fbk_modNames{n_fbkMods} = 'effort made';
                            switch fbkModel_Emade
                                case 1
                                    error('not ready yet');
                                    %             fbk_modVals(n_fbkMods,:) = ;
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % trial number
                        if fbkModel_trialN > 0
                            n_fbkMods = n_fbkMods + 1;
                            fbk_modNames{n_fbkMods} = 'trial number';
                            switch fbkModel_trialN
                                case 1
                                    fbk_modVals(n_fbkMods,:) = raw_or_z(trialN.(run_full_nm)(fbk_trial_idx));
                                case 2
                                    fbk_modVals(n_fbkMods,:) = raw_or_z(trialN_dEch.(run_full_nm)(fbk_trial_idx));
                                case 3
                                    fbk_modVals(n_fbkMods,:) = raw_or_z(trialN_dEnonDef_min_Edef.(run_full_nm)(fbk_trial_idx));
                                case 4
                                    fbk_modVals(n_fbkMods,:) = raw_or_z(trialN_dEnonDef.(run_full_nm)(fbk_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % confidence
                        if fbkModel_confidence > 0
                            n_fbkMods = n_fbkMods + 1;
                            fbk_modNames{n_fbkMods} = 'confidence';
                            switch fbkModel_confidence
                                case 1 % binary variable => no zscore
                                    fbk_modVals(n_choiceMods,:) = confidence.(run_full_nm)(fbk_trial_idx);
                                case {2,3,4} % confidence inferred by the model => ok to zscore
                                    fbk_modVals(n_choiceMods,:) = raw_or_z(confidence.(run_full_nm)(fbk_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        [matlabbatch] = First_level_loadEachCondition_DCM(matlabbatch, sub_idx, iCond,...
                            ['fbk_',RP_fbk_nm,'_',splitE_fbk_nm], modelFbkOnset, modelFbkDur,...
                            n_fbkMods, fbk_modNames, fbk_modVals,...
                            orth_vars, onsets_only_GLM);
                    end % run loop
                case 2 % all tasks independent but sessions pooled
                    for iTask = 1:nTasks
                        task_nm = tasks{iTask};
                        
                        % extract trial index for the current loop
                        switch RP_fbk_nm
                            case 'RP'
                                RPfilter_fbk = true(1,length(double(RP_var_binary.(task_nm))));
                            case 'R'
                                RPfilter_fbk = (RP_var_binary.(task_nm) == 1);
                            case 'P'
                                RPfilter_fbk = (RP_var_binary.(task_nm) == 0);
                        end
                        switch splitE_fbk_nm
                            case 'E'
                                Efilter_fbk = true(1,length(double(RP_var_binary.(task_nm))));
                            case 'E1'
                                Efilter_fbk = (E_varOption.(task_nm) == 1);
                            case 'E2'
                                Efilter_fbk = (E_varOption.(task_nm) == 2);
                            case 'E3'
                                Efilter_fbk = (E_varOption.(task_nm) == 3);
                            case 'Ech0'
                                Efilter_fbk = (E_chosen.(task_nm) == 0);
                            case 'Ech1'
                                Efilter_fbk = (E_chosen.(task_nm) == 1);
                            case 'Ech2'
                                Efilter_fbk = (E_chosen.(task_nm) == 2);
                            case 'Ech3'
                                Efilter_fbk = (E_chosen.(task_nm) == 3);
                            case 'lEch'
                                Efilter_fbk = (choice_hE.(task_nm) == 0);
                            case 'hEch'
                                Efilter_fbk = (choice_hE.(task_nm) == 1);
                        end
                        fbk_trial_idx = (RPfilter_fbk.*Efilter_fbk) == 1; % NEED to transform it into logical or will just focus on the first trial
                        
                        %% feedback onset
                        iCond = iCond + 1;
                        modelFbkOnset = onsets.fbkOnsets.(task_nm)(fbk_trial_idx);
                        % duration
                        switch fbkModel
                            case 'stick'
                                modelFbkDur = 0;
                            case 'boxcar'
                                modelFbkDur = durations.fbkDur.(task_nm)(fbk_trial_idx);
                        end
                        
                        %% feedback modulators
                        n_fbkMods = 0;
                        fbk_modNames = cell(1,1);
                        fbk_modVals = [];
                        
                        % win vs loss
                        if fbkModel_winVSloss > 0
                            if strcmp(RP_fbk_nm,'RP')
                                n_fbkMods = n_fbkMods + 1;
                                fbk_modNames{n_fbkMods} = 'win vs loss';
                                switch fbkModel_winVSloss % binary variable => no zscore
                                    case 1
                                        fbk_modVals(n_fbkMods,:) = win_vs_loss_fbk.(task_nm);
                                    otherwise
                                        error('not ready yet');
                                end
                            else
                                error('cannot split R and P trials and add a variable representing R/P trials') ;
                            end
                        end
                        
                        % choice = high effort
                        switch fbkModel_choicehE
                            case 0
                            case 1
                                n_fbkMods = n_fbkMods + 1;
                                fbk_modNames{n_fbkMods} = 'choice = high effort';
                                fbk_modVals(n_fbkMods,:) = choice_hE.(task_nm)(fbk_trial_idx); % binary variable => no zscore
                            case 2
                                n_fbkMods = n_fbkMods + 1;
                                fbk_modNames{n_fbkMods} = 'choice = high effort';
                                fbk_modVals(n_fbkMods,:) = choice_hE_bis.(task_nm)(fbk_trial_idx); % binary variable => no zscore
                            otherwise
                                error('not ready yet');
                        end
                        
                        % money chosen
                        if fbkModel_moneyObtained > 0
                            n_fbkMods = n_fbkMods + 1;
                            fbk_modNames{n_fbkMods} = 'money obtained';
                            switch fbkModel_moneyObtained
                                case 1 % money amount
                                    fbk_modVals(n_fbkMods,:) = raw_or_z(money_amount_obtained.(task_nm)(fbk_trial_idx));
                                case 2 % |money amount|
                                    fbk_modVals(n_fbkMods,:) = raw_or_z(abs(money_amount_obtained.(task_nm)(fbk_trial_idx)));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % effort performed
                        if fbkModel_Emade > 0
                            n_fbkMods = n_fbkMods + 1;
                            fbk_modNames{n_fbkMods} = 'effort made';
                            switch fbkModel_Emade
                                case 1
                                    error('not ready yet');
                                    %             fbk_modVals(n_fbkMods,:) = ;
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % trial number
                        if fbkModel_trialN > 0
                            n_fbkMods = n_fbkMods + 1;
                            fbk_modNames{n_fbkMods} = 'trial number';
                            switch fbkModel_trialN
                                case 1
                                    fbk_modVals(n_fbkMods,:) = raw_or_z(trialN.(task_nm)(fbk_trial_idx));
                                case 2
                                    fbk_modVals(n_fbkMods,:) = raw_or_z(trialN_dEch.(task_nm)(fbk_trial_idx));
                                case 3
                                    fbk_modVals(n_fbkMods,:) = raw_or_z(trialN_dEnonDef_min_Edef.(task_nm)(fbk_trial_idx));
                                case 4
                                    fbk_modVals(n_fbkMods,:) = raw_or_z(trialN_dEnonDef.(task_nm)(fbk_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        % confidence
                        if fbkModel_confidence > 0
                            n_fbkMods = n_fbkMods + 1;
                            fbk_modNames{n_fbkMods} = 'confidence';
                            switch fbkModel_confidence
                                case 1 % binary variable => no zscore
                                    fbk_modVals(n_choiceMods,:) = confidence.(task_nm)(fbk_trial_idx);
                                case {2,3,4} % confidence inferred by the model => ok to zscore
                                    fbk_modVals(n_choiceMods,:) = raw_or_z(confidence.(task_nm)(fbk_trial_idx));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        [matlabbatch] = First_level_loadEachCondition_DCM(matlabbatch, sub_idx, iCond,...
                            ['fbk_',RP_fbk_nm,'_',splitE_fbk_nm], modelFbkOnset, modelFbkDur,...
                            n_fbkMods, fbk_modNames, fbk_modVals,...
                            orth_vars, onsets_only_GLM);
                    end % task loop
                case {3,4,5} % all fixation crosses pooled across sessions
                    % extract trial index for the current loop
                    switch RP_fbk_nm
                        case 'RP'
                            RPfilter_fbk = true(1,length(double(RP_var_binary.allTrials)));
                        case 'R'
                            RPfilter_fbk = (RP_var_binary.allTrials == 1);
                        case 'P'
                            RPfilter_fbk = (RP_var_binary.allTrials == 0);
                    end
                    switch splitE_fbk_nm
                        case 'E'
                            Efilter_fbk = true(1,length(double(RP_var_binary.allTrials)));
                        case 'E1'
                            Efilter_fbk = (E_varOption.allTrials == 1);
                        case 'E2'
                            Efilter_fbk = (E_varOption.allTrials == 2);
                        case 'E3'
                            Efilter_fbk = (E_varOption.allTrials == 3);
                        case 'Ech0'
                            Efilter_fbk = (E_chosen.allTrials == 0);
                        case 'Ech1'
                            Efilter_fbk = (E_chosen.allTrials == 1);
                        case 'Ech2'
                            Efilter_fbk = (E_chosen.allTrials == 2);
                        case 'Ech3'
                            Efilter_fbk = (E_chosen.allTrials == 3);
                        case 'lEch'
                            Efilter_fbk = (choice_hE.allTrials == 0);
                        case 'hEch'
                            Efilter_fbk = (choice_hE.allTrials == 1);
                    end
                    fbk_trial_idx = (RPfilter_fbk.*Efilter_fbk) == 1; % NEED to transform it into logical or will just focus on the first trial
                    
                    %% feedback onset
                    iCond = iCond + 1;
                    modelFbkOnset = onsets.fbkOnsets.allTrials(fbk_trial_idx);
                    % duration
                    switch fbkModel
                        case 'stick'
                            modelFbkDur = 0;
                        case 'boxcar'
                            modelFbkDur = durations.fbkDur.allTrials(fbk_trial_idx);
                    end
                    
                    %% feedback modulators
                    n_fbkMods = 0;
                    fbk_modNames = cell(1,1);
                    fbk_modVals = [];
                    
                    % win vs loss
                    if fbkModel_winVSloss > 0
                        if strcmp(RP_fbk_nm,'RP')
                            n_fbkMods = n_fbkMods + 1;
                            fbk_modNames{n_fbkMods} = 'win vs loss';
                            switch fbkModel_winVSloss % binary variable => no zscore
                                case 1
                                    fbk_modVals(n_fbkMods,:) = win_vs_loss_fbk.allTrials;
                                otherwise
                                    error('not ready yet');
                            end
                        else
                            error('cannot split R and P trials and add a variable representing R/P trials') ;
                        end
                    end
                    
                    % choice = high effort
                    switch fbkModel_choicehE
                        case 0
                        case 1
                            n_fbkMods = n_fbkMods + 1;
                            fbk_modNames{n_fbkMods} = 'choice = high effort';
                            fbk_modVals(n_fbkMods,:) = choice_hE.allTrials(fbk_trial_idx); % binary variable => no zscore
                        case 2
                            n_fbkMods = n_fbkMods + 1;
                            fbk_modNames{n_fbkMods} = 'choice = high effort';
                            fbk_modVals(n_fbkMods,:) = choice_hE_bis.allTrials(fbk_trial_idx); % binary variable => no zscore
                        otherwise
                            error('not ready yet');
                    end
                    
                    % money chosen
                    if fbkModel_moneyObtained > 0
                        n_fbkMods = n_fbkMods + 1;
                        fbk_modNames{n_fbkMods} = 'money obtained';
                        switch fbkModel_moneyObtained
                            case 1 % money amount
                                fbk_modVals(n_fbkMods,:) = raw_or_z(money_amount_obtained.allTrials(fbk_trial_idx));
                            case 2 % |money amount|
                                fbk_modVals(n_fbkMods,:) = raw_or_z(abs(money_amount_obtained.allTrials(fbk_trial_idx)));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % effort performed
                    if fbkModel_Emade > 0
                        n_fbkMods = n_fbkMods + 1;
                        fbk_modNames{n_fbkMods} = 'effort made';
                        switch fbkModel_Emade
                            case 1
                                error('not ready yet');
                                %             fbk_modVals(n_fbkMods,:) = ;
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % trial number
                    if fbkModel_trialN > 0
                        n_fbkMods = n_fbkMods + 1;
                        fbk_modNames{n_fbkMods} = 'trial number';
                        switch fbkModel_trialN
                            case 1
                                fbk_modVals(n_fbkMods,:) = raw_or_z(trialN.allTrials(fbk_trial_idx));
                            case 2
                                fbk_modVals(n_fbkMods,:) = raw_or_z(trialN_dEch.allTrials(fbk_trial_idx));
                            case 3
                                fbk_modVals(n_fbkMods,:) = raw_or_z(trialN_dEnonDef_min_Edef.allTrials(fbk_trial_idx));
                            case 4
                                fbk_modVals(n_fbkMods,:) = raw_or_z(trialN_dEnonDef.allTrials(fbk_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    % confidence
                    if fbkModel_confidence > 0
                        n_fbkMods = n_fbkMods + 1;
                        fbk_modNames{n_fbkMods} = 'confidence';
                        switch fbkModel_confidence
                            case 1 % binary variable => no zscore
                                fbk_modVals(n_choiceMods,:) = confidence.allTrials(fbk_trial_idx);
                            case {2,3,4} % confidence inferred by the model => ok to zscore
                                fbk_modVals(n_choiceMods,:) = raw_or_z(confidence.allTrials(fbk_trial_idx));
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    [matlabbatch] = First_level_loadEachCondition_DCM(matlabbatch, sub_idx, iCond,...
                        ['fbk_',RP_fbk_nm,'_',splitE_fbk_nm], modelFbkOnset, modelFbkDur,...
                        n_fbkMods, fbk_modNames, fbk_modVals,...
                        orth_vars, onsets_only_GLM);
            end % DCM mode: info about how regressors concatenated or not across tasks or sessions
        end % E condition
    end % RP
end % model fbk period

end % function