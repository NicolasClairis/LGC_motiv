function[matlabbatch] = First_level_loadRegressors(matlabbatch, GLMprm, study_nm, sub_nm, sub_idx, iRun, jRun,...
    subBehaviorFolder, currRunBehaviorFileName, task_fullName, computerRoot)
% [matlabbatch] = First_level_loadRegressors(matlabbatch, GLMprm, study_nm, sub_nm, sub_idx, iRun, jRun,...
%     subBehaviorFolder, currRunBehaviorFileName, task_nm, computer_root)
%
% First_level_loadRegressors will load the regressors of interest for each
% task
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
% iRun: run number for the batch
%
% jRun: actual run number (for the subject) (will be equal to iRun if no
% bugs, but may be different if some runs are ignored for the current
% subject)
%
% subBehaviorFolder: subject behavioral folder name
%
% currRunBehaviorFilename: name of the file containing the behavioral data
%
% task_fullName: task name 'mental' or 'physical'
%
% computerRoot: root of where the data is (required to be able to load net
% value if the GLM asks for it)
%
% OUTPUTS
% matlabbatch: structure updated with the regressors of interest depending
% on GLMprm values
%
% See also First_level_batch.m

%% is it a GLM for onsets-only?
onsets_only_GLM = GLMprm.gal.onsets_only;

%% determine task to work on
switch task_fullName
    case 'physical'
        task_id = 'Ep';
        task_behavioral_id = 'physicalPerf';
    case 'mental'
        task_id = 'Em';
        task_behavioral_id = 'mentalE_perf';
end

%% load data
% warning('script needs to be updated to be able to model failed trials (choice or perf) separately');
behavioralDataStruct = load([subBehaviorFolder, currRunBehaviorFileName]);

%% extract all onsets of interest
T0 = behavioralDataStruct.onsets.T0;
if strcmp(study_nm,'fMRI_pilots') && ismember(subBehaviorFolder((end-30):end),...
        {[filesep,'fMRI_pilots',filesep,'pilot_s1',filesep,'behavior',filesep],...
        [filesep,'fMRI_pilots',filesep,'pilot_s2',filesep,'behavior',filesep]})
    % for pilots s1 & s2 where there was only 1 fixation cross before
    % choice, no cross before effort
    preChoiceCrossOnsets    = behavioralDataStruct.(task_behavioral_id).onsets.cross - T0; % only cross appearing BEFORE the start of each trial
    preEffortCrossOnsets = [];
else% for all other subjects
    preChoiceCrossOnsets    = behavioralDataStruct.(task_behavioral_id).onsets.preChoiceCross - T0; % only cross appearing BEFORE the start of each trial
    preEffortCrossOnsets    = behavioralDataStruct.(task_behavioral_id).onsets.preEffortCross - T0;
end
[allCrossesOnsets, allCrossesIdx] = sort([preChoiceCrossOnsets, preEffortCrossOnsets]);
whiteCrossOnsets        = [preChoiceCrossOnsets, behavioralDataStruct.onsets.finalCross - T0]; % add final cross (at the end of the experiment)
dispChoiceOptionOnsets  = behavioralDataStruct.(task_behavioral_id).onsets.dispChoiceOptions - T0;
choiceOnsets            = behavioralDataStruct.(task_behavioral_id).onsets.choice - T0;
dispChosenOnsets        = behavioralDataStruct.(task_behavioral_id).onsets.dispChoice - T0;
n_trials = length(dispChosenOnsets);

% extract trials where no choice was made
choiceMissedTrials = isnan(choiceOnsets);

% extract onsets when the effort was performed
EperfOnsets = NaN(1,n_trials);
for iTrial = 1:n_trials
    switch task_fullName
        case 'physical'
            EperfOnsets(iTrial) = behavioralDataStruct.(task_behavioral_id).onsets.effortPeriod{1,iTrial}.effort_phase - T0;
        case 'mental'
            EperfOnsets(iTrial) = behavioralDataStruct.(task_behavioral_id).onsets.effortPeriod{1,iTrial}.nb_1 - T0;
    end
end
% feedback onsets
fbkOnsets = behavioralDataStruct.(task_behavioral_id).onsets.fbk - T0;

%% durations of each event
if strcmp(study_nm,'fMRI_pilots') &&...
        ismember(subBehaviorFolder((end-30):end),...
        {[filesep,'fMRI_pilots',filesep,'pilot_s1',filesep,'behavior',filesep],...
        [filesep,'fMRI_pilots',filesep,'pilot_s2',filesep,'behavior',filesep],...
        [filesep,'fMRI_pilots',filesep,'pilot_s3',filesep,'behavior',filesep]})
    % duration not recorded yet in those pilots
    preChoiceCrossDur = dispChoiceOptionOnsets - preChoiceCrossOnsets;
    dispChoiceOptionsDur = dispChosenOnsets - dispChoiceOptionOnsets;
    if strcmp(study_nm,'fMRI_pilots') &&...
            ismember(subBehaviorFolder((end-30):end),...
            {[filesep,'fMRI_pilots',filesep,'pilot_s1',filesep,'behavior',filesep],...
            [filesep,'fMRI_pilots',filesep,'pilot_s2',filesep,'behavior',filesep]})
        % pilot s1 & s2 no fixation cross before effort
        dispChosenDur = EperfOnsets - dispChosenOnsets;
    else
        dispChosenDur = preEffortCrossOnsets - dispChosenOnsets;
        preEffortCrossDur = EperfOnsets - preEffortCrossOnsets;
    end
    EperfDur = fbkOnsets - EperfOnsets;
    fbkDur = whiteCrossOnsets(2:end) - fbkOnsets;
else
    preChoiceCrossDur = behavioralDataStruct.(task_behavioral_id).durations.preChoiceCross;
    dispChoiceOptionsDur = behavioralDataStruct.(task_behavioral_id).durations.dispChoiceOptions;
    dispChosenDur = behavioralDataStruct.(task_behavioral_id).durations.dispChoice;
    preEffortCrossDur = behavioralDataStruct.(task_behavioral_id).durations.preEffortCross;
    EperfDur = behavioralDataStruct.(task_behavioral_id).durations.effortPeriod;
    fbkDur = behavioralDataStruct.(task_behavioral_id).durations.fbk;
end
allCrossesDur = [preChoiceCrossDur, preEffortCrossDur];
allCrossesDur = allCrossesDur(allCrossesIdx);
choice_RT = choiceOnsets - dispChoiceOptionOnsets;
% choice_RT(isnan(choice_RT))=5;

%% extract regressors of interest
% loading R vs P trial
RP_var_binary = strcmp( behavioralDataStruct.(task_behavioral_id).choiceOptions.R_or_P, 'R');
RP_var = NaN(1, length(RP_var_binary));
RP_var(RP_var_binary == 1) = 1;
RP_var(RP_var_binary == 0) = -1;
R_trials = RP_var == 1;
P_trials = RP_var == -1;
% loading default option side
defaultSide = behavioralDataStruct.(task_behavioral_id).choiceOptions.default_LR;
highE_on_the_right = defaultSide == -1;
highE_on_the_left = defaultSide == 1;
% loading money choice
money_amount_left = behavioralDataStruct.(task_behavioral_id).choiceOptions.monetary_amount.left.*RP_var;
money_amount_right = behavioralDataStruct.(task_behavioral_id).choiceOptions.monetary_amount.right.*RP_var;
money_amount_sum = money_amount_left + money_amount_right;
money_amount_varOption = money_amount_left.*highE_on_the_left + money_amount_right.*highE_on_the_right;
abs_money_amount_varOption = abs(money_amount_varOption);
money_level_left_v0 = behavioralDataStruct.(task_behavioral_id).choiceOptions.R.left;
money_level_right_v0 = behavioralDataStruct.(task_behavioral_id).choiceOptions.R.right;
money_amount_fixedOption = money_amount_left.*highE_on_the_right + money_amount_right.*highE_on_the_left;
R_amount_varOption = money_amount_varOption.*R_trials;
P_amount_varOption = money_amount_varOption.*P_trials;

% fix money levels
[money_level_left, money_level_right] = First_level_fix_money_levels(study_nm, sub_nm, task_fullName,...
    money_level_left_v0, money_level_right_v0, RP_var);

% extract other relevant variables
money_level_varOption = (money_level_left.*highE_on_the_left +...
    money_level_right.*highE_on_the_right) + 4*R_trials; % add 4 to distinguish rewards and punishments
money_level_fixedOption = money_level_left.*highE_on_the_right +...
    money_level_right.*highE_on_the_left + 4*R_trials; % add 4 to distinguish rewards and punishments

% consider that levels of reward (R1 to R4) are equivalent to levels of punishments (P1 to P4) and pool them together
abs_money_level_varOption = (money_level_left.*highE_on_the_left +...
    money_level_right.*highE_on_the_right);

R_level_varOption = (money_level_left.*highE_on_the_left +...
    money_level_right.*highE_on_the_right).*R_trials;
P_level_varOption = (money_level_left.*highE_on_the_left +...
    money_level_right.*highE_on_the_right).*P_trials;
% loading effort choice
E_left = behavioralDataStruct.(task_behavioral_id).choiceOptions.E.left;
E_right = behavioralDataStruct.(task_behavioral_id).choiceOptions.E.right;
E_sum = E_left + E_right;
E_varOption = E_left.*highE_on_the_left + E_right.*highE_on_the_right;
E_fixedOption = E_left.*highE_on_the_right + E_right.*highE_on_the_left;
% loading chosen option
choice_LRandConf = behavioralDataStruct.(task_behavioral_id).choice;
% transform into one variable equal to -1 when choice = left and 1 when
% choice = right
choice_LR = choice_LRandConf;
choice_LR(choice_LRandConf == -2) = -1;
choice_LR(choice_LRandConf == 2) = 1;
choice_left = choice_LR == -1;
choice_right = choice_LR == 1;
% choice = high effort
run_nm = num2str(jRun); % careful: use jRun for the run name
[choice_hE] = extract_choice_hE(subBehaviorFolder, sub_nm, run_nm, task_fullName); % (0) for low E chosen and (+1) for high E chosen
choice_hE_bis = (choice_hE == 1) - (choice_hE == 0); % (-1) for low E chosen and (+1) for high E chosen
% loading effort variables
E_chosen = behavioralDataStruct.(task_behavioral_id).E_chosen; % 0 for low E chosen and 1/2/3 for high effort chosen
E_chosen_bis = E_varOption.*choice_hE_bis; % -1/-2/-3 for low E chosen and +1/+2/+3 for high effort chosen
E_unchosen = E_left.*choice_right + E_right.*choice_left;
E_chosen_min_E_unchosen = E_chosen - E_unchosen;
Ech_min_Efixed = E_chosen - E_fixedOption;

% money levels
money_level_chosen = money_level_left.*choice_left + money_level_right.*choice_right + 4*R_trials;% add 4 to distinguish rewards and punishments
money_level_unchosen = money_level_left.*choice_right + money_level_right.*choice_left + 4*R_trials;% add 4 to distinguish rewards and punishments
% consider that levels of reward (R1 to R4) are equivalent to levels of punishments (P1 to P4) and pool them together
abs_money_level_chosen = money_level_left.*choice_left + money_level_right.*choice_right;
abs_money_level_unchosen = money_level_left.*choice_right + money_level_right.*choice_left;
moneyChosen_min_moneyUnchosen_level = money_level_chosen - money_level_unchosen;
absMoneyChosen_min_moneyUnchosen_level = abs_money_level_chosen - abs_money_level_unchosen;
% money chosen - money fixed option
moneyChosen_min_moneyFixed_level = money_level_chosen - money_level_fixedOption;
R_level_chosen = (money_level_left.*choice_left + money_level_right.*choice_right).*R_trials;
P_level_chosen = (money_level_left.*choice_left + money_level_right.*choice_right).*P_trials;
R_level_unchosen = (money_level_left.*choice_right + money_level_right.*choice_left).*R_trials;
P_level_unchosen = (money_level_left.*choice_right + money_level_right.*choice_left).*P_trials;

% amount variables
money_amount_chosen = behavioralDataStruct.(task_behavioral_id).R_chosen.*RP_var;
money_amount_unchosen = money_amount_left.*choice_right + money_amount_right.*choice_left;
% money chosen - money unchosen option
moneyChosen_min_moneyFixed_amount = money_amount_chosen - money_amount_fixedOption;
moneyChosen_min_moneyUnchosen_amount = money_amount_chosen - money_amount_unchosen;
absMoneyChosen_min_moneyUnchosen_amount = abs(money_amount_chosen - money_amount_unchosen);
abs_money_amount_chosen = abs(money_amount_chosen);
abs_money_amount_unchosen = abs(money_amount_unchosen);
R_amount_chosen = money_amount_chosen.*R_trials;
P_amount_chosen = money_amount_chosen.*P_trials;

% interaction terms
money_level_x_E_varOption = money_level_varOption.*E_varOption;
money_level_x_E_chosen = money_level_chosen.*E_varOption;
R_level_x_E_varOption = R_level_varOption.*E_varOption;
R_level_x_E_chosen = R_level_chosen.*E_chosen;
P_level_x_E_varOption = P_level_varOption.*E_varOption;
P_level_x_E_chosen = P_level_chosen.*E_chosen;

switch task_fullName
    case 'physical'
        [latency_tmp,...
            AUC_tmp, forcePeak_tmp, AUC_overshoot_tmp,...
            AUC_N_tmp, forcePeak_N_tmp, AUC_overshoot_N_tmp] = extract_grip_force(subBehaviorFolder, sub_nm, run_nm);
        % time to start squeezing
        latency = latency_tmp.allTrials;
        % force variables in Volts
        AUC = AUC_tmp.allTrials;
        forcePeak = forcePeak_tmp.allTrials;
        AUC_overshoot = AUC_overshoot_tmp.allTrials;
        % force variables converted in Newtons
        AUC_N = AUC_N_tmp.allTrials;
        forcePeak_N = forcePeak_N_tmp.allTrials;
        AUC_overshoot_N = AUC_overshoot_N_tmp.allTrials;
        % extract fatigue
        fatigue = NaN(1,n_trials);
        for iTrial = 1:n_trials
            if iTrial == 1
                fatigue(iTrial) = 0;
            else
                fatigue(iTrial) = fatigue(iTrial-1) + AUC(iTrial-1);
            end
        end
    case 'mental'
        [~, n_correct_tmp, n_errors_tmp, RT_avg_tmp,...
            ~,...
            ~,...
            efficacy_with2first_tmp,...
            efficacy_pureNback_tmp,...
            efficacy_bis_with2first_tmp,...
            efficacy_bis_pureNback_tmp,...
            latency_tmp] = extract_mental_perf(subBehaviorFolder, sub_nm, run_nm);
        n_correct = n_correct_tmp.allTrials;
        n_errors = n_errors_tmp.allTrials;
        RT_avg = RT_avg_tmp.allTrials;
        efficacy_with2first = efficacy_with2first_tmp.allTrials;
        efficacy_pureNback = efficacy_pureNback_tmp.allTrials;
        efficacy_bis_with2first = efficacy_bis_with2first_tmp.allTrials;
        efficacy_bis_pureNback = efficacy_bis_pureNback_tmp.allTrials;
        [prevEfficacy_with2first,...
            prevEfficacy_pureNback,...
            prevEfficacy_bis_with2first,...
            prevEfficacy_bis_pureNback] = deal(NaN(1,n_trials));
        for iTrial = 1:n_trials
            if iTrial == 1
                prevEfficacy_with2first(iTrial) = 0;
                prevEfficacy_pureNback(iTrial) = 0;
                prevEfficacy_bis_with2first(iTrial) = 0;
                prevEfficacy_bis_pureNback(iTrial) = 0;
            else
                prevEfficacy_with2first(iTrial) = efficacy_with2first(iTrial-1);
                prevEfficacy_pureNback(iTrial) = efficacy_pureNback(iTrial-1);
                prevEfficacy_bis_with2first(iTrial) = efficacy_bis_with2first(iTrial-1);
                prevEfficacy_bis_pureNback(iTrial) = efficacy_bis_pureNback(iTrial-1);
            end
        end
        latency = latency_tmp.allTrials;
end

% loading feedback
money_amount_obtained = behavioralDataStruct.(task_behavioral_id).gain; % could be different from reward chosen (in case of failure) but mostly similar
win_vs_loss_fbk = money_amount_obtained > 0;
% load run name
switch jRun
    case {1,2}
        run_nm_bis = 'run1';
    case {3,4}
        run_nm_bis = 'run2';
    otherwise
        error('case not ready yet: net value in the model and more than 4 runs.');
end
% load net value and confidence
RPconds = {'R','P','RP'};
EsplitConditions = {'E','E1','E2','E3',...
    'Ech0','Ech1','Ech2','Ech3',...
    'lEch','hEch'};
for iRP = 1:length(RPconds)
    RP_nm = RPconds{iRP};
    
    for iE = 1:length(EsplitConditions)
        Esplit_nm = EsplitConditions{iE};
        
        % load net value
        if GLMprm.choice.(task_id).(RP_nm).(Esplit_nm).NV_chosen > 0 ||...
                GLMprm.choice.(task_id).(RP_nm).(Esplit_nm).NV_varOption > 0 ||...
                GLMprm.chosen.(task_id).(RP_nm).(Esplit_nm).NV_chosen > 0 ||...
                GLMprm.chosen.(task_id).(RP_nm).(Esplit_nm).NV_varOption > 0 ||...
                GLMprm.Eperf.(task_id).(RP_nm).(Esplit_nm).NV_chosen > 0 ||...
                GLMprm.Eperf.(task_id).(RP_nm).(Esplit_nm).NV_varOption > 0
            
            % extract NV model name
            if GLMprm.choice.(task_id).(RP_nm).(Esplit_nm).NV_chosen > 0 ||...
                    GLMprm.choice.(task_id).(RP_nm).(Esplit_nm).NV_varOption > 0
                NV_mdl_nm = GLMprm.choice.(task_id).(RP_nm).(Esplit_nm).NV_mdl;
            elseif GLMprm.chosen.(task_id).(RP_nm).(Esplit_nm).NV_chosen > 0 ||...
                    GLMprm.chosen.(task_id).(RP_nm).(Esplit_nm).NV_varOption > 0
                NV_mdl_nm = GLMprm.chosen.(task_id).(RP_nm).(Esplit_nm).NV_mdl;
            elseif GLMprm.Eperf.(task_id).(RP_nm).(Esplit_nm).NV_chosen > 0 ||...
                    GLMprm.Eperf.(task_id).(RP_nm).(Esplit_nm).NV_varOption > 0
                NV_mdl_nm = GLMprm.Eperf.(task_id).(RP_nm).(Esplit_nm).NV_mdl;
            end
            % extract net value
            if strcmp(NV_mdl_nm(1:4),'mdl_') % classic model
                % load net value
                [~, modelledDataStruct] = logitfit_choices(computerRoot, study_nm, sub_nm,...
                    0, 'levels', 6, 6);
                NV_ch_min_unch = modelledDataStruct.NV_chosen.(task_id).(NV_mdl_nm).(run_nm_bis);
                NV_varOption = modelledDataStruct.NV_varOption.(task_id).(NV_mdl_nm).(run_nm_bis);
                pChoice_hE = modelledDataStruct.pChoice_hE.(task_id).(NV_mdl_nm).(run_nm_bis);
            elseif strcmp(NV_mdl_nm(1:14),'bayesianModel_') % bayesian model
                bayesianMdl_nm = strrep(NV_mdl_nm,'bayesianModel','mdl');
                gitResultsFolder = [fullfile('C:','Users','clairis','Desktop',...
                    'GitHub','LGC_motiv','LGC_Motiv_results',study_nm,'bayesian_modeling'),filesep];
                [NV_ch_min_unch, NV_varOption,~,pChoice_hE,...
                    NV_varOption_plus_bias, NV_ch_min_unch_with_bias] = extract_bayesian_mdl(gitResultsFolder, subBehaviorFolder,...
                    sub_nm, run_nm, task_fullName, bayesianMdl_nm);
            else
                error(['model with ',NV_mdl_nm,' not ready yet']);
            end
            
            % probability of chosen option
            pChosen = pChoice_hE.*(choice_hE == 1)' + (1 - pChoice_hE).*(choice_hE == 0)';
        end % net value
        
        % load confidence
        if GLMprm.choice.(task_id).(RP_nm).(Esplit_nm).confidence > 0 ||...
                GLMprm.chosen.(task_id).(RP_nm).(Esplit_nm).confidence > 0 ||...
                GLMprm.Eperf.(task_id).(RP_nm).(Esplit_nm).confidence > 0 ||...
                GLMprm.fbk.(task_id).(RP_nm).(Esplit_nm).confidence > 0
            
            % use confidence based on ratings
            if GLMprm.choice.(task_id).(RP_nm).(Esplit_nm).confidence == 1 ||...
                    GLMprm.chosen.(task_id).(RP_nm).(Esplit_nm).confidence == 1 ||...
                    GLMprm.Eperf.(task_id).(RP_nm).(Esplit_nm).confidence == 1 ||...
                    GLMprm.fbk.(task_id).(RP_nm).(Esplit_nm).confidence == 1
                confidence = abs(choice_LRandConf) == 2; % 0 when low confidence and 1 when high confidence
            else
                % extract NV model name
                if GLMprm.choice.(task_id).(RP_nm).(Esplit_nm).confidence > 1
                    conf_mdl_nm = GLMprm.choice.(task_id).(RP_nm).(Esplit_nm).conf_mdl;
                elseif GLMprm.chosen.(task_id).(RP_nm).(Esplit_nm).confidence > 1
                    conf_mdl_nm = GLMprm.chosen.(task_id).(RP_nm).(Esplit_nm).conf_mdl;
                elseif GLMprm.Eperf.(task_id).(RP_nm).(Esplit_nm).confidence > 1
                    conf_mdl_nm = GLMprm.Eperf.(task_id).(RP_nm).(Esplit_nm).conf_mdl;
                elseif GLMprm.fbk.(task_id).(RP_nm).(Esplit_nm).confidence > 1
                    conf_mdl_nm = GLMprm.fbk.(task_id).(RP_nm).(Esplit_nm).conf_mdl;
                end
                % extract confidence
                if strcmp(conf_mdl_nm(1:4),'mdl_') % classic model
                    % load inferred confidence
                    [~, modelledDataStruct] = logitfit_choices(computerRoot, study_nm, sub_nm,...
                        0, 'levels', 6, 6);
                    confidence_highE = modelledDataStruct.confidenceFitted_highE.(conf_mdl_nm).(task_id).(run_nm_bis);
                    pChoice = modelledDataStruct.choicesFitted.(conf_mdl_nm).(task_id).(run_nm_bis);
                elseif strcmp(conf_mdl_nm(1:14),'bayesianModel_') % bayesian model
                    bayesianMdl_nm = strrep(conf_mdl_nm,'bayesianModel','mdl');
                    gitResultsFolder = [fullfile('C:','Users','clairis','Desktop',...
                        'GitHub','LGC_motiv','LGC_Motiv_results',study_nm,'bayesian_modeling'),filesep];
                    [~, ~, confidence_highE, pChoice] = extract_bayesian_mdl(gitResultsFolder, subBehaviorFolder,...
                        sub_nm, run_nm, task_fullName, bayesianMdl_nm);
                else
                    error(['model with ',conf_mdl_nm,' not ready yet']);
                end
                
                % define confidence proxy to use
                if GLMprm.choice.(task_id).(RP_nm).(Esplit_nm).confidence == 2 ||...
                        GLMprm.chosen.(task_id).(RP_nm).(Esplit_nm).confidence == 2 ||...
                        GLMprm.Eperf.(task_id).(RP_nm).(Esplit_nm).confidence == 2 ||...
                        GLMprm.fbk.(task_id).(RP_nm).(Esplit_nm).confidence == 2
                    confidence = confidence_highE;
                elseif GLMprm.choice.(task_id).(RP_nm).(Esplit_nm).confidence == 3 ||...
                        GLMprm.chosen.(task_id).(RP_nm).(Esplit_nm).confidence == 3 ||...
                        GLMprm.Eperf.(task_id).(RP_nm).(Esplit_nm).confidence == 3 ||...
                        GLMprm.fbk.(task_id).(RP_nm).(Esplit_nm).confidence == 3
                    confidence_left = (pChoice.*highE_on_the_left + (1 - pChoice).*highE_on_the_right - 0.5).^2;
                    confidence = confidence_left./0.25; % normalise to be between 0 and 1
                elseif GLMprm.choice.(task_id).(RP_nm).(Esplit_nm).confidence == 4 ||...
                        GLMprm.chosen.(task_id).(RP_nm).(Esplit_nm).confidence == 4 ||...
                        GLMprm.Eperf.(task_id).(RP_nm).(Esplit_nm).confidence == 4 ||...
                        GLMprm.fbk.(task_id).(RP_nm).(Esplit_nm).confidence == 4
                    confidence = ((pChoice - 0.5).^2)./0.25;
                end
                
            end % which confidence to use
        end % confidence
    end % Effort conditions
end % R/P/RP loop

% extract fatigue
n_theoreticalTrialsPerRun = 54;
trialN = (1:n_theoreticalTrialsPerRun) - 1; % trial number - 1 to reflect fatigue level
trialN_dEch = trialN.*E_chosen_min_E_unchosen;
trialN_dEnonDef_min_Edef = trialN.*(E_varOption - E_fixedOption);
trialN_dEnonDef = trialN.*E_varOption;

%% remove trials where no choice was performed
if sum(choiceMissedTrials) > 0
    % note: allCross not cleaned from missed choices on purpose as it would
    % make no sense
    
    % extract onsets and durations of missed trials
    preChoiceCrossOnsets_missedOnsets = preChoiceCrossOnsets(choiceMissedTrials);
    dispChoiceOption_missedOnsets = dispChoiceOptionOnsets(choiceMissedTrials);
    choice_missedOnsets = choiceOnsets(choiceMissedTrials);
    dispChosen_missedOnsets = dispChosenOnsets(choiceMissedTrials);
    preEffortCross_missedOnsets = preEffortCrossOnsets(choiceMissedTrials);
    Eperf_missedOnsets = EperfOnsets(choiceMissedTrials);
    fbk_missedOnsets = fbkOnsets(choiceMissedTrials);
    % durations
    preChoiceCross_missedDur = preChoiceCrossDur(choiceMissedTrials);
    dispChoiceOptions_missedDur = dispChoiceOptionsDur(choiceMissedTrials);
    dispChosen_missedDur = dispChosenDur(choiceMissedTrials);
    preEffortCross_missedDur = preEffortCrossDur(choiceMissedTrials);
    Eperf_missedDur = EperfDur(choiceMissedTrials);
    fbk_missedDur = fbkDur(choiceMissedTrials);
    
    % onsets
    preChoiceCrossOnsets(choiceMissedTrials) = [];
    dispChoiceOptionOnsets(choiceMissedTrials) = [];
    choiceOnsets(choiceMissedTrials) = [];
    dispChosenOnsets(choiceMissedTrials) = [];
    preEffortCrossOnsets(choiceMissedTrials) = [];
    EperfOnsets(choiceMissedTrials) = [];
    fbkOnsets(choiceMissedTrials) = [];
    % durations
    preChoiceCrossDur(choiceMissedTrials) = [];
    dispChoiceOptionsDur(choiceMissedTrials) = [];
    dispChosenDur(choiceMissedTrials) = [];
    preEffortCrossDur(choiceMissedTrials) = [];
    EperfDur(choiceMissedTrials) = [];
    fbkDur(choiceMissedTrials) = [];
    % regressors
    R_trials(choiceMissedTrials) = [];
    P_trials(choiceMissedTrials) = [];
    defaultSide(choiceMissedTrials) = [];
    RP_var_binary(choiceMissedTrials) = [];
    choice_hE(choiceMissedTrials) = [];
    choice_hE_bis(choiceMissedTrials) = [];
    money_level_chosen(choiceMissedTrials) = [];
    R_level_chosen(choiceMissedTrials) = [];
    P_level_chosen(choiceMissedTrials) = [];
    R_level_unchosen(choiceMissedTrials) = [];
    P_level_unchosen(choiceMissedTrials) = [];
    money_amount_left(choiceMissedTrials) = [];
    money_amount_right(choiceMissedTrials) = [];
    money_amount_sum(choiceMissedTrials) = [];
    money_amount_varOption(choiceMissedTrials) = [];
    R_amount_varOption(choiceMissedTrials) = [];
    P_amount_varOption(choiceMissedTrials) = [];
    money_amount_fixedOption(choiceMissedTrials) = [];
    money_level_left(choiceMissedTrials) = [];
    money_level_right(choiceMissedTrials) = [];
    money_level_varOption(choiceMissedTrials) = [];
    R_level_varOption(choiceMissedTrials) = [];
    P_level_varOption(choiceMissedTrials) = [];
    abs_money_amount_varOption(choiceMissedTrials) = [];
    abs_money_level_varOption(choiceMissedTrials) = [];
    E_left(choiceMissedTrials) = [];
    E_right(choiceMissedTrials) = [];
    E_sum(choiceMissedTrials) = [];
    E_varOption(choiceMissedTrials) = [];
    money_level_x_E_varOption(choiceMissedTrials) = [];
    money_level_x_E_chosen(choiceMissedTrials) = [];
    R_level_x_E_varOption(choiceMissedTrials) = [];
    R_level_x_E_chosen(choiceMissedTrials) = [];
    P_level_x_E_varOption(choiceMissedTrials) = [];
    P_level_x_E_chosen(choiceMissedTrials) = [];
    choice_LRandConf(choiceMissedTrials) = [];
    choice_LR(choiceMissedTrials) = [];
    if exist('confidence','var') && ~isempty(confidence)
        confidence(choiceMissedTrials) = [];
    end
    E_chosen(choiceMissedTrials) = [];
    E_chosen_bis(choiceMissedTrials) = [];
    E_unchosen(choiceMissedTrials) = [];
    E_chosen_min_E_unchosen(choiceMissedTrials) = [];
    Ech_min_Efixed(choiceMissedTrials) = [];
    money_amount_chosen(choiceMissedTrials) = [];
    R_amount_chosen(choiceMissedTrials) = [];
    P_amount_chosen(choiceMissedTrials) = [];
    abs_money_amount_chosen(choiceMissedTrials) = [];
    abs_money_level_chosen(choiceMissedTrials) = [];
    abs_money_amount_unchosen(choiceMissedTrials) = [];
    abs_money_level_unchosen(choiceMissedTrials) = [];
    money_amount_unchosen(choiceMissedTrials) = [];
    money_level_unchosen(choiceMissedTrials) = [];
    moneyChosen_min_moneyUnchosen_amount(choiceMissedTrials) = [];
    absMoneyChosen_min_moneyUnchosen_amount(choiceMissedTrials) = [];
    moneyChosen_min_moneyUnchosen_level(choiceMissedTrials) = [];
    absMoneyChosen_min_moneyUnchosen_level(choiceMissedTrials) = [];
    moneyChosen_min_moneyFixed_amount(choiceMissedTrials) = [];
    moneyChosen_min_moneyFixed_level(choiceMissedTrials) = [];
    money_amount_obtained(choiceMissedTrials) = [];
    win_vs_loss_fbk(choiceMissedTrials) = [];
    choice_RT(choiceMissedTrials) = [];
    if exist('NV_mdl_nm','var') &&...
            (strcmp(NV_mdl_nm(1:4),'mdl_') ||...
            strcmp(NV_mdl_nm(1:14),'bayesianModel_'))
        NV_ch_min_unch(choiceMissedTrials) = [];
        NV_ch_min_unch_with_bias(choiceMissedTrials) = [];
        NV_varOption(choiceMissedTrials) = [];
        NV_varOption_plus_bias(choiceMissedTrials) = [];
        pChoice_hE(choiceMissedTrials) = [];
        pChosen(choiceMissedTrials) = [];
    end
    trialN(choiceMissedTrials) = [];
    trialN_dEch(choiceMissedTrials) = [];
    trialN_dEnonDef_min_Edef(choiceMissedTrials) = [];
    trialN_dEnonDef(choiceMissedTrials) = [];
    switch task_fullName
        case 'physical'
            AUC(choiceMissedTrials) = [];
            forcePeak(choiceMissedTrials) = [];
            AUC_overshoot(choiceMissedTrials) = [];
            AUC_N(choiceMissedTrials) = [];
            forcePeak_N(choiceMissedTrials) = [];
            AUC_overshoot_N(choiceMissedTrials) = [];
            fatigue(choiceMissedTrials) = [];
        case 'mental'
            n_correct(choiceMissedTrials) = [];
            n_errors(choiceMissedTrials) = [];
            RT_avg(choiceMissedTrials) = [];
            efficacy_with2first(choiceMissedTrials) = [];
            efficacy_pureNback(choiceMissedTrials) = [];
            efficacy_bis_with2first(choiceMissedTrials) = [];
            efficacy_bis_pureNback(choiceMissedTrials) = [];
            prevEfficacy_with2first(choiceMissedTrials) = [];
            prevEfficacy_pureNback(choiceMissedTrials) = [];
            prevEfficacy_bis_with2first(choiceMissedTrials) = [];
            prevEfficacy_bis_pureNback(choiceMissedTrials) = [];
    end
    latency(choiceMissedTrials) = [];
end

% zscored R/P ignoring P/R trials (respectively)
[z_R_level_chosen, z_P_level_chosen,...
    z_R_level_unchosen, z_P_level_unchosen,...
    z_R_amount_chosen, z_P_amount_chosen,...
    z_R_level_varOption, z_P_level_varOption,...
    z_R_amount_varOption, z_P_amount_varOption] = deal(zeros(size(R_level_chosen)));
z_R_level_chosen(R_trials) = zscore(R_level_chosen(R_trials));
z_P_level_chosen(P_trials) = zscore(P_level_chosen(P_trials));
z_R_level_unchosen(R_trials) = zscore(R_level_unchosen(R_trials));
z_P_level_unchosen(P_trials) = zscore(P_level_unchosen(P_trials));
z_R_amount_chosen(R_trials) = zscore(R_amount_chosen(R_trials));
z_P_amount_chosen(P_trials) = zscore(P_amount_chosen(P_trials));
z_R_level_varOption(R_trials) = zscore(R_level_varOption(R_trials));
z_P_level_varOption(P_trials) = zscore(P_level_varOption(P_trials));
z_R_amount_varOption(R_trials) = zscore(R_amount_varOption(R_trials));
z_P_amount_varOption(P_trials) = zscore(P_amount_varOption(P_trials));

% extract trials where no performance was achieved
perfOkTrials = ~isnan(latency);

%% load the batch according to GLMprm variables
% load general GLM parameters
orth_vars = GLMprm.gal.orth_vars;
zPerRun = GLMprm.gal.zPerRun;
switch zPerRun
    case 0 % raw data = no change
        raw_or_z = @(x) x;
    case 1 % zscore all variables per run
        raw_or_z = @(x) zscore(x);
end

%% GLM parameters for choice Reward/Punishment pool
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

%% GLM parameters for choice effort pool
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

%% initialize conditions
iCond = 0;

%% all fixation cross
allCrossesModel = GLMprm.model_onset.(task_id).allCrosses;
if ismember(allCrossesModel,{'stick','boxcar'})
    iCond = iCond + 1;
    switch allCrossesModel
        case 'stick'
            modelAllCrossesdur = 0;
        case 'boxcar'
            modelAllCrossesdur = allCrossesDur;
    end
    
    %% pre-choice cross modulators
    n_allCrossesMods = 0;
    allCrosses_modNames = cell(1,1);
    allCrosses_modVals = [];
    
    %% load all
    [matlabbatch] = First_level_loadEachCondition(matlabbatch, sub_idx, iRun, iCond,...
        'fixation cross', allCrossesOnsets, modelAllCrossesdur,...
        n_allCrossesMods, allCrosses_modNames, allCrosses_modVals,...
        orth_vars, onsets_only_GLM);
end

%% fixation cross before choice
preChoiceCrossModel = GLMprm.model_onset.(task_id).preChoiceCross;
preChoiceCrossModel_RT = GLMprm.preChoiceCross.(task_id).RT;
preChoiceCrossModel_choiceHighE = GLMprm.preChoiceCross.(task_id).choiceHighE;
preChoiceCrossModel_E_chosen    = GLMprm.preChoiceCross.(task_id).E_chosen;
if ismember(preChoiceCrossModel,{'stick','boxcar','boxcar_bis'})
    iCond = iCond + 1;
    switch preChoiceCrossModel
        case 'stick'
            modelPreChoiceCrossdur = 0;
        case 'boxcar'
            modelPreChoiceCrossdur = preChoiceCrossDur;
        case 'boxcar_bis' % go from fixation cross display until choice is made
            modelPreChoiceCrossdur = preChoiceCrossDur + choice_RT;
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
                preChoiceCross_modVals(n_preChoiceCrossMods,:) = raw_or_z(choice_RT);
            case 5
                preChoiceCross_modVals(n_preChoiceCrossMods,:) = zscore(choice_RT);
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
            preChoiceCross_modVals(n_preChoiceCrossMods,:) = choice_hE; % binary variable => no zscore
        case 2
            n_preChoiceCrossMods = n_preChoiceCrossMods + 1;
            preChoiceCross_modNames{n_preChoiceCrossMods} = 'choice = high effort';
            preChoiceCross_modVals(n_preChoiceCrossMods,:) = choice_hE_bis; % binary variable => no zscore
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
                preChoiceCross_modVals(n_preChoiceCrossMods,:) = raw_or_z(E_chosen);
            case 3
                preChoiceCross_modVals(n_preChoiceCrossMods,:) = raw_or_z(E_chosen_bis);
            case 4
                preChoiceCross_modVals(n_preChoiceCrossMods,:) = zscore(E_chosen);
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
                preChoiceCross_modVals(n_preChoiceCrossMods,:) = raw_or_z(choice_RT);
            case 2
                preChoiceCross_modVals(n_preChoiceCrossMods,:) = zscore(choice_RT);
            otherwise
                error('not ready yet');
        end
    end
    
    %% load all
    [matlabbatch] = First_level_loadEachCondition(matlabbatch, sub_idx, iRun, iCond,...
        'preChoice fixation cross', preChoiceCrossOnsets, modelPreChoiceCrossdur,...
        n_preChoiceCrossMods, preChoiceCross_modNames, preChoiceCross_modVals,...
        orth_vars, onsets_only_GLM);
end


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
            
            % extract trial index for the current loop
            switch RP_dispChoice_nm
                case 'RP'
                    RPfilter_dispChoice = true(1,length(double(RP_var_binary)));
                case 'R'
                    RPfilter_dispChoice = (RP_var_binary == 1);
                case 'P'
                    RPfilter_dispChoice = (RP_var_binary == 0);
            end
            switch splitE_dispChoice_nm
                case 'E'
                    Efilter_dispChoice = true(1,length(double(RP_var_binary)));
                case 'E1'
                    Efilter_dispChoice = (E_varOption == 1);
                case 'E2'
                    Efilter_dispChoice = (E_varOption == 2);
                case 'E3'
                    Efilter_dispChoice = (E_varOption == 3);
                case 'Ech0'
                    Efilter_dispChoice = (E_chosen == 0);
                case 'Ech1'
                    Efilter_dispChoice = (E_chosen == 1);
                case 'Ech2'
                    Efilter_dispChoice = (E_chosen == 2);
                case 'Ech3'
                    Efilter_dispChoice = (E_chosen == 3);
                case 'lEch'
                    Efilter_dispChoice = (choice_hE == 0);
                case 'hEch'
                    Efilter_dispChoice = (choice_hE == 1);
            end
            choice_trial_idx = (RPfilter_dispChoice.*Efilter_dispChoice) == 1; % NEED to transform it into logical or will just focus on the first trial
            
            %% choice onset
            iCond = iCond + 1;
            modelChoiceOnset = dispChoiceOptionOnsets(choice_trial_idx);
            % duration
            switch choiceModel
                case 'stick'
                    modelChoiceDur = 0;
                case 'boxcar'
                    modelChoiceDur = dispChoiceOptionsDur(choice_trial_idx);
                case 'boxcar_bis'
                    modelChoiceDur = dispChoiceOptionsDur(choice_trial_idx) + dispChosenDur(choice_trial_idx);
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
                        choice_modVals(n_choiceMods,:) = raw_or_z(choice_RT(choice_trial_idx));
                    case 5
                        choice_modVals(n_choiceMods,:) = zscore(choice_RT(choice_trial_idx));
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
                        choice_modVals(n_choiceMods,:) = raw_or_z(NV_ch_min_unch(choice_trial_idx));
                    case 5
                        choice_modNames{n_choiceMods} = 'p(chosen)';
                        choice_modVals(n_choiceMods,:) = raw_or_z(pChosen(choice_trial_idx));
                    case 6
                        choice_modNames{n_choiceMods} = 'NVch-NVunch';
                        choice_modVals(n_choiceMods,:) = raw_or_z(NV_ch_min_unch_with_bias(choice_trial_idx));
                    otherwise
                        error('not ready yet');
                end
            end
            
            % reward vs punishments
            if choiceModel_RP == 1
                if strcmp(RP_dispChoice_nm,'RP')
                    n_choiceMods = n_choiceMods + 1;
                    choice_modNames{n_choiceMods} = 'R vs P';
                    choice_modVals(n_choiceMods,:) = RP_var_binary(choice_trial_idx); % binary variable => no zscore
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
                    choice_modVals(n_choiceMods,:) = choice_hE(choice_trial_idx); % binary variable => no zscore
                case 2
                    n_choiceMods = n_choiceMods + 1;
                    choice_modNames{n_choiceMods} = 'choice = high effort';
                    choice_modVals(n_choiceMods,:) = choice_hE_bis(choice_trial_idx); % binary variable => no zscore
                otherwise
                    error('case not ready yet');
            end
            
            % reward for the non-default option
            if choiceModel_R_varOption > 0
                n_choiceMods = n_choiceMods + 1;
                choice_modNames{n_choiceMods} = 'R non-default';
                switch choiceModel_R_varOption
                    case 1 % R high effort option (amount)
                        choice_modVals(n_choiceMods,:) = raw_or_z(R_amount_varOption(choice_trial_idx));
                    case 2 % R high effort option (level)
                        choice_modVals(n_choiceMods,:) = raw_or_z(R_level_varOption(choice_trial_idx));
                    case 3 % z(R high effort option) (amount)
                        choice_modVals(n_choiceMods,:) = raw_or_z(z_R_amount_varOption(choice_trial_idx));
                    case 4 % z(R high effort option) (level)
                        choice_modVals(n_choiceMods,:) = raw_or_z(z_R_level_varOption(choice_trial_idx));
                end
            end
            
            % reward for the chosen option
            if choiceModel_R_chosen > 0
                n_choiceMods = n_choiceMods + 1;
                choice_modNames{n_choiceMods} = 'R chosen';
                switch choiceModel_R_chosen
                    case 1 % R chosen option amount
                        choice_modVals(n_choiceMods,:) = raw_or_z(R_amount_chosen(choice_trial_idx));
                    case 2 % R chosen option level (0/1/2/3)
                        choice_modVals(n_choiceMods,:) = raw_or_z(R_level_chosen(choice_trial_idx));
                    case 3 % z(R chosen option) amount
                        choice_modVals(n_choiceMods,:) = raw_or_z(z_R_amount_chosen(choice_trial_idx));
                    case 4 % z(R chosen option) level (0/1/2/3)
                        choice_modVals(n_choiceMods,:) = raw_or_z(z_R_level_chosen(choice_trial_idx));
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
                        choice_modVals(n_choiceMods,:) = raw_or_z(P_amount_varOption(choice_trial_idx));
                    case 2 % R high effort option (level)
                        choice_modVals(n_choiceMods,:) = raw_or_z(P_level_varOption(choice_trial_idx));
                    case 3 % z(R high effort option) (amount)
                        choice_modVals(n_choiceMods,:) = raw_or_z(z_P_amount_varOption(choice_trial_idx));
                    case 4 % z(R high effort option) (level)
                        choice_modVals(n_choiceMods,:) = raw_or_z(z_P_level_varOption(choice_trial_idx));
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
                        choice_modVals(n_choiceMods,:) = raw_or_z(P_amount_chosen(choice_trial_idx));
                    case 2 % R high effort option (level)
                        choice_modVals(n_choiceMods,:) = raw_or_z(P_level_chosen(choice_trial_idx));
                    case 3 % z(R high effort option) (amount)
                        choice_modVals(n_choiceMods,:) = raw_or_z(z_P_amount_chosen(choice_trial_idx));
                    case 4 % z(R high effort option) (level)
                        choice_modVals(n_choiceMods,:) = raw_or_z(z_P_level_chosen(choice_trial_idx));
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
                        choice_modVals(n_choiceMods,:) = raw_or_z(money_amount_left(choice_trial_idx));
                    case 2
                        choice_modVals(n_choiceMods,:) = raw_or_z(abs(money_amount_left(choice_trial_idx)));
                    case 3
                        choice_modVals(n_choiceMods,:) = raw_or_z(money_level_left(choice_trial_idx));
                    case 4
                        choice_modVals(n_choiceMods,:) = raw_or_z(abs(money_level_left(choice_trial_idx)));
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
                        choice_modVals(n_choiceMods,:) = raw_or_z(money_amount_right(choice_trial_idx));
                    case 2
                        choice_modVals(n_choiceMods,:) = raw_or_z(abs(money_amount_right(choice_trial_idx)));
                    case 3
                        choice_modVals(n_choiceMods,:) = raw_or_z(money_level_right(choice_trial_idx));
                    case 4
                        choice_modVals(n_choiceMods,:) = raw_or_z(abs(money_level_right(choice_trial_idx)));
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
                        choice_modVals(n_choiceMods,:) = raw_or_z(money_amount_chosen(choice_trial_idx));
                    case 2
                        choice_modVals(n_choiceMods,:) = raw_or_z(abs_money_amount_chosen(choice_trial_idx));
                    case 3
                        choice_modVals(n_choiceMods,:) = raw_or_z(money_level_chosen(choice_trial_idx));
                    case 4
                        choice_modVals(n_choiceMods,:) = raw_or_z(abs_money_level_chosen(choice_trial_idx));
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
                        choice_modVals(n_choiceMods,:) = raw_or_z(money_amount_unchosen(choice_trial_idx));
                    case 2
                        choice_modVals(n_choiceMods,:) = raw_or_z(abs_money_amount_unchosen(choice_trial_idx));
                    case 3
                        choice_modVals(n_choiceMods,:) = raw_or_z(money_level_unchosen(choice_trial_idx));
                    case 4
                        choice_modVals(n_choiceMods,:) = raw_or_z(abs_money_level_unchosen(choice_trial_idx));
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
                        choice_modVals(n_choiceMods,:) = raw_or_z(money_amount_varOption(choice_trial_idx));
                    case 2
                        choice_modVals(n_choiceMods,:) = raw_or_z(abs_money_amount_varOption(choice_trial_idx));
                    case 3
                        choice_modVals(n_choiceMods,:) = raw_or_z(money_level_varOption(choice_trial_idx));
                    case 4
                        choice_modVals(n_choiceMods,:) = raw_or_z(abs_money_level_varOption(choice_trial_idx));
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
                        choice_modVals(n_choiceMods,:) = raw_or_z(moneyChosen_min_moneyUnchosen_amount(choice_trial_idx));
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
                        choice_modVals(n_choiceMods,:) = raw_or_z(moneyChosen_min_moneyFixed_amount(choice_trial_idx));
                    case 2
                        choice_modVals(n_choiceMods,:) = raw_or_z(moneyChosen_min_moneyFixed_level(choice_trial_idx));
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
                        choice_modVals(n_choiceMods,:) = raw_or_z(money_amount_sum(choice_trial_idx));
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
                        choice_modVals(n_choiceMods,:) = raw_or_z(E_left(choice_trial_idx));
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
                        choice_modVals(n_choiceMods,:) = raw_or_z(E_right(choice_trial_idx));
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
                        choice_modVals(n_choiceMods,:) = raw_or_z(E_chosen(choice_trial_idx));
                    case 3
                        choice_modVals(n_choiceMods,:) = raw_or_z(E_chosen_bis(choice_trial_idx));
                    case 4
                        choice_modVals(n_choiceMods,:) = zscore(E_chosen(choice_trial_idx));
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
                        choice_modVals(n_choiceMods,:) = raw_or_z(E_unchosen(choice_trial_idx));
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
                        choice_modVals(n_choiceMods,:) = raw_or_z(E_varOption(choice_trial_idx));
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
                        choice_modVals(n_choiceMods,:) = raw_or_z(E_chosen_min_E_unchosen(choice_trial_idx));
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
                        choice_modVals(n_choiceMods,:) = raw_or_z(Ech_min_Efixed(choice_trial_idx));
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
                        choice_modVals(n_choiceMods,:) = raw_or_z(E_sum(choice_trial_idx));
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
                        choice_modVals(n_choiceMods,:) = raw_or_z(money_level_x_E_varOption(choice_trial_idx));
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
                        choice_modVals(n_choiceMods,:) = raw_or_z(money_level_x_E_chosen(choice_trial_idx));
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
                        choice_modVals(n_choiceMods,:) = raw_or_z(R_level_x_E_varOption(choice_trial_idx));
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
                        choice_modVals(n_choiceMods,:) = raw_or_z(R_level_x_E_chosen(choice_trial_idx));
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
                        choice_modVals(n_choiceMods,:) = raw_or_z(P_level_x_E_varOption(choice_trial_idx));
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
                        choice_modVals(n_choiceMods,:) = raw_or_z(P_level_x_E_chosen(choice_trial_idx));
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
                        choice_modVals(n_choiceMods,:) = raw_or_z(NV_ch_min_unch(choice_trial_idx));
                    case 2
                        choice_modNames{n_choiceMods} = 'p(chosen)';
                        choice_modVals(n_choiceMods,:) = raw_or_z(pChosen(choice_trial_idx));
                    case 3
                        choice_modNames{n_choiceMods} = 'NVch-NVunch';
                        choice_modVals(n_choiceMods,:) = raw_or_z(NV_ch_min_unch_with_bias(choice_trial_idx));
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
                        choice_modVals(n_choiceMods,:) = raw_or_z(NV_varOption(choice_trial_idx));
                    case 2
                        choice_modNames{n_choiceMods} = '|delta NV high E - low E|';
                        choice_modVals(n_choiceMods,:) = raw_or_z(abs(NV_varOption(choice_trial_idx)));
                    case 3
                        choice_modNames{n_choiceMods} = 'p(choice=hE)';
                        choice_modVals(n_choiceMods,:) = raw_or_z(pChoice_hE(choice_trial_idx));
                    case 4
                        choice_modNames{n_choiceMods} = 'delta NV high E - low E + bias';
                        choice_modVals(n_choiceMods,:) = raw_or_z(NV_varOption_plus_bias(choice_trial_idx));
                    case 5
                        choice_modNames{n_choiceMods} = '|delta NV high E - low E + bias|';
                        choice_modVals(n_choiceMods,:) = raw_or_z(abs(NV_varOption_plus_bias(choice_trial_idx)));
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
                        choice_modVals(n_choiceMods,:) = raw_or_z(NV_varOption(choice_trial_idx));
                    case 2
                        choice_modNames{n_choiceMods} = '|delta NV high E - low E|';
                        choice_modVals(n_choiceMods,:) = raw_or_z(abs(NV_varOption(choice_trial_idx)));
                    case 3
                        choice_modNames{n_choiceMods} = 'p(choice=hE)';
                        choice_modVals(n_choiceMods,:) = raw_or_z(pChoice_hE(choice_trial_idx));
                    case 4
                        choice_modNames{n_choiceMods} = 'delta NV high E - low E + bias';
                        choice_modVals(n_choiceMods,:) = raw_or_z(NV_varOption_plus_bias(choice_trial_idx));
                    case 5
                        choice_modNames{n_choiceMods} = '|delta NV high E - low E + bias|';
                        choice_modVals(n_choiceMods,:) = raw_or_z(abs(NV_varOption_plus_bias(choice_trial_idx)));
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
                            choice_modVals(n_choiceMods,:) = raw_or_z(AUC(choice_trial_idx));
                        case 2
                            choice_modVals(n_choiceMods,:) = raw_or_z(AUC_overshoot(choice_trial_idx));
                        case 3
                            choice_modVals(n_choiceMods,:) = raw_or_z(AUC_N(choice_trial_idx));
                        case 4
                            choice_modVals(n_choiceMods,:) = raw_or_z(AUC_overshoot_N(choice_trial_idx));
                        case 5
                            choice_modVals(n_choiceMods,:) = zscore(AUC(choice_trial_idx));
                        case 6
                            choice_modVals(n_choiceMods,:) = zscore(AUC_overshoot(choice_trial_idx));
                        case 7
                            choice_modVals(n_choiceMods,:) = zscore(AUC_N(choice_trial_idx));
                        case 8
                            choice_modVals(n_choiceMods,:) = zscore(AUC_overshoot_N(choice_trial_idx));
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
                            choice_modVals(n_choiceMods,:) = raw_or_z(fatigue(choice_trial_idx));
                        case 2
                            choice_modVals(n_choiceMods,:) = zscore(fatigue(choice_trial_idx));
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
                            choice_modVals(n_choiceMods,:) = raw_or_z(E_chosen(choice_trial_idx).*fatigue(choice_trial_idx));
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
                            choice_modVals(n_choiceMods,:) = raw_or_z(efficacy_with2first(choice_trial_idx));
                        case 2
                            choice_modVals(n_choiceMods,:) = raw_or_z(efficacy_pureNback(choice_trial_idx));
                        case 3
                            choice_modVals(n_choiceMods,:) = raw_or_z(efficacy_bis_with2first(choice_trial_idx));
                        case 4
                            choice_modVals(n_choiceMods,:) = raw_or_z(efficacy_bis_pureNback(choice_trial_idx));
                        case 5
                            choice_modVals(n_choiceMods,:) = zscore(efficacy_with2first(choice_trial_idx));
                        case 6
                            choice_modVals(n_choiceMods,:) = zscore(efficacy_pureNback(choice_trial_idx));
                        case 7
                            choice_modVals(n_choiceMods,:) = zscore(efficacy_bis_with2first(choice_trial_idx));
                        case 8
                            choice_modVals(n_choiceMods,:) = zscore(efficacy_bis_pureNback(choice_trial_idx));
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
                            choice_modVals(n_choiceMods,:) = raw_or_z(prevEfficacy_with2first(choice_trial_idx));
                        case 2
                            choice_modVals(n_choiceMods,:) = raw_or_z(prevEfficacy_pureNback(choice_trial_idx));
                        case 3
                            choice_modVals(n_choiceMods,:) = raw_or_z(prevEfficacy_bis_with2first(choice_trial_idx));
                        case 4
                            choice_modVals(n_choiceMods,:) = raw_or_z(prevEfficacy_bis_pureNback(choice_trial_idx));
                        case 5
                            choice_modVals(n_choiceMods,:) = zscore(prevEfficacy_with2first(choice_trial_idx));
                        case 6
                            choice_modVals(n_choiceMods,:) = zscore(prevEfficacy_pureNback(choice_trial_idx));
                        case 7
                            choice_modVals(n_choiceMods,:) = zscore(prevEfficacy_bis_with2first(choice_trial_idx));
                        case 8
                            choice_modVals(n_choiceMods,:) = zscore(prevEfficacy_bis_pureNback(choice_trial_idx));
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
                            choice_modVals(n_choiceMods,:) = raw_or_z(E_chosen(choice_trial_idx).*prevEfficacy_with2first(choice_trial_idx));
                        case 2
                            choice_modVals(n_choiceMods,:) = raw_or_z(E_chosen(choice_trial_idx).*prevEfficacy_pureNback(choice_trial_idx));
                        case 3
                            choice_modVals(n_choiceMods,:) = raw_or_z(E_chosen(choice_trial_idx).*prevEfficacy_bis_with2first(choice_trial_idx));
                        case 4
                            choice_modVals(n_choiceMods,:) = raw_or_z(E_chosen(choice_trial_idx).*prevEfficacy_bis_pureNback(choice_trial_idx));
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
                        choice_modVals(n_choiceMods,:) = raw_or_z(trialN(choice_trial_idx));
                    case 2
                        choice_modVals(n_choiceMods,:) = raw_or_z(trialN_dEch(choice_trial_idx));
                    case 3
                        choice_modVals(n_choiceMods,:) = raw_or_z(trialN_dEnonDef_min_Edef(choice_trial_idx));
                    case 4
                        choice_modVals(n_choiceMods,:) = raw_or_z(trialN_dEnonDef(choice_trial_idx));
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
                        choice_modVals(n_choiceMods,:) = confidence(choice_trial_idx);
                    case {2,3,4} % confidence inferred by the model => ok to zscore
                        choice_modVals(n_choiceMods,:) = raw_or_z(confidence(choice_trial_idx));
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
                        choice_modVals(n_choiceMods,:) = raw_or_z(choice_RT(choice_trial_idx));
                    case 2
                        choice_modVals(n_choiceMods,:) = zscore(choice_RT(choice_trial_idx));
                    otherwise
                        error('not ready yet');
                end
            end
            
            [matlabbatch] = First_level_loadEachCondition(matlabbatch, sub_idx, iRun, iCond,...
                ['choice_',RP_dispChoice_nm,'_',splitE_dispChoice_nm], modelChoiceOnset, modelChoiceDur,...
                n_choiceMods, choice_modNames, choice_modVals,...
                orth_vars, onsets_only_GLM);
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
            
            % extract trial index for the current loop
            switch RP_dispChosen_nm
                case 'RP'
                    RPfilter_dispChosen = true(1,length(double(RP_var_binary)));
                case 'R'
                    RPfilter_dispChosen = (RP_var_binary == 1);
                case 'P'
                    RPfilter_dispChosen = (RP_var_binary == 0);
            end
            switch splitE_dispChosen_nm
                case 'E'
                    Efilter_dispChosen = true(1,length(double(RP_var_binary)));
                case 'E1'
                    Efilter_dispChosen = (E_varOption == 1);
                case 'E2'
                    Efilter_dispChosen = (E_varOption == 2);
                case 'E3'
                    Efilter_dispChosen = (E_varOption == 3);
                case 'Ech0'
                    Efilter_dispChosen = (E_chosen == 0);
                case 'Ech1'
                    Efilter_dispChosen = (E_chosen == 1);
                case 'Ech2'
                    Efilter_dispChosen = (E_chosen == 2);
                case 'Ech3'
                    Efilter_dispChosen = (E_chosen == 3);
                case 'lEch'
                    Efilter_dispChosen = (choice_hE == 0);
                case 'hEch'
                    Efilter_dispChosen = (choice_hE == 1);
            end
            chosen_trial_idx = (RPfilter_dispChosen.*Efilter_dispChosen) == 1; % NEED to transform it into logical or will just focus on the first trial
            
            %% chosen onset
            iCond = iCond + 1;
            modelChosenOnset = dispChosenOnsets(chosen_trial_idx);
            % duration
            switch chosenModel
                case 'stick'
                    modelChosenDur = 0;
                case 'boxcar' % duration displaying the chosen option
                    modelChosenDur = dispChosenDur(chosen_trial_idx);
                case 'boxcar_bis' % duration going form display of chosen option
                    % until the end of the exertion of the effort
                    modelChosenDur = dispChosenDur(chosen_trial_idx) +...
                        preEffortCrossDur(chosen_trial_idx) +...
                        EperfDur(chosen_trial_idx);
                case 'boxcar_ter' % duration going form display of chosen option
                    % until the end of the effort preparation cross
                    % (entailing the whole effort preparation period)
                    modelChosenDur = dispChosenDur(chosen_trial_idx) +...
                        preEffortCrossDur(chosen_trial_idx);
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
                        chosen_modVals(n_chosenMods,:) = raw_or_z(choice_RT(chosen_trial_idx));
                    case 5
                        chosen_modVals(n_chosenMods,:) = zscore(choice_RT(chosen_trial_idx));
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
                        chosen_modVals(n_chosenMods,:) = raw_or_z(NV_ch_min_unch(chosen_trial_idx));
                    case 5
                        chosen_modNames{n_chosenMods} = 'p(chosen)';
                        chosen_modVals(n_chosenMods,:) = raw_or_z(pChosen(chosen_trial_idx));
                    case 6
                        chosen_modNames{n_chosenMods} = 'NVch-NVunch';
                        chosen_modVals(n_chosenMods,:) = raw_or_z(NV_ch_min_unch_with_bias(chosen_trial_idx));
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
                        chosen_modVals(n_chosenMods,:) = RP_var_binary(chosen_trial_idx);
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
                    chosen_modVals(n_chosenMods,:) = choice_hE(chosen_trial_idx); % binary variable => no zscore
                case 2
                    n_chosenMods = n_chosenMods + 1;
                    chosen_modNames{n_chosenMods} = 'choice = high effort';
                    chosen_modVals(n_chosenMods,:) = choice_hE_bis(chosen_trial_idx); % binary variable => no zscore
                otherwise
                    error('not ready yet');
            end
            
            % reward for the non-default option
            if chosenModel_R_varOption > 0
                n_chosenMods = n_chosenMods + 1;
                chosen_modNames{n_chosenMods} = 'R non-default';
                switch chosenModel_R_varOption
                    case 1 % R high effort option (amount)
                        chosen_modVals(n_chosenMods,:) = raw_or_z(R_amount_varOption(chosen_trial_idx));
                    case 2 % R high effort option (level)
                        chosen_modVals(n_chosenMods,:) = raw_or_z(R_level_varOption(chosen_trial_idx));
                    case 3 % z(R high effort option) (amount)
                        chosen_modVals(n_chosenMods,:) = raw_or_z(z_R_amount_varOption(chosen_trial_idx));
                    case 4 % z(R high effort option) (level)
                        chosen_modVals(n_chosenMods,:) = raw_or_z(z_R_level_varOption(chosen_trial_idx));
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
                        chosen_modVals(n_chosenMods,:) = raw_or_z(R_amount_chosen(chosen_trial_idx));
                    case 2 % R high effort option (level)
                        chosen_modVals(n_chosenMods,:) = raw_or_z(R_level_chosen(chosen_trial_idx));
                    case 3 % z(R high effort option) (amount)
                        chosen_modVals(n_chosenMods,:) = raw_or_z(z_R_amount_chosen(chosen_trial_idx));
                    case 4 % z(R high effort option) (level)
                        chosen_modVals(n_chosenMods,:) = raw_or_z(z_R_level_chosen(chosen_trial_idx));
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
                        chosen_modVals(n_chosenMods,:) = raw_or_z(P_amount_varOption(chosen_trial_idx));
                    case 2 % R high effort option (level)
                        chosen_modVals(n_chosenMods,:) = raw_or_z(P_level_varOption(chosen_trial_idx));
                    case 3 % z(R high effort option) (amount)
                        chosen_modVals(n_chosenMods,:) = raw_or_z(z_P_amount_varOption(chosen_trial_idx));
                    case 4 % z(R high effort option) (level)
                        chosen_modVals(n_chosenMods,:) = raw_or_z(z_P_level_varOption(chosen_trial_idx));
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
                        chosen_modVals(n_chosenMods,:) = raw_or_z(P_amount_chosen(chosen_trial_idx));
                    case 2 % R high effort option (level)
                        chosen_modVals(n_chosenMods,:) = raw_or_z(P_level_chosen(chosen_trial_idx));
                    case 3 % z(R high effort option) (amount)
                        chosen_modVals(n_chosenMods,:) = raw_or_z(z_P_amount_chosen(chosen_trial_idx));
                    case 4 % z(R high effort option) (level)
                        chosen_modVals(n_chosenMods,:) = raw_or_z(z_P_level_chosen(chosen_trial_idx));
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
                        chosen_modVals(n_chosenMods,:) = raw_or_z(money_amount_left(chosen_trial_idx));
                    case 2
                        chosen_modVals(n_chosenMods,:) = raw_or_z(abs(money_amount_left(chosen_trial_idx)));
                    case 3
                        chosen_modVals(n_chosenMods,:) = raw_or_z(money_level_left(chosen_trial_idx));
                    case 4
                        chosen_modVals(n_chosenMods,:) = raw_or_z(abs(money_level_chosen(chosen_trial_idx)));
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
                        chosen_modVals(n_chosenMods,:) = raw_or_z(money_amount_right(chosen_trial_idx));
                    case 2
                        chosen_modVals(n_chosenMods,:) = raw_or_z(abs(money_amount_right(chosen_trial_idx)));
                    case 3
                        chosen_modVals(n_chosenMods,:) = raw_or_z(money_level_right(chosen_trial_idx));
                    case 4
                        chosen_modVals(n_chosenMods,:) = raw_or_z(abs(money_level_right(chosen_trial_idx)));
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
                        chosen_modVals(n_chosenMods,:) = raw_or_z(money_amount_chosen(chosen_trial_idx));
                    case 2
                        chosen_modVals(n_chosenMods,:) = raw_or_z(abs_money_amount_chosen(chosen_trial_idx));
                    case 3
                        chosen_modVals(n_chosenMods,:) = raw_or_z(money_level_chosen(chosen_trial_idx));
                    case 4
                        chosen_modVals(n_chosenMods,:) = raw_or_z(abs_money_level_chosen(chosen_trial_idx));
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
                        chosen_modVals(n_chosenMods,:) = raw_or_z(money_amount_unchosen(chosen_trial_idx));
                    case 2
                        chosen_modVals(n_chosenMods,:) = raw_or_z(abs_money_amount_unchosen(chosen_trial_idx));
                    case 3
                        chosen_modVals(n_chosenMods,:) = raw_or_z(money_level_unchosen(chosen_trial_idx));
                    case 4
                        chosen_modVals(n_chosenMods,:) = raw_or_z(abs_money_level_unchosen(chosen_trial_idx));
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
                        chosen_modVals(n_chosenMods,:) = raw_or_z(money_amount_varOption(chosen_trial_idx));
                    case 2
                        chosen_modVals(n_chosenMods,:) = raw_or_z(abs_money_amount_varOption(chosen_trial_idx));
                    case 3
                        chosen_modVals(n_chosenMods,:) = raw_or_z(money_level_varOption(chosen_trial_idx));
                    case 4
                        chosen_modVals(n_chosenMods,:) = raw_or_z(abs_money_level_varOption(chosen_trial_idx));
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
                        chosen_modVals(n_chosenMods,:) = raw_or_z(moneyChosen_min_moneyUnchosen_amount(chosen_trial_idx));
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
                        chosen_modVals(n_chosenMods,:) = raw_or_z(moneyChosen_min_moneyFixed_amount(chosen_trial_idx));
                    case 2
                        chosen_modVals(n_chosenMods,:) = raw_or_z(moneyChosen_min_moneyFixed_level(chosen_trial_idx));
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
                        chosen_modVals(n_chosenMods,:) = raw_or_z(money_amount_sum(chosen_trial_idx));
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
                        chosen_modVals(n_chosenMods,:) = raw_or_z(E_left(chosen_trial_idx));
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
                        chosen_modVals(n_chosenMods,:) = raw_or_z(E_right(chosen_trial_idx));
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
                        chosen_modVals(n_chosenMods,:) = raw_or_z(E_chosen(chosen_trial_idx));
                    case 3
                        chosen_modVals(n_chosenMods,:) = raw_or_z(E_chosen_bis(chosen_trial_idx));
                    case 4
                        chosen_modVals(n_chosenMods,:) = zscore(E_chosen(chosen_trial_idx));
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
                        chosen_modVals(n_chosenMods,:) = raw_or_z(E_unchosen(chosen_trial_idx));
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
                        chosen_modVals(n_chosenMods,:) = E_varOption(chosen_trial_idx); % binary variable => no zscore
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
                        chosen_modVals(n_chosenMods,:) = raw_or_z(E_chosen_min_E_unchosen(chosen_trial_idx));
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
                        chosen_modVals(n_chosenMods,:) = raw_or_z(Ech_min_Efixed(chosen_trial_idx));
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
                        chosen_modVals(n_chosenMods,:) = raw_or_z(E_sum(chosen_trial_idx));
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
                        chosen_modVals(n_chosenMods,:) = raw_or_z(money_level_x_E_varOption(chosen_trial_idx));
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
                        chosen_modVals(n_chosenMods,:) = raw_or_z(money_level_x_E_chosen(chosen_trial_idx));
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
                        chosen_modVals(n_chosenMods,:) = raw_or_z(R_level_x_E_varOption(chosen_trial_idx));
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
                        chosen_modVals(n_chosenMods,:) = raw_or_z(R_level_x_E_chosen(chosen_trial_idx));
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
                        chosen_modVals(n_chosenMods,:) = raw_or_z(P_level_x_E_varOption(chosen_trial_idx));
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
                        chosen_modVals(n_chosenMods,:) = raw_or_z(P_level_x_E_chosen(chosen_trial_idx));
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
                        chosen_modVals(n_chosenMods,:) = raw_or_z(NV_ch_min_unch(chosen_trial_idx));
                    case 2
                        chosen_modNames{n_chosenMods} = 'p(chosen)';
                        chosen_modVals(n_chosenMods,:) = raw_or_z(pChosen(chosen_trial_idx));
                    case 3
                        chosen_modNames{n_chosenMods} = 'NVch-NVunch';
                        chosen_modVals(n_chosenMods,:) = raw_or_z(NV_ch_min_unch_with_bias(chosen_trial_idx));
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
                        chosen_modVals(n_chosenMods,:) = raw_or_z(NV_varOption(chosen_trial_idx));
                    case 2
                        chosen_modNames{n_chosenMods} = '|delta NV high E - low E|';
                        chosen_modVals(n_chosenMods,:) = raw_or_z(abs(NV_varOption(chosen_trial_idx)));
                    case 3
                        chosen_modNames{n_chosenMods} = 'p(choice=hE)';
                        chosen_modVals(n_chosenMods,:) = raw_or_z(abs(pChoice_hE(chosen_trial_idx)));
                    case 4
                        chosen_modNames{n_chosenMods} = 'delta NV high E - low E + bias';
                        chosen_modVals(n_chosenMods,:) = raw_or_z(NV_varOption_plus_bias(chosen_trial_idx));
                    case 5
                        chosen_modNames{n_chosenMods} = '|delta NV high E - low E + bias|';
                        chosen_modVals(n_chosenMods,:) = raw_or_z(abs(NV_varOption_plus_bias(chosen_trial_idx)));
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
                        chosen_modVals(n_chosenMods,:) = raw_or_z(NV_varOption(chosen_trial_idx));
                    case 2
                        chosen_modNames{n_chosenMods} = '|delta NV high E - low E|';
                        chosen_modVals(n_chosenMods,:) = raw_or_z(abs(NV_varOption(chosen_trial_idx)));
                    case 3
                        chosen_modNames{n_chosenMods} = 'p(choice=hE)';
                        chosen_modVals(n_chosenMods,:) = raw_or_z(abs(pChoice_hE(chosen_trial_idx)));
                    case 4
                        chosen_modNames{n_chosenMods} = 'delta NV high E - low E + bias';
                        chosen_modVals(n_chosenMods,:) = raw_or_z(NV_varOption_plus_bias(chosen_trial_idx));
                    case 5
                        chosen_modNames{n_chosenMods} = '|delta NV high E - low E + bias|';
                        chosen_modVals(n_chosenMods,:) = raw_or_z(abs(NV_varOption_plus_bias(chosen_trial_idx)));
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
                            chosen_modVals(n_chosenMods,:) = raw_or_z(AUC(chosen_trial_idx));
                        case 2
                            chosen_modVals(n_chosenMods,:) = raw_or_z(AUC_overshoot(chosen_trial_idx));
                        case 3
                            chosen_modVals(n_chosenMods,:) = raw_or_z(AUC_N(chosen_trial_idx));
                        case 4
                            chosen_modVals(n_chosenMods,:) = raw_or_z(AUC_overshoot_N(chosen_trial_idx));
                        case 5
                            chosen_modVals(n_chosenMods,:) = zscore(AUC(chosen_trial_idx));
                        case 6
                            chosen_modVals(n_chosenMods,:) = zscore(AUC_overshoot(chosen_trial_idx));
                        case 7
                            chosen_modVals(n_chosenMods,:) = zscore(AUC_N(chosen_trial_idx));
                        case 8
                            chosen_modVals(n_chosenMods,:) = zscore(AUC_overshoot_N(chosen_trial_idx));
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
                            chosen_modVals(n_chosenMods,:) = raw_or_z(fatigue(chosen_trial_idx));
                        case 2
                            chosen_modVals(n_chosenMods,:) = zscore(fatigue(chosen_trial_idx));
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
                            chosen_modVals(n_chosenMods,:) = raw_or_z(E_chosen(chosen_trial_idx).*fatigue(chosen_trial_idx));
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
                            chosen_modVals(n_chosenMods,:) = raw_or_z(efficacy_with2first(chosen_trial_idx));
                        case 2
                            chosen_modVals(n_chosenMods,:) = raw_or_z(efficacy_pureNback(chosen_trial_idx));
                        case 3
                            chosen_modVals(n_chosenMods,:) = raw_or_z(efficacy_bis_with2first(chosen_trial_idx));
                        case 4
                            chosen_modVals(n_chosenMods,:) = raw_or_z(efficacy_bis_pureNback(chosen_trial_idx));
                        case 5
                            chosen_modVals(n_chosenMods,:) = zscore(efficacy_with2first(chosen_trial_idx));
                        case 6
                            chosen_modVals(n_chosenMods,:) = zscore(efficacy_pureNback(chosen_trial_idx));
                        case 7
                            chosen_modVals(n_chosenMods,:) = zscore(efficacy_bis_with2first(chosen_trial_idx));
                        case 8
                            chosen_modVals(n_chosenMods,:) = zscore(efficacy_bis_pureNback(chosen_trial_idx));
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
                            chosen_modVals(n_chosenMods,:) = raw_or_z(prevEfficacy_with2first(chosen_trial_idx));
                        case 2
                            chosen_modVals(n_chosenMods,:) = raw_or_z(prevEfficacy_pureNback(chosen_trial_idx));
                        case 3
                            chosen_modVals(n_chosenMods,:) = raw_or_z(prevEfficacy_bis_with2first(chosen_trial_idx));
                        case 4
                            chosen_modVals(n_chosenMods,:) = raw_or_z(prevEfficacy_bis_pureNback(chosen_trial_idx));
                        case 5
                            chosen_modVals(n_chosenMods,:) = zscore(prevEfficacy_with2first(chosen_trial_idx));
                        case 6
                            chosen_modVals(n_chosenMods,:) = zscore(prevEfficacy_pureNback(chosen_trial_idx));
                        case 7
                            chosen_modVals(n_chosenMods,:) = zscore(prevEfficacy_bis_with2first(chosen_trial_idx));
                        case 8
                            chosen_modVals(n_chosenMods,:) = zscore(prevEfficacy_bis_pureNback(chosen_trial_idx));
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
                            chosen_modVals(n_chosenMods,:) = raw_or_z(Echosen(chosen_trial_idx).*prevEfficacy_with2first(chosen_trial_idx));
                        case 2
                            chosen_modVals(n_chosenMods,:) = raw_or_z(Echosen(chosen_trial_idx).*prevEfficacy_pureNback(chosen_trial_idx));
                        case 3
                            chosen_modVals(n_chosenMods,:) = raw_or_z(Echosen(chosen_trial_idx).*prevEfficacy_bis_with2first(chosen_trial_idx));
                        case 4
                            chosen_modVals(n_chosenMods,:) = raw_or_z(Echosen(chosen_trial_idx).*prevEfficacy_bis_pureNback(chosen_trial_idx));
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
                        chosen_modVals(n_chosenMods,:) = raw_or_z(trialN(chosen_trial_idx));
                    case 2
                        chosen_modVals(n_chosenMods,:) = raw_or_z(trialN_dEch(chosen_trial_idx));
                    case 3
                        chosen_modVals(n_chosenMods,:) = raw_or_z(trialN_dEnonDef_min_Edef(chosen_trial_idx));
                    case 4
                        chosen_modVals(n_chosenMods,:) = raw_or_z(trialN_dEnonDef(chosen_trial_idx));
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
                        chosen_modVals(n_chosenMods,:) = confidence(chosen_trial_idx); % binary variable => no zscore
                    case {2,3,4}
                        chosen_modVals(n_chosenMods,:) = raw_or_z(confidence(chosen_trial_idx)); % confidence inferred by the model
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
                        chosen_modVals(n_chosenMods,:) = raw_or_z(choice_RT(chosen_trial_idx));
                    case 2
                        chosen_modVals(n_chosenMods,:) = zscore(choice_RT(chosen_trial_idx));
                    otherwise
                        error('not ready yet');
                end
            end
            
            [matlabbatch] = First_level_loadEachCondition(matlabbatch, sub_idx, iRun, iCond,...
                ['dispChosen_',RP_dispChosen_nm,'_',splitE_dispChosen_nm], modelChosenOnset, modelChosenDur,...
                n_chosenMods, chosen_modNames, chosen_modVals,...
                orth_vars, onsets_only_GLM);
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
            
            % extract trial index for the current loop
            switch RP_preEcross_nm
                case 'RP'
                    RPfilter_preEcross = true(1,length(double(RP_var_binary)));
                case 'R'
                    RPfilter_preEcross = (RP_var_binary == 1);
                case 'P'
                    RPfilter_preEcross = (RP_var_binary == 0);
            end
            switch splitE_preEcross_nm
                case 'E'
                    Efilter_preEcross = true(1,length(double(RP_var_binary)));
                case 'E1'
                    Efilter_preEcross = (E_varOption == 1);
                case 'E2'
                    Efilter_preEcross = (E_varOption == 2);
                case 'E3'
                    Efilter_preEcross = (E_varOption == 3);
                case 'Ech0'
                    Efilter_preEcross = (E_chosen == 0);
                case 'Ech1'
                    Efilter_preEcross = (E_chosen == 1);
                case 'Ech2'
                    Efilter_preEcross = (E_chosen == 2);
                case 'Ech3'
                    Efilter_preEcross = (E_chosen == 3);
                case 'lEch'
                    Efilter_preEcross = (choice_hE == 0);
                case 'hEch'
                    Efilter_preEcross = (choice_hE == 1);
            end
            preEcross_trial_idx = (RPfilter_preEcross.*Efilter_preEcross) == 1; % NEED to transform it into logical or will just focus on the first trial
            
            %% pre-effort cross onset
            iCond = iCond + 1;
            modelpreEcrossOnset = preEffortCrossOnsets(preEcross_trial_idx);
            % duration
            switch preEffortCrossModel
                case 'stick'
                    modelPreEffortCrossdur = 0;
                case 'boxcar'
                    modelPreEffortCrossdur = preEffortCrossDur(preEcross_trial_idx);
                case 'boxcar_bis'
                    modelPreEffortCrossdur = preEffortCrossDur(preEcross_trial_idx) +...
                        EperfDur(preEcross_trial_idx);
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
                    preEcross_modVals(n_preEcrossMods,:) = choice_hE(preEcross_trial_idx); % binary variable => no zscore
                case 2
                    n_preEcrossMods = n_preEcrossMods + 1;
                    preEcross_modNames{n_preEcrossMods} = 'choice = high effort';
                    preEcross_modVals(n_preEcrossMods,:) = choice_hE_bis(preEcross_trial_idx); % binary variable => no zscore
                otherwise
                        error('not ready yet');
            end
            
            % money chosen
            if preEcrossModel_money_chosen > 0
                n_preEcrossMods = n_preEcrossMods + 1;
                preEcross_modNames{n_preEcrossMods} = 'money chosen';
                switch preEcrossModel_money_chosen
                    case 1
                        preEcross_modVals(n_preEcrossMods,:) = raw_or_z(money_amount_chosen(preEcross_trial_idx));
                    case 2
                        preEcross_modVals(n_preEcrossMods,:) = raw_or_z(abs_money_amount_chosen(preEcross_trial_idx));
                    case 3
                        preEcross_modVals(n_preEcrossMods,:) = raw_or_z(money_level_chosen(preEcross_trial_idx));
                    case 4
                        preEcross_modVals(n_preEcrossMods,:) = raw_or_z(abs_money_level_chosen(preEcross_trial_idx));
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
                        preEcross_modVals(n_preEcrossMods,:) = raw_or_z(E_chosen(preEcross_trial_idx));
                    case 3
                        preEcross_modVals(n_preEcrossMods,:) = raw_or_z(E_chosen_bis(preEcross_trial_idx));
                    case 4
                        preEcross_modVals(n_preEcrossMods,:) = zscore(E_chosen(preEcross_trial_idx));
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
                            preEcross_modVals(n_preEcrossMods,:) = raw_or_z(forcePeak(preEcross_trial_idx));
                        case 2
                            preEcross_modVals(n_preEcrossMods,:) = raw_or_z(forcePeak_N(preEcross_trial_idx));
                        case 3
                            preEcross_modVals(n_preEcrossMods,:) = zscore(forcePeak(preEcross_trial_idx));
                        case 4
                            preEcross_modVals(n_preEcrossMods,:) = zscore(forcePeak_N(preEcross_trial_idx));
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
                        preEcross_modVals(n_preEcrossMods,:) = raw_or_z(NV_ch_min_unch(preEcross_trial_idx));
                    case 2
                        preEcross_modNames{n_preEcrossMods} = 'p(chosen)';
                        preEcross_modVals(n_preEcrossMods,:) = raw_or_z(pChosen(preEcross_trial_idx));
                    case 3
                        preEcross_modNames{n_preEcrossMods} = 'NVch-NVunch';
                        preEcross_modVals(n_preEcrossMods,:) = raw_or_z(NV_ch_min_unch_with_bias(preEcross_trial_idx));
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
                        preEcross_modVals(n_preEcrossMods,:) = raw_or_z(NV_varOption(preEcross_trial_idx));
                    case 2
                        preEcross_modNames{n_preEcrossMods} = '|delta NV high E - low E|';
                        preEcross_modVals(n_preEcrossMods,:) = raw_or_z(abs(NV_varOption(preEcross_trial_idx)));
                    case 3
                        preEcross_modNames{n_preEcrossMods} = 'p(choice=hE)';
                        preEcross_modVals(n_preEcrossMods,:) = raw_or_z(abs(pChoice_hE(preEcross_trial_idx)));
                    case 4
                        preEcross_modNames{n_preEcrossMods} = 'delta NV high E - low E + bias';
                        preEcross_modVals(n_preEcrossMods,:) = raw_or_z(NV_varOption_plus_bias(preEcross_trial_idx));
                    case 5
                        preEcross_modNames{n_preEcrossMods} = '|delta NV high E - low E + bias|';
                        preEcross_modVals(n_preEcrossMods,:) = raw_or_z(abs(NV_varOption_plus_bias(preEcross_trial_idx)));
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
                        preEcross_modVals(n_preEcrossMods,:) = raw_or_z(NV_varOption(preEcross_trial_idx));
                    case 2
                        preEcross_modNames{n_preEcrossMods} = '|delta NV high E - low E|';
                        preEcross_modVals(n_preEcrossMods,:) = raw_or_z(abs(NV_varOption(preEcross_trial_idx)));
                    case 3
                        preEcross_modNames{n_preEcrossMods} = 'p(choice=hE)';
                        preEcross_modVals(n_preEcrossMods,:) = raw_or_z(abs(pChoice_hE(preEcross_trial_idx)));
                    case 4
                        preEcross_modNames{n_preEcrossMods} = 'delta NV high E - low E + bias';
                        preEcross_modVals(n_preEcrossMods,:) = raw_or_z(NV_varOption_plus_bias(preEcross_trial_idx));
                    case 5
                        preEcross_modNames{n_preEcrossMods} = '|delta NV high E - low E + bias|';
                        preEcross_modVals(n_preEcrossMods,:) = raw_or_z(abs(NV_varOption_plus_bias(preEcross_trial_idx)));
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
                        preEcross_modVals(n_preEcrossMods,:) = raw_or_z(latency(preEcross_trial_idx));
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
                        preEcross_modVals(n_preEcrossMods,:) = raw_or_z(trialN(preEcross_trial_idx));
                    case 2
                        preEcross_modVals(n_preEcrossMods,:) = raw_or_z(trialN_dEch(preEcross_trial_idx));
                    case 3
                        preEcross_modVals(n_preEcrossMods,:) = raw_or_z(trialN_dEnonDef_min_Edef(preEcross_trial_idx));
                    case 4
                        preEcross_modVals(n_preEcrossMods,:) = raw_or_z(trialN_dEnonDef(preEcross_trial_idx));
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
                        preEcross_modVals(n_preEcrossMods,:) = confidence(preEcross_trial_idx);
                    case {2,3,4} % confidence inferred by the model => ok to zscore
                        preEcross_modVals(n_preEcrossMods,:) = raw_or_z(confidence(preEcross_trial_idx));
                    otherwise
                        error('not ready yet');
                end
            end
            
            [matlabbatch] = First_level_loadEachCondition(matlabbatch, sub_idx, iRun, iCond,...
                ['preEffort fixation cross',RP_preEcross_nm,'_',splitE_preEcross_nm], modelpreEcrossOnset, modelPreEffortCrossdur,...
                n_preEcrossMods, preEcross_modNames, preEcross_modVals,...
                orth_vars, onsets_only_GLM);
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
            
            % extract trial index for the current loop
            switch RP_Eperf_nm
                case 'RP'
                    RPfilter_Eperf = true(1,length(double(RP_var_binary)));
                case 'R'
                    RPfilter_Eperf = (RP_var_binary == 1);
                case 'P'
                    RPfilter_Eperf = (RP_var_binary == 0);
            end
            switch splitE_Eperf_nm
                case 'E'
                    Efilter_Eperf = true(1,length(double(RP_var_binary)));
                case 'E1'
                    Efilter_Eperf = (E_varOption == 1);
                case 'E2'
                    Efilter_Eperf = (E_varOption == 2);
                case 'E3'
                    Efilter_Eperf = (E_varOption == 3);
                case 'Ech0'
                    Efilter_Eperf = (E_chosen == 0);
                case 'Ech1'
                    Efilter_Eperf = (E_chosen == 1);
                case 'Ech2'
                    Efilter_Eperf = (E_chosen == 2);
                case 'Ech3'
                    Efilter_Eperf = (E_chosen == 3);
                case 'lEch'
                    Efilter_Eperf = (choice_hE == 0);
                case 'hEch'
                    Efilter_Eperf = (choice_hE == 1);
            end
            switch onsets_only_GLM
                case 0
                    Eperf_trial_idx = (RPfilter_Eperf.*Efilter_Eperf.*perfOkTrials) == 1; % NEED to transform it into logical or will just focus on the first trial
                case 1 % keep same number of trials as in choices and feedback to avoid problems in the analysis
                    Eperf_trial_idx = (RPfilter_Eperf.*Efilter_Eperf) == 1;
            end
            
            %% Effort performance onset
            iCond = iCond + 1;
            modelEperfOnset = EperfOnsets(Eperf_trial_idx);
            % duration
            switch EperfModel
                case 'stick'
                    modelEperfDur = 0;
                case 'boxcar'
                    modelEperfDur = EperfDur(Eperf_trial_idx);
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
                    Eperf_modVals(n_EperfMods,:) = choice_hE(Eperf_trial_idx); % binary variable => no zscore
                case 2
                    n_EperfMods = n_EperfMods + 1;
                    Eperf_modNames{n_EperfMods} = 'choice = high effort';
                    Eperf_modVals(n_EperfMods,:) = choice_hE_bis(Eperf_trial_idx); % binary variable => no zscore
                otherwise
                    error('not ready yet');
            end
            
            % money chosen
            if EperfModel_money_chosen > 0
                n_EperfMods = n_EperfMods + 1;
                Eperf_modNames{n_EperfMods} = 'money chosen';
                switch EperfModel_money_chosen
                    case 1
                        Eperf_modVals(n_EperfMods,:) = raw_or_z(money_amount_chosen(Eperf_trial_idx));
                    case 2
                        Eperf_modVals(n_EperfMods,:) = raw_or_z(abs_money_amount_chosen(Eperf_trial_idx));
                    case 3
                        Eperf_modVals(n_EperfMods,:) = raw_or_z(money_level_chosen(Eperf_trial_idx));
                    case 4
                        Eperf_modVals(n_EperfMods,:) = raw_or_z(abs_money_level_chosen(Eperf_trial_idx));
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
                        Eperf_modVals(n_EperfMods,:) = raw_or_z(E_chosen(Eperf_trial_idx));
                    case 3
                        Eperf_modVals(n_EperfMods,:) = raw_or_z(E_chosen_bis(Eperf_trial_idx));
                    case 4
                        Eperf_modVals(n_EperfMods,:) = zscore(E_chosen(Eperf_trial_idx));
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
                            Eperf_modVals(n_EperfMods,:) = raw_or_z(forcePeak(Eperf_trial_idx));
                        case 2
                            Eperf_modVals(n_EperfMods,:) = raw_or_z(forcePeak_N(Eperf_trial_idx));
                        case 3
                            Eperf_modVals(n_EperfMods,:) = zscore(forcePeak(Eperf_trial_idx));
                        case 4
                            Eperf_modVals(n_EperfMods,:) = zscore(forcePeak_N(Eperf_trial_idx));
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
                            Eperf_modVals(n_EperfMods,:) = raw_or_z(AUC(Eperf_trial_idx));
                        case 2
                            Eperf_modVals(n_EperfMods,:) = raw_or_z(AUC_overshoot(Eperf_trial_idx));
                        case 3
                            Eperf_modVals(n_EperfMods,:) = raw_or_z(AUC_N(Eperf_trial_idx));
                        case 4
                            Eperf_modVals(n_EperfMods,:) = raw_or_z(AUC_overshoot_N(Eperf_trial_idx));
                        case 5
                            Eperf_modVals(n_EperfMods,:) = zscore(AUC(Eperf_trial_idx));
                        case 6
                            Eperf_modVals(n_EperfMods,:) = zscore(AUC_overshoot(Eperf_trial_idx));
                        case 7
                            Eperf_modVals(n_EperfMods,:) = zscore(AUC_N(Eperf_trial_idx));
                        case 8
                            Eperf_modVals(n_EperfMods,:) = zscore(AUC_overshoot_N(Eperf_trial_idx));
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
                            Eperf_modVals(n_EperfMods,:) = raw_or_z(efficacy_with2first(Eperf_trial_idx));
                        case 2
                            Eperf_modVals(n_EperfMods,:) = raw_or_z(efficacy_pureNback(Eperf_trial_idx));
                        case 3
                            Eperf_modVals(n_EperfMods,:) = raw_or_z(efficacy_bis_with2first(Eperf_trial_idx));
                        case 4
                            Eperf_modVals(n_EperfMods,:) = raw_or_z(efficacy_bis_pureNback(Eperf_trial_idx));
                        case 5
                            Eperf_modVals(n_EperfMods,:) = zscore(efficacy_with2first(Eperf_trial_idx));
                        case 6
                            Eperf_modVals(n_EperfMods,:) = zscore(efficacy_pureNback(Eperf_trial_idx));
                        case 7
                            Eperf_modVals(n_EperfMods,:) = zscore(efficacy_bis_with2first(Eperf_trial_idx));
                        case 8
                            Eperf_modVals(n_EperfMods,:) = zscore(efficacy_bis_pureNback(Eperf_trial_idx));
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
                            Eperf_modVals(n_EperfMods,:) = raw_or_z(RT_avg(Eperf_trial_idx));
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
                            Eperf_modVals(n_EperfMods,:) = raw_or_z(n_correct(Eperf_trial_idx));
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
                            Eperf_modVals(n_EperfMods,:) = raw_or_z(n_errors(Eperf_trial_idx));
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
                        Eperf_modVals(n_EperfMods,:) = raw_or_z(NV_ch_min_unch(Eperf_trial_idx));
                    case 2
                        Eperf_modNames{n_EperfMods} = 'p(chosen)';
                        Eperf_modVals(n_EperfMods,:) = raw_or_z(pChosen(Eperf_trial_idx));
                    case 3
                        Eperf_modNames{n_EperfMods} = 'NVch-NVunch';
                        Eperf_modVals(n_EperfMods,:) = raw_or_z(NV_ch_min_unch_with_bias(Eperf_trial_idx));
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
                        Eperf_modVals(n_EperfMods,:) = raw_or_z(NV_varOption(Eperf_trial_idx));
                    case 2
                        Eperf_modNames{n_EperfMods} = '|delta NV high E - low E|';
                        Eperf_modVals(n_EperfMods,:) = raw_or_z(abs(NV_varOption(Eperf_trial_idx)));
                    case 3
                        Eperf_modNames{n_EperfMods} = 'p(choice=hE)';
                        Eperf_modVals(n_EperfMods,:) = raw_or_z(abs(pChoice_hE(Eperf_trial_idx)));
                    case 4
                        Eperf_modNames{n_EperfMods} = 'delta NV high E - low E + bias';
                        Eperf_modVals(n_EperfMods,:) = raw_or_z(NV_varOption_plus_bias(Eperf_trial_idx));
                    case 5
                        Eperf_modNames{n_EperfMods} = '|delta NV high E - low E + bias|';
                        Eperf_modVals(n_EperfMods,:) = raw_or_z(abs(NV_varOption_plus_bias(Eperf_trial_idx)));
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
                        Eperf_modVals(n_EperfMods,:) = raw_or_z(NV_varOption(Eperf_trial_idx));
                    case 2
                        Eperf_modNames{n_EperfMods} = '|delta NV high E - low E|';
                        Eperf_modVals(n_EperfMods,:) = raw_or_z(abs(NV_varOption(Eperf_trial_idx)));
                    case 3
                        Eperf_modNames{n_EperfMods} = 'p(choice=hE)';
                        Eperf_modVals(n_EperfMods,:) = raw_or_z(abs(pChoice_hE(Eperf_trial_idx)));
                    case 4
                        Eperf_modNames{n_EperfMods} = 'delta NV high E - low E + bias';
                        Eperf_modVals(n_EperfMods,:) = raw_or_z(NV_varOption_plus_bias(Eperf_trial_idx));
                    case 5
                        Eperf_modNames{n_EperfMods} = '|delta NV high E - low E + bias|';
                        Eperf_modVals(n_EperfMods,:) = raw_or_z(abs(NV_varOption_plus_bias(Eperf_trial_idx)));
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
                        Eperf_modVals(n_EperfMods,:) = raw_or_z(latency(Eperf_trial_idx));
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
                                Eperf_modVals(n_EperfMods,:) = raw_or_z(fatigue(Eperf_trial_idx));
                            case 2
                                Eperf_modVals(n_EperfMods,:) = zscore(fatigue(Eperf_trial_idx));
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
                                Eperf_modVals(n_EperfMods,:) = raw_or_z(E_chosen(Eperf_trial_idx).*fatigue(Eperf_trial_idx));
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
                                Eperf_modVals(n_EperfMods,:) = raw_or_z(prevEfficacy_with2first(Eperf_trial_idx));
                            case 2
                                Eperf_modVals(n_EperfMods,:) = raw_or_z(prevEfficacy_pureNback(Eperf_trial_idx));
                            case 3
                                Eperf_modVals(n_EperfMods,:) = raw_or_z(prevEfficacy_bis_with2first(Eperf_trial_idx));
                            case 4
                                Eperf_modVals(n_EperfMods,:) = raw_or_z(prevEfficacy_bis_pureNback(Eperf_trial_idx));
                            case 5
                                Eperf_modVals(n_EperfMods,:) = zscore(prevEfficacy_with2first(Eperf_trial_idx));
                            case 6
                                Eperf_modVals(n_EperfMods,:) = zscore(prevEfficacy_pureNback(Eperf_trial_idx));
                            case 7
                                Eperf_modVals(n_EperfMods,:) = zscore(prevEfficacy_bis_with2first(Eperf_trial_idx));
                            case 8
                                Eperf_modVals(n_EperfMods,:) = zscore(prevEfficacy_bis_pureNback(Eperf_trial_idx));
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
                                Eperf_modVals(n_EperfMods,:) = raw_or_z(E_chosen(Eperf_trial_idx).*prevEfficacy_with2first(Eperf_trial_idx));
                            case 2
                                Eperf_modVals(n_EperfMods,:) = raw_or_z(E_chosen(Eperf_trial_idx).*prevEfficacy_pureNback(Eperf_trial_idx));
                            case 3
                                Eperf_modVals(n_EperfMods,:) = raw_or_z(E_chosen(Eperf_trial_idx).*prevEfficacy_bis_with2first(Eperf_trial_idx));
                            case 4
                                Eperf_modVals(n_EperfMods,:) = raw_or_z(E_chosen(Eperf_trial_idx).*prevEfficacy_bis_pureNback(Eperf_trial_idx));
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
                        Eperf_modVals(n_EperfMods,:) = raw_or_z(trialN(Eperf_trial_idx));
                    case 2
                        Eperf_modVals(n_EperfMods,:) = raw_or_z(trialN_dEch(Eperf_trial_idx));
                    case 3
                        Eperf_modVals(n_EperfMods,:) = raw_or_z(trialN_dEnonDef_min_Edef(Eperf_trial_idx));
                    case 4
                        Eperf_modVals(n_EperfMods,:) = raw_or_z(trialN_dEnonDef(Eperf_trial_idx));
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
                        Eperf_modVals(n_EperfMods,:) = confidence(Eperf_trial_idx);
                    case {2,3,4} % confidence inferred by the model => ok to zscore
                        Eperf_modVals(n_EperfMods,:) = raw_or_z(confidence(Eperf_trial_idx));
                    otherwise
                        error('not ready yet');
                end
            end
            
            [matlabbatch] = First_level_loadEachCondition(matlabbatch, sub_idx, iRun, iCond,...
                ['Eperf_',RP_Eperf_nm,'_',splitE_Eperf_nm], modelEperfOnset, modelEperfDur,...
                n_EperfMods, Eperf_modNames, Eperf_modVals,...
                orth_vars, onsets_only_GLM);
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
            
            % extract trial index for the current loop
            switch RP_fbk_nm
                case 'RP'
                    RPfilter_fbk = true(1,length(double(RP_var_binary)));
                case 'R'
                    RPfilter_fbk = (RP_var_binary == 1);
                case 'P'
                    RPfilter_fbk = (RP_var_binary == 0);
            end
            switch splitE_fbk_nm
                case 'E'
                    Efilter_fbk = true(1,length(double(RP_var_binary)));
                case 'E1'
                    Efilter_fbk = (E_varOption == 1);
                case 'E2'
                    Efilter_fbk = (E_varOption == 2);
                case 'E3'
                    Efilter_fbk = (E_varOption == 3);
                case 'Ech0'
                    Efilter_fbk = (E_chosen == 0);
                case 'Ech1'
                    Efilter_fbk = (E_chosen == 1);
                case 'Ech2'
                    Efilter_fbk = (E_chosen == 2);
                case 'Ech3'
                    Efilter_fbk = (E_chosen == 3);
                case 'lEch'
                    Efilter_fbk = (choice_hE == 0);
                case 'hEch'
                    Efilter_fbk = (choice_hE == 1);
            end
            fbk_trial_idx = (RPfilter_fbk.*Efilter_fbk) == 1; % NEED to transform it into logical or will just focus on the first trial
            
            %% feedback onset
            iCond = iCond + 1;
            modelFbkOnset = fbkOnsets(fbk_trial_idx);
            % duration
            switch fbkModel
                case 'stick'
                    modelFbkDur = 0;
                case 'boxcar'
                    modelFbkDur = fbkDur(fbk_trial_idx);
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
                            fbk_modVals(n_fbkMods,:) = win_vs_loss_fbk;
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
                    fbk_modVals(n_fbkMods,:) = choice_hE(fbk_trial_idx); % binary variable => no zscore
                case 2
                    n_fbkMods = n_fbkMods + 1;
                    fbk_modNames{n_fbkMods} = 'choice = high effort';
                    fbk_modVals(n_fbkMods,:) = choice_hE_bis(fbk_trial_idx); % binary variable => no zscore
                otherwise
                    error('not ready yet');
            end
            
            % money chosen
            if fbkModel_moneyObtained > 0
                n_fbkMods = n_fbkMods + 1;
                fbk_modNames{n_fbkMods} = 'money obtained';
                switch fbkModel_moneyObtained
                    case 1 % money amount
                        fbk_modVals(n_fbkMods,:) = raw_or_z(money_amount_obtained(fbk_trial_idx));
                    case 2 % |money amount|
                        fbk_modVals(n_fbkMods,:) = raw_or_z(abs(money_amount_obtained(fbk_trial_idx)));
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
                        fbk_modVals(n_fbkMods,:) = raw_or_z(trialN(fbk_trial_idx));
                    case 2
                        fbk_modVals(n_fbkMods,:) = raw_or_z(trialN_dEch(fbk_trial_idx));
                    case 3
                        fbk_modVals(n_fbkMods,:) = raw_or_z(trialN_dEnonDef_min_Edef(fbk_trial_idx));
                    case 4
                        fbk_modVals(n_fbkMods,:) = raw_or_z(trialN_dEnonDef(fbk_trial_idx));
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
                        fbk_modVals(n_choiceMods,:) = confidence(fbk_trial_idx);
                    case {2,3,4} % confidence inferred by the model => ok to zscore
                        fbk_modVals(n_choiceMods,:) = raw_or_z(confidence(fbk_trial_idx));
                    otherwise
                        error('not ready yet');
                end
            end
            
            [matlabbatch] = First_level_loadEachCondition(matlabbatch, sub_idx, iRun, iCond,...
                ['fbk_',RP_fbk_nm,'_',splitE_fbk_nm], modelFbkOnset, modelFbkDur,...
                n_fbkMods, fbk_modNames, fbk_modVals,...
                orth_vars, onsets_only_GLM);
        end % E condition
    end % RP
end % model fbk period

end % function