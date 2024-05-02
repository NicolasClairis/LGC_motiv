function[] = computational_mdl(mdl_n)
% [] = computational_mdl(mdl_n)
% computational_mdl will extract the corresponding data for the model
% defined in input. All the data will be saved in the following folder:
% C:\Users\clairis\Desktop\GitHub\LGC_motiv\LGC_Motiv_results\study1\bayesian_modeling
% under the labels "bayesian_deltaNV_data.mat", "bayesian_pChoice_data.mat"
% and "behavioral_prm.mat".
%
% INPUTS
% mdl_n: model number


%% check inputs
if ~exist('mdl_n','var') || isempty(mdl_n)
    mdl_n = input('What is the model number?');
end

%% define subjects
study_nm = 'study1';
[subject_id, NS] = LGCM_subject_selection(study_nm, 'behavior');

%% define working directories
root = 'E:';

%% main parameters
n_trialsPerSession = 54;
nRuns = 4;
n_totalTrials = n_trialsPerSession.*nRuns;
% apply multisession to each block. 
is_multisession = true;

% load model parameters
[mdl_prm] = computational_mdl_prm(mdl_n);
binary_answers = mdl_prm.binary_answers;

%% define bayesian evolution and observation functions
% no evolution function with current models
n_F_prm = 0;
n_hiddenStates = 0;
f_function = [];

% observation function
switch mdl_n
    case 1
        g_function = @g_observation_mdl1;
    case 2
        g_function = @g_observation_mdl2;
    case 3
        g_function = @g_observation_mdl3;
    case 4
        g_function = @g_observation_mdl4;
    case 5
        g_function = @g_observation_mdl5;
    case 6
        g_function = @g_observation_mdl6;
    otherwise
        error(['Define g function for model ',num2str(mdl_nm)]);
end

%% initialize variables of interest
[] = deal(NaN());
%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    subBehaviorFolder = [fullfile(root,['CID',sub_nm],'behavior'),filesep];
    
    % initialize variables of interest for this subject
    [deltaR, deltaP, deltaE,...
        Fp, currEff, prevEff,...
        choice_binary, choice_with_conf,...
        Ep_or_Em_trials] = deal(NaN(n_totalTrials,1));
    
    % extract relevant runs
    [runs, n_runs] = runs_definition(study_nm, sub_nm, condition);
    % extract corresponding data for each run and pool all in one big
    % vector
    for iR = 1:n_runs
        % run-related informations
        jR = runs.runsToKeep(iR);
        run_nm = num2str(jR);
        run_task_nm = runs.tasks{iR};
        switch run_task_nm
            case 'Ep'
                task_fullName = 'physical';
            case 'Em'
                task_fullName = 'mental';
        end
        run_trial_idx = (1:n_trialsPerSession) + n_trialsPerSession.*(jR - 1);
        
        % extract variables of interest
        [deltaR(run_trial_idx)] = extract_deltaR_money(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        [deltaP(run_trial_idx)] = extract_deltaP_money(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        [deltaE(run_trial_idx)] = extract_hE_level(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        switch run_task_nm
            case 'Ep'
                [sumPrevAUC_N, sumPrevAUC] = extract_physical_fatigue(subBehaviorFolder, sub_nm, run_nm, task_fullName);
                Fp(run_trial_idx) = sumPrevAUC_N;
                currEff(run_trial_idx) = 0;
                prevEff(run_trial_idx) = 0;
            case 'Em'
                Fp(run_trial_idx) = 0;
                [~, ~, ~, ~,...
                    ~,...
                    ~,...
                    efficacy_with2first,...
                    efficacy_pureNback,...
                    efficacy_bis_with2first,...
                    efficacy_bis_pureNback, ~,...
                    efficacy_ter_with2first,...
                    efficacy_ter_pureNback] = extract_mental_perf(subBehaviorFolder, sub_nm, run_nm);
                currEff(run_trial_idx) = efficacy_ter_with2first;
                [prevEfficacy_with2first,...
                    prevEfficacy_pureNback,...
                    prevEfficacy_bis_with2first,...
                    prevEfficacy_bis_pureNback,...
                    prevEfficacy_ter_with2first,...
                    prevEfficacy_ter_pureNback] = extract_mental_previous_efficacy(subBehaviorFolder, sub_nm, run_nm, task_fullName);
                prevEff(run_trial_idx) = prevEfficacy_ter_with2first;
        end
        
        % choice output
        [choice_binary(run_trial_idx)] = extract_choice_hE(subBehaviorFolder,...
            sub_nm, run_nm, task_fullName);
        [choice_with_conf(run_trial_idx)] = extract_choice_hE_bis(subBehaviorFolder,...
            sub_nm, run_nm, task_fullName);
        
        % define task
        switch run_task_nm
            case 'Ep'
                Ep_or_Em_trials(run_trial_idx) = 1;
            case 'Em'
                Ep_or_Em_trials(run_trial_idx) = 0;
        end % task
    end % run loop
    
    %% adjustement of the variables to facilitate VBA modeling
    % monetary amounts in cents
    deltaR = deltaR.*100;
    deltaP = deltaP.*100;
    % physical fatigue in smaller range
    Fp = Fp./1000;
    
    %% pool everybody in var
    var = [deltaR, deltaP, deltaE, Ep_or_Em_trials, Fp, currEff, prevEff];
    
    %% define if output is binary or 4-levels
    switch binary_answers
        case false
            all_choices = choice_with_conf;
        case true
            all_choices = choice_binary;
    end
    
    %% compute model for this subject
    [posterior,out] = VBA_NLStateSpaceModel(all_choices, var, f_function, g_function, options.dim, options);
end % subject loop

%% save results

end % function