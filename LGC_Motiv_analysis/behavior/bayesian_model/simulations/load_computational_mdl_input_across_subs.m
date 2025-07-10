function[vars, mean_var, var_names, AUC_perEch, currEff_perEch, mean_AUC_perEch, mean_currEff_perEch] = load_computational_mdl_input_across_subs()
% [vars, mean_var, var_names, AUC_perEch, currEff_perEch, mean_AUC_perEch, mean_currEff_perEch] = load_computational_mdl_input_across_subs()
% load_computational_mdl_input_across_subs will go through all subjects of
% the experiments defined in input and will extract their data in order to
% build a big matrix with all the variables and to average the inputs 
% across subjects
%
% OUTPUTS
% vars: (7 input variables)*(216 trials)*(NS) matrix including the input
% variables for the model across all these subjects
%
% mean_var: (7 input variables)*(216 trials) matrix averaged across
% subjects of the var matrix
%
% var_names: list of the names corresponding to each of the 7 input
% variables
%
% AUC_perEch: (n_Ech_levels*NS) matrix with integral of the force exerted
% on each physical effort trial depending on the effort level and subject
%
% currEff_perEch: (n_Ech_levels*NS) matrix with efficiency on each mental
% effort trial depending on the effort level and subject
%
% mean_AUC_perEch & mean_currEff_perEch: same as two previous variables
% but averaged across subjects

%% define subjects to consider (same as the ones included in computational_mdl)
study_nm = 'study1';
condition1 = 'behavior_noSatTaskSub_noSatRun_lenient'; % by default, include all behavioral sessions except those where behavior was saturated
[subject_id, NS] = LGCM_subject_selection(study_nm, condition1);
condition2 = 'behavior_noSatTaskSub'; % this will allow to extract the information regarding the inputs for all trials, even though runs will be excluded from the analysis

%% define working directories
root = 'E:';

%% general parameters
nTrialsPerRun = 54;
nRuns = 4;
nTotalTrials = nTrialsPerRun.*nRuns;
var_names = {'dR','dP','dE','EpEm','Fp','currEff','prevEff'};
n_vars = length(var_names);
vars = NaN(n_vars, nTotalTrials, NS);

Ech_levels = 0:3;
n_Ech = length(Ech_levels);
[AUC_perEch, currEff_perEch] = deal(NaN(n_Ech, NS));

%% loop through subjects to load the data
for iS = 1:NS
    sub_nm = subject_id{iS};
    subBehaviorFolder = [fullfile(root,study_nm,['CID',sub_nm],'behavior'),filesep];
    
    % input variables
    [deltaR, deltaP, deltaE,...
        AUC, Fp, currEff, prevEff,...
        Ep_or_Em_trials,...
        choice_hE] = deal(NaN(1,nTotalTrials));
    ok_trials = true(1,nTotalTrials);
    
    % extract relevant runs
    [runs_ok] = runs_definition(study_nm, sub_nm, condition1); % allows to identify runs saturated to remove from analysis
    [allRuns, n_allRuns] = runs_definition(study_nm, sub_nm, condition2); % allows to extract all task inputs to still compute choice prediction
    % extract corresponding data for each run and pool all in one big
    % vector
    for iR = 1:n_allRuns
        % run-related informations
        jR = allRuns.runsToKeep(iR);
        run_nm = num2str(jR);
        run_task_nm = allRuns.tasks{iR};
        switch run_task_nm
            case 'Ep'
                task_fullName = 'physical';
            case 'Em'
                task_fullName = 'mental';
        end
        run_trial_idx = (1:nTrialsPerRun) + nTrialsPerRun.*(jR - 1);

        % extract variables of interest
        [deltaR(run_trial_idx)] = extract_deltaR_money(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        [deltaP(run_trial_idx)] = extract_deltaP_money(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        [deltaE(run_trial_idx)] = extract_hE_level(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        [choice_hE(run_trial_idx)] = extract_choice_hE(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        switch run_task_nm
            case 'Ep'
                [sumPrevAUC_N, sumPrevAUC] = extract_physical_fatigue(subBehaviorFolder, sub_nm, run_nm, task_fullName);
                Fp(run_trial_idx) = sumPrevAUC;
                [~,AUC_tmp] = extract_grip_force(subBehaviorFolder, sub_nm, run_nm);
                AUC(run_trial_idx) = AUC_tmp.allTrials;
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
                currEff(run_trial_idx) = efficacy_ter_with2first.allTrials;
%                 % to reproduce Arthur's mistake:
%                 currEff(run_trial_idx(1)) = 0;
%                 currEff(run_trial_idx(2:end)) = efficacy_ter_with2first.allTrials(2:end);
                [prevEfficacy_with2first,...
                    prevEfficacy_pureNback,...
                    prevEfficacy_bis_with2first,...
                    prevEfficacy_bis_pureNback,...
                    prevEfficacy_ter_with2first,...
                    prevEfficacy_ter_pureNback] = extract_mental_previous_efficacy(subBehaviorFolder, sub_nm, run_nm, task_fullName);
                prevEff(run_trial_idx) = prevEfficacy_ter_with2first;
        end
        
        % define task
        switch run_task_nm
            case 'Ep'
                Ep_or_Em_trials(run_trial_idx) = 1;
            case 'Em'
                Ep_or_Em_trials(run_trial_idx) = 0;
        end % task
        
        %% filter bad trials or sessions
        % remove trials where no choice was made
        [choice_tmp] = extract_choice_hE(subBehaviorFolder,...
            sub_nm, run_nm, task_fullName);
        ok_trials(run_trial_idx(isnan(choice_tmp))) = false;
        % keep these trials if run not saturated, otherwise ignore it from
        % the analysis
        if ~ismember(jR, runs_ok.runsToKeep)
            ok_trials(run_trial_idx) = false;
        end
    end % run loop
    
    %% variable range adjustements
    % monetary amounts in cents
    deltaR = deltaR.*100;
    deltaP = deltaP.*100;
    % physical fatigue in smaller range
    Fp = Fp./1000;
    
    %% pool everybody in var
    vars(:,:,iS) = [deltaR; deltaP; deltaE; Ep_or_Em_trials; Fp; currEff; prevEff];
    % remove bad trials
    bad_trials = ok_trials == false;
    vars(:,bad_trials,iS) = NaN;
    
    %% extract AUC and Efficiency per effort level
    for iEch = Ech_levels
        jEch = iEch + 1;
        switch iEch
            case 0 % low E chosen
                Ech_trial_idx = choice_hE == 0;
            otherwise % high E chosen
                Ech_trial_idx = deltaE == iEch;
        end
        % physical force exerted
        Ep_trials = (Ep_or_Em_trials == 1).*(Ech_trial_idx) == 1;
        AUC_perEch(jEch,iS) = mean(AUC(Ep_trials),2,'omitnan');
        % mental efficiency
        Em_trials = (Ep_or_Em_trials == 0).*(Ech_trial_idx) == 1;
        currEff_perEch(jEch,iS) = mean(currEff(Em_trials),2,'omitnan');
    end
end % subject loop

%% average the data across subjects
mean_var = mean(vars,3,'omitnan');
% replace Ep_or_Em_trials to fix it (alternating between Ep and Em across
% sessions) instead of weird average
mean_var(strcmp(var_names,'EpEm'),:) = [ones(1,nTrialsPerRun), zeros(1,nTrialsPerRun), ones(1,nTrialsPerRun), zeros(1,nTrialsPerRun)];

% average also AUC and currEff
mean_AUC_perEch = mean(AUC_perEch,2,'omitnan');
mean_currEff_perEch = mean(currEff_perEch,2,'omitnan');

end % function