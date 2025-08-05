function[prm, mdl_quality, subject_id, NS, choices_raw, choices_pred] = computational_mdl(mdl_n)
% [prm, mdl_quality, subject_id, NS, choices_raw, choices_pred] = computational_mdl(mdl_n)
% computational_mdl will extract the corresponding data for the model
% defined in input. All the data will be saved in the following folder:
% C:\Users\clairis\Desktop\GitHub\LGC_motiv\LGC_Motiv_results\study1\bayesian_modeling
% under the labels "bayesian_deltaNV_data.mat", "bayesian_pChoice_data.mat"
% and "behavioral_prm.mat".
%
% INPUTS
% mdl_n: model number (see computational_mdl_prm.m for actual model parameters corresponding to this model number)
%
% OUTPUTS
% prm: structure with individual parameters for the model selected
%
% mdl_quality: structure with variables reflecting model quality (R², BIC,
% AIC, free energy, etc.)
%
% subject_id: list of subjects included
%
% NS: total number of subjects
%
% choices_raw: nTrials*NS matrix with raw choices across sessions and
% subjects
%
% choices_pred: nTrials*NS matrix with raw choices across sessions and
% subjects


%% check inputs
if exist('mdl_n','var')
    [mdl_n, mdl_n_nm] = which_bayesian_mdl_n(mdl_n);
else
    [mdl_n, mdl_n_nm] = which_bayesian_mdl_n;
end

%% define subjects
study_nm = 'study1';
condition1 = 'behavior_noSatTaskSub_noSatRun_lenient'; % by default, include all behavioral sessions except those where behavior was saturated
[subject_id, NS] = LGCM_subject_selection(study_nm, condition1);
condition2 = 'behavior_noSatTaskSub'; % this will allow to extract the information regarding the inputs for all trials, even though runs will be excluded from the analysis

%% define working directories
root = 'F:';
root_saveFolder = fullfile('C:','Users','Nicolas Clairis','Documents'); % pc-specific pathway to update if script launched elsewhere
saveFolder = fullfile(root_saveFolder,'GitHub','LGC_motiv','LGC_Motiv_results','study1','bayesian_modeling');

%% main parameters
nTrialsPerRun = 54;
nRuns = 4;
nTotalTrials = nTrialsPerRun.*nRuns;
% apply multisession (ie VBA will allow different noise values in the 
% estimation of the parameters for each block considering that some blocks 
% may be more trustful than others for the estimation of the parameters)
is_multisession = true;

% load model parameters
[mdl_prm] = computational_mdl_prm(mdl_n);
binary_answers = mdl_prm.binary_answers;

%% define bayesian evolution and observation functions
% no evolution function with current models
n_F_prm = 0;
n_hiddenStates = 0;
f_evol_function = [];
g_obs_function = @g_observation_mdl;

%% initialize variables of interest
% model parameters
G_parameters = mdl_prm.G_prm_names;
n_G_prm = length(G_parameters);
for iP = 1:n_G_prm
    G_prm_nm = G_parameters{iP};
    [old_prm.(G_prm_nm), prm.(G_prm_nm)] = deal(NaN(1,NS));
end

% model quality
[mdl_quality.R2,...
    mdl_quality.AIC,...
    mdl_quality.BIC,...
    mdl_quality.free_energy] = deal(NaN(1, NS));

% output
[choices_raw, choices_pred,...
    dV_pred] = deal(NaN(nTotalTrials, NS));

%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    sub_nm_bis = ['CID',sub_nm];
    subBehaviorFolder = [fullfile(root,study_nm,['CID',sub_nm],'behavior'),filesep];
    
    % initialize variables of interest for this subject
    options = struct;
    % input variables
    [deltaR, deltaP, deltaE,...
        Fp, currEff, prevEff,...
        choice_binary, choice_with_conf,...
        Ep_or_Em_trials] = deal(NaN(1,nTotalTrials));
    % indication of trials to exclude from the analysis (when no choice was made)
    options.isYout = ones(1,nTotalTrials);
    
    % extract data individually
    [choices_raw_perSub.(sub_nm_bis), choices_pred_perSub.(sub_nm_bis),...
        dV_pred_perSub.(sub_nm_bis)] = deal(NaN(nTotalTrials,1));
    
    % extract relevant runs
    [runs_ok, n_runs_ok] = runs_definition(study_nm, sub_nm, condition1); % allows to identify runs saturated to remove from analysis
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
        
        % keep trials when run not saturated (isYout = 0),
        % otherwise ignore the run from the analysis (isYout = 1) to 
        % improve the model estimation
        if ismember(jR, runs_ok.runsToKeep)
            options.isYout(run_trial_idx) = 0;
        end

        % extract variables of interest
        [deltaR(run_trial_idx)] = extract_deltaR_money(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        [deltaP(run_trial_idx)] = extract_deltaP_money(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        [deltaE(run_trial_idx)] = extract_hE_level(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        switch run_task_nm
            case 'Ep'
                [sumPrevAUC_N, sumPrevAUC] = extract_physical_fatigue(subBehaviorFolder, sub_nm, run_nm, task_fullName);
                Fp(run_trial_idx) = sumPrevAUC;
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
    
    % remove trials where no choice was made
    options.isYout(isnan(choice_binary)) = 1;
    
    %% adjustement of the variables to facilitate VBA modeling
    % monetary amounts in cents
    deltaR = deltaR.*100;
    deltaP = deltaP.*100;
    % physical fatigue in smaller range
    Fp = Fp./1000;
    
    %% pool everybody in var
    vars = [deltaR; deltaP; deltaE; Ep_or_Em_trials; Fp; currEff; prevEff];
    var_names = {'dR','dP','dE','EpEm','Fp','currEff','prevEff'};
    % control no NaNs remaining as that would make VBA crash
    if sum(sum(isnan(vars))) > 0
        switch sub_nm
            case '040' % run 3 and run 4 not executed => replace with zeros to avoid VBA crash, but trials will be ignored through isYout
                vars(:,109:216) = 0;
            otherwise
                % identify variable and trial which are problematic and report it
                iVar = 1;
                while sum(isnan(vars(:,iVar))) == 0 && iVar < size(vars,2)
                    iVar = iVar + 1;
                    error(['var contains a NaN in trial ',num2str(find(isnan(vars(:,iVar)))),' for the variable ',var_names{iVar}]);
                end

                % only in case where var is equal to NaN for saturated
                % runs:
                % % replace all NaN trials by zero for all variables to avoid
                % % VBA crash, trials should be ignored through isYout
                % % function anyway
                % var(isnan(var)) = 0;
        end % subject filter
    end % NaN filter
    
    %% define if output is binary or 4-levels
    switch binary_answers
        case false % 4-level choices: 0/0.25/0.75/1
            all_choices = choice_with_conf;
        case true % binary choices: 0/1 = low/high E chosen
            all_choices = choice_binary;
    end
    choices_raw(:,iS) = all_choices; % store data for output
    choices_raw_perSub.(sub_nm_bis)(:) = choices_raw(:,iS);
    % since NaN makes the toolbox crash, and that we ignore NaN choices using isYout in the
    % options, we redefine the no choice with a -1 idx, to take it into account later
    all_choices(isnan(all_choices)) = -1;
    
    %% extract main informations regarding the model
    % unless you want to follow the evolution of the fit, avoid displaying
    % VBA fitting procedure because it significantly slows down the whole
    % procedure
    options.DisplayWin = 0; % display window during inversion (general recommendation: avoid as it will take much time but good to try the first time to verify the data fits well)
    options.verbose = 0; % display text during inversion (1) or not (0)
    options.GnFigs = 0; % flag to display (=1) or not (=0) the Gauss-Newton inner loops display figures {0} by default
    
    % define nature of the output (binary/non-binary)
    switch binary_answers
        case false
            options.sources.type = 0; % output is not binary
        case true
            options.sources.type = 1; % output is binary
    end
    % since you define options.sources.type, you also have to define
    % options.sources.out or VBA_check will crash but can be left empty if
    % you don't define options.sources.type as VBA will automatically fill
    % these fields together.
    options.sources.out = 1:size(all_choices,1); % vector varying with the number of output variables (in principle = 1 in this design)
    
    % priors on the parameters
    % no hidden states, nor F function => no priors
    options.priors.muX0 = [];
    options.priors.SigmaX0 = [];
    options.priors.muTheta = [];
    options.priors.SigmaTheta = [];
    % all parameters are in G function with mean centered at 0 and sigma of
    % 100 to have enough space
    options.priors.muPhi = zeros(n_G_prm, 1);
    options.priors.SigmaPhi = eye(n_G_prm).*100;
    % summary for dim
    dim = struct('n',n_hiddenStates,'n_t',nTotalTrials,...
        'n_theta',n_F_prm,'n_phi',n_G_prm);
    options.dim = dim;
    % number of hidden states (n)
    % nb of trials (n_t),
    % nb of theta parameters in f function (n_theta)
    % and nb of Phi parameters in G function (n_phi)
    
    % Arthur specific modifications:
    % data to force the number of iteration of the model. seems ignored
    options.MinIter = 3; % minimum number of VB iterations {1} by default
%     options.MaxIter = 10; % maximum number of VB iterations {32} by
%     default, Arthur constrained it to 10 which allows to go faster during testing, but I think it's better to
%     leave the default option as recommended by Jules.
    options.TolFun  = 1e-7 ; % minimum absolute increase of the free energy {2e-2} by default
    
    % multisession
    switch is_multisession
        case true
            options.multisession.split = repmat(nTrialsPerRun,1,nRuns); % split in 4 equal sessions
            % fix the parameters to be equal across sessions
            options.multisession.fixed.phi = 1:n_G_prm;
        otherwise
            error(['is_multisession = ',num2str(is_multisession),' not ready yet.']);
    end
    
    % information about model parameters (positivity constraint, etc.)
    options.inG.mdl_prm = mdl_prm;
    options.inG.var_names = var_names;
    
    %% compute model for this subject
    [posterior, out] = VBA_NLStateSpaceModel(all_choices, vars, f_evol_function, g_obs_function, dim, options);
    
    %% extract variables of interest
    % model prediction for choices
    choices_pred(:,iS) = out.suffStat.gx;
    choices_pred_perSub.(sub_nm_bis)(:) = choices_pred(:,iS);

    % parameters
    for iP = 1:n_G_prm
        G_prm_nm = G_parameters{iP};
        switch mdl_prm.pos.(G_prm_nm)
            case false % no constraint on the parameter
                prm.(G_prm_nm)(iS) = posterior.muPhi(iP);
            case true % positivity constraint => transform parameter accordingly
                [prm.(G_prm_nm)(iS)] = fn_for_posterior(posterior.muPhi(iP), posterior.SigmaPhi(iP,iP), 'pos2'); % 'pos2' refers to the log(1+exp(X)) transformation
                old_prm.(G_prm_nm)(iS) = log(1 + exp(posterior.muPhi(iP))); % wrong but this is what we did initially => need to control if it's ok
        end
    end

    % extract difference in net value
    switch mdl_n
        case {1,2} % no bias
            dV_pred(:,iS) = -log((1-out.suffStat.gx)./out.suffStat.gx);
        case {3,4,5,6} % choice bias
            dV_pred(:,iS) = -prm.kBias(iS) - log((1-out.suffStat.gx)./out.suffStat.gx);
    end
    dV_pred_perSub.(sub_nm_bis)(:) = dV_pred(:,iS);

    % extract the data separately for each run
    for iR = 1:nRuns
        run_name = ['run',num2str(iR)];
        run_trials_idx = (1:nTrialsPerRun) + nTrialsPerRun.*(iR - 1);
        choices_raw_perSub_perRun.(sub_nm_bis).(run_name) = choices_raw_perSub.(sub_nm_bis)(run_trials_idx);
        choices_pred_perSub_perRun.(sub_nm_bis).(run_name) = choices_pred_perSub.(sub_nm_bis)(run_trials_idx);
        dV_pred_perSub_perRun.(sub_nm_bis).(run_name) = dV_pred_perSub.(sub_nm_bis)(run_trials_idx);
    end % run loop
    
    % model quality
    mdl_quality.R2(iS) = out.fit.R2;
    mdl_quality.AIC(iS) = out.fit.AIC;
    mdl_quality.BIC(iS) = out.fit.BIC;
    mdl_quality.free_energy(iS) = out.F;
    mdl_quality.RMSE(iS) = sqrt(mean( ((choices_raw(:,iS) - choices_pred(:,iS)).^2),'omitnan'));
    mdl_quality.MAE(iS) = mean( abs(choices_raw(:,iS) - choices_pred(:,iS)),'omitnan');
    
    % display information about where the script is
    disp(['Subject ',num2str(iS),'/',num2str(NS)]);
end % subject loop

%% save results
save([saveFolder,filesep,'bayesian_model_',mdl_n_nm,'_results.mat'],...
    'prm', 'mdl_quality', 'subject_id', 'NS',...
    'choices_raw', 'choices_pred','dV_pred',...
    'choices_raw_perSub','choices_pred_perSub','dV_pred_perSub',...
    'choices_raw_perSub_perRun','choices_pred_perSub_perRun','dV_pred_perSub_perRun',...
    'old_prm');
end % function