function[pval, a_path, b_path, c_path, c_prime_path, subject_id,...
    input_f_input_bin, ROI_f_input_bin, output_f_input_bin,...
    ROI_f_ROI_bin, output_f_ROI_bin,...
    input_prm_nm, ROI_short_nm, output_prm_nm] = mediation_behavior_fMRI_behavior(study_nm, condition, subject_id)
%[pval, a_path, b_path, c_path, c_prime_path, subject_id,...
%     input_f_input_bin, ROI_f_input_bin, output_f_input_bin,...
%     ROI_f_ROI_bin, output_f_ROI_bin,...
%     input_prm_nm, ROI_short_nm, output_prm_nm] = mediation_behavior_fMRI_behavior(study_nm, condition, subject_id)
% script aiming at checking whether fMRI activity deconvoluted and 
% extracted trial/trial mediates the impact of task variables to behavior.
% In particular, we are interested in looking whether the degree of
% uncertainty related to the task variables trial/trial has an impact on
% RT mediated by dmPFC.
%
% INPUTS
% study_nm: study name ('study1'/'study2')
%
% condition: condition to use for the subjects (important to know which
% runs to include or not) (will be asked if left empty)
%
% subject_id: list of subjects to include (will be asked if left empty)
%
% OUTPUTS
% pval: p.value for each path of the mediation
%
% a_path: beta for path going from behavior to fMRI for each subject
% 
% b_path: beta for path going from fMRI to behavior for each subject
%
% c_path: beta for direct path from behavior input to behavior output for
% each subject
%
% c_prime_path: beta for direct path from behavior input to behavior output
% for each subject after taking into account the mediator
%
% subject_id: list of subjects included
%
% input_f_input_bin: bin of behavioral input as of function of itself
% 
% ROI_f_input_bin: bin of ROI activity as of function of behavioral input
% 
% output_f_input_bin: bin of behavioral output as of function of behavioral input
% 
% ROI_f_ROI_bin: bin of ROI activity as of function of itself
%     
% output_f_ROI_bin: bin of behavioral output as of function of ROI activity
%
% input_prm_nm: input parameter name
%
% ROI_short_nm: ROI short name
%
% output_prm_nm: output parameter name

%% study by default
if ~exist('study_nm','var') || isempty(study_nm)
    study_nm = 'study1';
end

%% working directories
computerRoot = LGCM_root_paths;
studyBehaviorFolder = [computerRoot, filesep, study_nm, filesep];

%% selection of participants
if ~exist('condition','var') || isempty(condition)
    condition = subject_condition;
end
if ~exist('subject_id','var') || isempty(subject_id)
    [subject_id, NS] = LGCM_subject_selection(study_nm, condition);
else
    NS = length(subject_id);
end

%% extract ROI activity for all subjects
[ROI_trial_b_trial] = extract_ROI_betas_onsets_only(computerRoot,...
    study_nm, subject_id, condition);
% define which ROI, and which time period is of interest to you
% define ROI
ROI_names = fieldnames(ROI_trial_b_trial);
ROI_subList = ROI_trial_b_trial.subject_id;
ROI_names(strcmp(ROI_names,'subject_id')) = [];
if length(ROI_names) > 1
    which_ROI = listdlg('PromptString','Which ROI?','ListString',ROI_names);
    ROI_nm = ROI_names{which_ROI};
else
    ROI_nm = ROI_names;
end
ROI_short_nm = inputdlg('ROI short name?');
ROI_short_nm = ROI_short_nm{1};
% define task
task_names = {'Ep','Em','EpEmPool'};
which_task = listdlg('PromptString','Which task?','ListString',task_names);
task_to_look = task_names{which_task};
% define time period
timePeriods = fieldnames(ROI_trial_b_trial.(ROI_nm{1}).Ep.run1);
which_timePeriod = listdlg('PromptString','Which time phase of the trial?',...
    'listString',timePeriods);
timePeriod_nm = timePeriods{which_timePeriod};

%% select parameters of interest
potential_input_prm = {'NV_hE','NV_ch','uncertainty',...
    'E_level','money_level','deltaMoney_level'};
potential_output_prm = {'RT','uncertainty_rtg','choice_hE'};

which_input = listdlg('PromptString','please select input parameter',...
    'ListString',potential_input_prm);
input_prm_nm = potential_input_prm{which_input};
which_output = listdlg('PromptString','please select output parameter',...
    'ListString',potential_output_prm);
output_prm_nm = potential_output_prm{which_output};

%% if a parameter is based on modeling, you need to decide which model to use
% bayesian across tasks or simple model each task separately?
if ismember(input_prm_nm,{'NV_hE','NV_ch','uncertainty'})
    needModeling = 1;
    [mdlType, mdlN] = behavioral_model_selection;
    % warning until this gets fixed
    if strcmp(mdlType,'bayesian')
        error('uncertainty/NV extraction from bayesian model not ready yet.');
    end
else
    needModeling = 0;
end

%% extract the parameters of interest
nRuns = 4;
nTrialsPerRun = 54;
nTrials = nTrialsPerRun*nRuns;
% input/mediator/output
[input_prm, ROI_mediator, output_prm] = deal(NaN(nTrials, NS));
% bin variables
nBins = 3;
[input_f_input_bin, ROI_f_input_bin, output_f_input_bin,...
    ROI_f_ROI_bin, output_f_ROI_bin] = deal(NaN(nBins, NS));
% paths
[a_path, b_path, c_path, c_prime_path] = deal(NaN(1,NS));

for iS = 1:NS
    sub_nm = subject_id{iS};
    subBehaviorFolder = [studyBehaviorFolder, 'CID',sub_nm, filesep, 'behavior',filesep];
    
    % extract confidence based on the model
    if needModeling == 1
        if strcmp(mdlType,'simple')
            [~, dataInferred] = logitfit_choices(computerRoot, study_nm, sub_nm,...
                0, 'levels', 6, 6);
        end
    end
    
    % extract runs
    [runsStruct] = runs_definition(study_nm, sub_nm, 'behavior');
    okRuns = runsStruct.runsToKeep;
    taskNames = runsStruct.tasks;
    jRun = 0;
    for iRun = 1:length(okRuns)
        kRun = okRuns(iRun);
        run_nm = num2str(kRun);
        jRun = jRun + 1;
        task_nm_tmp = taskNames{jRun};
        runTrials_idx = (1:nTrialsPerRun) + nTrialsPerRun*(kRun-1);
        switch task_nm_tmp
            case 'Em'
                task_fullName = 'mental';
            case 'Ep'
                task_fullName = 'physical';
        end
        % define which task session it is
        switch kRun
            case {1,2}
                run_nm_bis = ['run',num2str(1)];
            case {3,4}
                run_nm_bis = ['run',num2str(2)];
        end
        
        %% load the data
        behaviorStruct_tmp = load([subBehaviorFolder,...
            'CID',sub_nm,'_session',run_nm,'_',task_fullName,...
            '_task.mat']);
        choiceOptions_tmp = behaviorStruct_tmp.choice_opt;
        switch task_nm_tmp
            case 'Em'
                choiceAndPerf_tmp = behaviorStruct_tmp.mentalE_perf;
            case 'Ep'
                choiceAndPerf_tmp = behaviorStruct_tmp.physicalPerf;
        end
        
        %% default side
        defaultSide_tmp = choiceOptions_tmp.default_LR;
        %% extract R or P
        RP_var_tmp = strcmp(choiceOptions_tmp.R_or_P,'R');
        
        %% effort level
        E_highE_tmp = (choiceOptions_tmp.E.left).*(defaultSide_tmp == 1) +...
            (choiceOptions_tmp.E.right).*(defaultSide_tmp == -1);
        
        %% high-effort money amount
        money_hE_tmp = ((choiceOptions_tmp.monetary_amount.left).*(defaultSide_tmp == 1) +...
            (choiceOptions_tmp.monetary_amount.right).*(defaultSide_tmp == -1)).*((RP_var_tmp == 1) - (RP_var_tmp == 0));
        money_lE_tmp = ((choiceOptions_tmp.monetary_amount.left).*(defaultSide_tmp == -1) +...
            (choiceOptions_tmp.monetary_amount.right).*(defaultSide_tmp == 1)).*((RP_var_tmp == 1) - (RP_var_tmp == 0));
        
        %% delta between high and low effort options
        deltaMoney_tmp = money_hE_tmp - money_lE_tmp;
        
        %% RT
        onsets_tmp = behaviorStruct_tmp.onsets;
        switch task_nm_tmp
            case 'Ep'
                onsets_tmp = behaviorStruct_tmp.physicalPerf.onsets;
                choice_LR_tmp = behaviorStruct_tmp.physicalPerf.choice;
            case 'Em'
                onsets_tmp = behaviorStruct_tmp.mentalE_perf.onsets;
                choice_LR_tmp = behaviorStruct_tmp.mentalE_perf.choice;
        end
        RT_tmp = onsets_tmp.choice - onsets_tmp.dispChoiceOptions;
        
        %% confidence rating
        uncertaintyRtg_tmp = NaN(1,length(choice_LR_tmp));
        uncertaintyRtg_tmp(abs(choice_LR_tmp) == 2) = 0; % high confidence = low uncertainty
        uncertaintyRtg_tmp(abs(choice_LR_tmp) == 1) = 1; % low confidence = high uncertainty
        
        %% confidence inferred by the model
        if needModeling == 1 
            switch mdlType
                case 'simple'
                    NV_hE_tmp = dataInferred.NV_varOption.(task_nm_tmp).(['mdl_',mdlN]).(run_nm_bis);
                    NV_ch_tmp = dataInferred.NV_chosen.(task_nm_tmp).(['mdl_',mdlN]).(run_nm_bis);
                    uncertainty_tmp = - dataInferred.confidenceFitted.(['mdl_',mdlN]).(run_nm_bis); % revert sign to transform confidence into uncertainty
                otherwise
                    error('not ready yet');
            end
        end
        
        %% choice
        choice_LR_tmp = choiceAndPerf_tmp.choice;
        % remove confidence info from choice:
        choice_LR_tmp(choice_LR_tmp == 2) = 1;
        choice_LR_tmp(choice_LR_tmp == -2) = -1;
        % extract high effort choice
        choice_highE_tmp = NaN(1,length(choice_LR_tmp));
        choice_highE_tmp(choice_LR_tmp == -defaultSide_tmp) = 1;
        choice_highE_tmp(choice_LR_tmp == defaultSide_tmp) = 0;
        
        %% extract input behavioral variable
        switch input_prm_nm
            case 'NV_hE'
                input_prm(runTrials_idx, iS) = NV_hE_tmp;
            case 'NV_ch'
                input_prm(runTrials_idx, iS) = NV_ch_tmp;
            case 'uncertainty'
                input_prm(runTrials_idx, iS) = uncertainty_tmp;
            case 'E_level'
                input_prm(runTrials_idx, iS) = E_highE_tmp;
            case 'money_level'
                input_prm(runTrials_idx, iS) = money_hE_tmp;
            case 'deltaMoney_level'
                input_prm(runTrials_idx, iS) = deltaMoney_tmp;
            otherwise
                error(['input = ',input_prm_nm,' not ready yet']);
        end
        
        %% extract fMRI ROI mediator
        if strcmp(task_to_look,'EpEmPool') ||...
                (strcmp(task_to_look, task_nm))
            ROI_mediator(runTrials_idx, iS) = ROI_trial_b_trial.(ROI_nm{1}).(task_nm_tmp).(run_nm_bis).(timePeriod_nm)(:, iS);
        end
        %% extract output behavioral variable
        switch output_prm_nm
            case 'RT'
                output_prm(runTrials_idx, iS) = RT_tmp;
            case 'uncertainty_rtg'
                output_prm(runTrials_idx, iS) = uncertaintyRtg_tmp;
            case 'choice_hE'
                output_prm(runTrials_idx, iS) = choice_highE_tmp;
            otherwise
                error(['output = ',output_prm_nm,' not ready yet']);
        end
    end % run loop
    
    %% do the mediation between input, fMRI and output per subject
    dispIndivData = 0;
    [a_path(iS), b_path(iS), c_path(iS), c_prime_path(iS)] = mediation(input_prm(:,iS), ROI_mediator(:,iS), output_prm(:,iS),...
        input_prm_nm, ROI_short_nm, output_prm_nm, dispIndivData);
    
    %% do the bins per subject
    [ROI_f_input_bin(:,iS), input_f_input_bin(:,iS)] = do_bin2(...
        ROI_mediator(:,iS), input_prm(:,iS), nBins, 0);
    output_f_input_bin(:,iS) = do_bin2(...
        output_prm(:,iS), input_prm(:,iS), nBins, 0);
    [output_f_ROI_bin(:,iS), ROI_f_ROI_bin(:,iS)] = do_bin2(...
        output_prm(:,iS), ROI_mediator(:,iS), nBins, 0);
end % subject loop

%% average
[m_input_f_input_bin,...
    sem_input_f_input_bin,...
    sd_input_f_input_bin] = mean_sem_sd(input_f_input_bin, 2);
[m_ROI_f_input_bin,...
    sem_ROI_f_input_bin,...
    sd_ROI_f_input_bin] = mean_sem_sd(ROI_f_input_bin, 2);
[m_output_f_input_bin,...
    sem_output_f_input_bin,...
    sd_output_f_input_bin] = mean_sem_sd(output_f_input_bin, 2);
[m_ROI_f_ROI_bin,...
    sem_ROI_f_ROI_bin,...
    sd_ROI_f_ROI_bin] = mean_sem_sd(ROI_f_ROI_bin, 2);
[m_output_f_ROI_bin,...
    sem_output_f_ROI_bin,...
    sd_output_f_ROI_bin] = mean_sem_sd(output_f_ROI_bin, 2);

%% test the mediation betas between input, fMRI and output across participants
[~,pval.a] = ttest(a_path);
[~,pval.b] = ttest(b_path);
[~,pval.c] = ttest(c_path);
[~,pval.c_prime] = ttest(c_prime_path);

%% figure showing different bins
pSize = 50;
lWidth = 3;
% start figure
fig;

% ROI = f(input)
subplot(1,3,1);
hdl = errorbar(m_input_f_input_bin,...
    m_ROI_f_input_bin,...
    sem_ROI_f_input_bin);
hdl.LineWidth = lWidth;
xlabel(input_prm_nm);
ylabel(ROI_short_nm);
legend_size(pSize);

% output = f(ROI)
subplot(1,3,2);
hdl = errorbar(m_ROI_f_ROI_bin,...
    m_output_f_ROI_bin,...
    sem_output_f_ROI_bin);
hdl.LineWidth = lWidth;
legend_size(pSize);
xlabel(ROI_short_nm);
ylabel(output_prm_nm);

% output = f(input)
subplot(1,3,3);
hdl = errorbar(m_output_f_input_bin,...
    m_output_f_input_bin,...
    sem_output_f_input_bin);
hdl.LineWidth = lWidth;
legend_size(pSize);
xlabel(input_prm_nm);
ylabel(output_prm_nm);
end % function