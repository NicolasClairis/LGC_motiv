%% look whether within each input parameter bin, the activity of the selected ROI 
% determines which choices should be made.
%

%% number of bins
nBins = 6;

%% study by default
if ~exist('study_nm','var') || isempty(study_nm)
    study_nm = 'study1';
end

%% working directories
computerRoot = LGCM_root_paths;
studyBehaviorFolder = [computerRoot, filesep, study_nm, filesep];

%% selection of participants
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% extract ROI activity for all subjects
[ROI_trial_b_trial] = extract_ROI_betas_onsets_only(computerRoot,...
    study_nm, subject_id, condition);
% define which ROI, and which time period is of interest to you
% define ROI
ROI_names = fieldnames(ROI_trial_b_trial);
ROI_subList = ROI_trial_b_trial.subject_id;
ROI_names(strcmp(ROI_names,'subject_id')) = [];
if length(ROI_names) > 1
    error(['There should be only 1 ROI selected, not ',num2str(length(ROI_names))])
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
potential_input_prm = {'NV','uncertainty','E_level','money_level',...
    'deltaMoney_level'};
which_input = listdlg('PromptString','please select input parameter',...
    'ListString',potential_input_prm);
input_prm_nm = potential_input_prm{which_input};

%% select which model to use (if relevant)
if ismember(input_prm_nm,{'NV','uncertainty'})
    needModeling = true;
    [mdlType, mdlN] = behavioral_model_selection;
else
    needModeling = false;
end

%% extract the parameters of interest
nRuns = 4;
nTrialsPerRun = 54;
nTrials = nTrialsPerRun*nRuns;
% input/mediator/output
[input_prm, ROI_mediator, choice_hE] = deal(NaN(nTrials, NS));

% bin variables
[input_f_input_low_ROI, choice_f_input_low_ROI,...
    input_f_input_high_ROI, choice_f_input_high_ROI] = deal(NaN(nBins, NS));

%% loop through subjects
for iS = 1:NS
sub_nm = subject_id{iS};
    subBehaviorFolder = [studyBehaviorFolder, 'CID',sub_nm, filesep, 'behavior',filesep];

    % extract confidence based on the model
    if needModeling == true
        switch mdlType
            case 'simple'
            [~, dataInferred] = logitfit_choices(computerRoot, study_nm, sub_nm,...
                0, 'levels', 6, 6);
        case 'bayesian'
            error('please update for bayesian model');
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
                    uncertainty_tmp = - dataInferred.confidenceFitted.(['mdl_',mdlN]).(run_nm_bis); % revert sign to transform confidence into uncertainty
                otherwise
                    error('not ready yet');
            end
        end
        
        %% choice
        choice_LR_tmp = choiceAndPerf_tmp.choice;
        choice_highE_tmp = NaN(1,length(choice_LR_tmp));
        choice_highE_tmp(choice_LR_tmp == -defaultSide_tmp) = 1;
        choice_highE_tmp(choice_LR_tmp == defaultSide_tmp) = 0;
        
        %% extract input behavioral variable
        switch input_prm_nm
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
        choice_hE(runTrials_idx, iS) = choice_highE_tmp;
    end % run loop
    

    %% extract the bins
    
end % subject loop