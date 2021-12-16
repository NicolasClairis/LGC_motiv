function[betas, RTstruct] = fit_RT(computerRoot, study_nm, sub_nm,...
    figDisp, dispMoneyOrLevels, n_NV_bins, n_trialN_bins)
% [betas, RTstruct] = fit_RT(computerRoot, study_nm, sub_nm,...
%       figDisp, dispMoneyOrLevels, n_NV_bins, n_trialN_bins)
% logitfit_RT will perform a linear regression on the reaction times (RT)
% of the participants
%
% INPUTS
% computerRoot: pathway where data is
% 
% study_nm: study name
%
% sub_nm: subject number id 'XXX'
%
% figDisp: display individual figure (1) or not (0)
%
% dispMoneyOrLevels: display actual money ('money') or reward levels
% ('levels')
%
% n_NV_bins: number of bins for net value
%
% n_trialN_bins: number of bins for fatigue
%
% OUTPUTS
% betas: structure with betas
%
% RTstruct: structure with bins for RT

%% if root not defined => ask for it
if ~exist('computerRoot','var') || isempty(computerRoot)
    computerRoot = LGCM_root_paths;
end

%% working directories
subBehaviorFolder = [computerRoot, filesep, study_nm, filesep,...
    'CID',sub_nm, filesep, 'behavior', filesep];

%% by default, display individual figure
if ~exist('figDisp','var') || isempty(figDisp)
    figDisp = 1;
    disp(['figDisp was not defined in the inputs so that by default ',...
        'figures are displayed for each individual.']);
end

%% by default, display monetary levels instead of actual monetary amounts
if ~exist('dispMoneyOrLevels','var') || isempty(dispMoneyOrLevels)
    dispMoneyOrLevels = 'levels';
end

%% extract runs
[runsStruct] = runs_definition(study_nm, sub_nm, 'behavior');
nRuns = length(runsStruct.tasks);
runs_Ep = strcmp(runsStruct.tasks,'Ep');
runs_Em = strcmp(runsStruct.tasks,'Em');

%% define R/P/E values
money_levels = [-3, -2, -1, 1, 2, 3];
nMoneyLevels = length(money_levels);
[actualMoney_values.Ep,...
    actualMoney_values.Em]   = deal(NaN(1,nMoneyLevels));
E_levels = [1, 2, 3];
nELevels = length(E_levels);

if ~exist('n_NV_bins','var') || isempty(n_NV_bins)
    n_NV_bins = 6;
end

if ~exist('n_trialN_bins','var') || isempty(n_trialN_bins)
    n_trialN_bins = 6;
end

%% load choices fit
[betas, choices] = logitfit_choices(computerRoot, study_nm, sub_nm,...
    0, dispMoneyOrLevels, n_NV_bins, n_trialN_bins);

%% initialize variables of interest
nTrialsPerRun = 54;
nTrialsPerRPConditionPerRun = nTrialsPerRun/2;
nRunsPerTask = 2;
nTrials = nTrialsPerRun*nRunsPerTask;
nTrialsPerRPCond = nTrialsPerRPConditionPerRun*nRunsPerTask;

%% loop through physical and mental
for iPM = 1:2
    switch iPM
        case 1
            task_id = 'Ep';
            task_fullName = 'physical';
        case 2
            task_id = 'Em';
            task_fullName = 'mental';
    end
    
    [RT_perTrial.(task_id),...
        NV.(task_id).mdl1,...
        NV.(task_id).mdl2,...
        NV.(task_id).mdl3,...
        NV.(task_id).mdl4,...
        trialN.(task_id),...
        R_nonDef_valuePerTrial.(task_id),...
        R_default_valuePerTrial.(task_id),...
        P_default_valuePerTrial.(task_id),...
        P_nonDef_valuePerTrial.(task_id),...
        money_default_valuePerTrial.(task_id),...
        money_nonDef_valuePerTrial.(task_id),...
        money_nonDef_levelPerTrial.(task_id),...
        E_nonDef_levelPerTrial.(task_id),...
        E_default_levelPerTrial.(task_id),...
        confidence_perTrial.(task_id)] = deal(NaN(nTrials,1));
    
    
    jRun = 0;
    for iRun = 1:nRuns
        runToInclude = 0;
        switch task_id
            case 'Ep'
                if runs_Ep(iRun) == 1
                    jRun = jRun + 1;
                    runToInclude = 1;
                end
            case 'Em'
                if runs_Em(iRun) == 1
                    jRun = jRun + 1;
                    runToInclude = 1;
                end
        end
        
        if runToInclude == 1
            runTrials_idx = (1:nTrialsPerRun) + nTrialsPerRun*(jRun-1);
            runTrialsRP_idx = (1:nTrialsPerRPCond) + nTrialsPerRPCond*(jRun-1);
            
            %% load data
            behaviorStruct_tmp = load([subBehaviorFolder,...
                'CID',sub_nm,'_session',num2str(iRun),'_',task_fullName,...
                '_task.mat']);
            choiceOptions_tmp = behaviorStruct_tmp.choice_opt;
            onsets_tmp = behaviorStruct_tmp.onsets;
            switch task_id
                case 'Ep'
                    onsets_tmp = behaviorStruct_tmp.physicalPerf.onsets;
                case 'Em'
                    onsets_tmp = behaviorStruct_tmp.mentalE_perf.onsets;
            end
            RT_tmp = onsets_tmp.choice - onsets_tmp.dispChoiceOptions;
            
        end % run filter
        
    end % run loop
    
    %% figures
    if figDisp == 1
        pSize = 30;
        lWidth = 3;
        lWidth_borders = 1;
        
    end % figure display
        
end % physical/mental loop

%% extract output
