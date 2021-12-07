function[betas, choices] = logitfit_choices(computerRoot, study_nm, sub_nm, figDisp)
% 
% logitfit_choices will perform a logistic regression on the choices
% performed by the participants
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
% OUTPUTS
% betas: structure with betas
%
% choices: structure with bins for choices

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

%% extract runs
[runsStruct] = runs_definition(study_nm, sub_nm, 'behavior');

%% define R/P/E values
P_levels = [-3, -2, -1];
R_levels = [1, 2, 3];
E_levels = [1, 2, 3];

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
    
    for iRun = 1:runs.nb_runs.(task_id)
        
        %% load data
        behaviorStruct = load([subBehaviorFolder,...
            'CID',sub_nm,'_session',num2str(iRun),'_',task_fullName,...
            '_task.mat']);
        choiceOptions = behaviorStruct.choice_opt;
        switch task_id
            case 'Ep'
                choice_LR = behaviorStruct.physicalPerf.choice;
            case 'Em'
                choice_LR = behaviorStruct.mentalE_perf.choice;
        end
        % remove confidence information
        choice_LR(choice_LR == -2) = -1;
        choice_LR(choice_LR == 2) = 1;
        
        % load actual values for the specific subject
        P_vals = [R_money.P_1, R_money.P_2, R_money.P_3];
        R_vals = [R_money.R_1, R_money.R_2, R_money.R_3];
        
        % deltas
        deltaR_Ep = ;
        deltaEp = ;
        deltaP_Ep = ;
        
        deltaR_Em = ;
        deltaEm = ;
        deltaP_Em = ;
        
        choice_nonDef = choice_LR == ;
        
    end % run loop
    
    %% perform the fit
    
    
    %% figures
    if figDisp == 1
        %% display choice = f(net value)
        fig;
        
    end
    
end % physical/mental

end % function