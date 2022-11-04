function[hE_level] = extract_hE_level(subBehaviorFolder, sub_nm, run_nm, task_fullName)
% [hE_level] = extract_choice_hE(subBehaviorFolder, sub_nm, run_nm, task_fullName)
%
% INPUTS
% subBehaviorFolder: folder where data is stored
%
% sub_nm: string with subject CID name
%
% run_nm: string with run name
%
% task_fullName: task full name 'mental'/'physical'
%
% OUTPUTS
% hE_level: 1*nTrials vector with information about the effort level 
% associated to the high effort option for the current study, subject and
% run.

%% load data
behaviorStruct = load([subBehaviorFolder,...
    'CID',sub_nm,'_session',run_nm,'_',task_fullName,...
    '_task.mat']);
choiceOptions = behaviorStruct.choice_opt;

%% default side
defaultSide = choiceOptions.default_LR;
%% extract high effort level
hE_level = (choiceOptions.E.left).*(defaultSide == 1) +...
            (choiceOptions.E.right).*(defaultSide == -1);
        
end % function