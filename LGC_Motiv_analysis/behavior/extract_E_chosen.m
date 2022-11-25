function[E_chosen] = extract_E_chosen(subBehaviorFolder, sub_nm, run_nm, task_fullName)
% [E_chosen] = extract_E_chosen(subBehaviorFolder, sub_nm, run_nm, task_fullName)
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
% E_chosen: 1*nTrials vector with information about the effort level 
% chosen for the current study, subject and run.

%% load data
behaviorStruct = load([subBehaviorFolder,...
    'CID',sub_nm,'_session',run_nm,'_',task_fullName,...
    '_task.mat']);
choiceOptions = behaviorStruct.choice_opt;
switch task_fullName
    case 'mental'
        choiceAndPerf = behaviorStruct.mentalE_perf;
    case 'physical'
        choiceAndPerf = behaviorStruct.physicalPerf;
end
%% choice was left or right?
choice_LR = choiceAndPerf.choice;
% remove confidence info from choice:
choice_LR(choice_LR == 2) = 1;
choice_LR(choice_LR == -2) = -1;
%% extract high effort level
E_chosen = (choiceOptions.E.left).*(choice_LR == -1) +...
            (choiceOptions.E.right).*(choice_LR == 1);
        
end % function