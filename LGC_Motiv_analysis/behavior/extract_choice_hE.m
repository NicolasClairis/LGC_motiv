function[choice_highE] = extract_choice_hE(subBehaviorFolder,...
    sub_nm, run_nm, task_fullName)
% [choice_highE] = extract_choice_hE(subBehaviorFolder,...
%   sub_nm, run_nm, task_fullName)
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
% choice_highE: 1*nTrials vector with information about when subject sub_nm
% selected the high effort option (=1), the low effort option (=0) or no
% option at all (=NaN) for the given task (task_fullName) and run (run_nm).

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
%% default side
defaultSide = choiceOptions.default_LR;
%% choice
choice_LR = choiceAndPerf.choice;
% remove confidence info from choice:
choice_LR(choice_LR == 2) = 1;
choice_LR(choice_LR == -2) = -1;
% extract high effort choice
choice_highE = NaN(1,length(choice_LR));
choice_highE(choice_LR == -defaultSide) = 1;
choice_highE(choice_LR == defaultSide) = 0;
end % function