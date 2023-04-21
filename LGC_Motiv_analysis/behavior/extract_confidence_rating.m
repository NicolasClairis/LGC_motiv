function[conf_rtg] = extract_confidence_rating(subBehaviorFolder, sub_nm, run_nm, task_fullName)
% [conf_rtg] = extract_confidence_rating(subBehaviorFolder, sub_nm, run_nm, task_fullName)
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
% conf_rtg: 1*nTrials vector with information about the confidence rating
% for the current study, subject and run. Confidence is encoded as a binary
% variable (+1 when confidence was high, and 0 when confidence was low, NaN whenever no rating was provided).

%% load data
behaviorStruct = load([subBehaviorFolder,...
    'CID',sub_nm,'_session',run_nm,'_',task_fullName,...
    '_task.mat']);
switch task_fullName
    case 'mental'
        choiceAndPerf = behaviorStruct.mentalE_perf;
    case 'physical'
        choiceAndPerf = behaviorStruct.physicalPerf;
end
%% choice was left (-1/-2) or right (+1/+2)? low (-1/+1) or high (-2/+2) conf?
choice_LR = choiceAndPerf.choice;
%% extract confidence rating
conf_rtg = NaN(1,length(choice_LR));
conf_rtg(abs(choice_LR) == 2) = 1;
conf_rtg(abs(choice_LR) == 1) = 0;
        
end % function