function[RT] = extract_RT(subBehaviorFolder, sub_nm, run_nm, task_fullName)
% [RT] = extract_RT(subBehaviorFolder, sub_nm, run_nm, task_fullName)
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
% RT: 1*nTrials vector with information about the reaction time (RT) for
% the choices of the current study, subject and run.

%% load data
behaviorStruct = load([subBehaviorFolder,...
    'CID',sub_nm,'_session',run_nm,'_',task_fullName,...
    '_task.mat']);
switch task_fullName
    case 'physical'
        onsets_tmp = behaviorStruct.physicalPerf.onsets;
    case 'mental'
        onsets_tmp = behaviorStruct.mentalE_perf.onsets;
end
RT = onsets_tmp.choice - onsets_tmp.dispChoiceOptions;
        
end % function