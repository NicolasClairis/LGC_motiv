function[perf] = extract_perf(subBehaviorFolder, sub_nm, run_nm, task_fullName)
% [perf] = extract_perf(subBehaviorFolder, sub_nm, run_nm)
%
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
% perf: performance in terms of percentage of the circle that was completed
% by the end of the trial
%
% N. Clairis - november 2022

%% load the data
behaviorStruct = load([subBehaviorFolder,...
    'CID',sub_nm,'_session',run_nm,'_',task_fullName,'_task_behavioral_tmp.mat']);
perf = behaviorStruct.summary.percentagePerf.*100;


end % function