function[sumPrevAUC] = extract_physical_fatigue(subBehaviorFolder, sub_nm, run_nm, task_fullName)
% [sumPrevAUC] = extract_physical_fatigue(subBehaviorFolder, sub_nm, run_nm, task_fullName)
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
% sumPrevAUC: 1*nTrials vector with information about the cumulated effort
% performed until trial n.

%% sanity check
if ~strcmp(task_fullName,'physical')
    error(['task_fullName = ',task_fullName,' instead of ''physical'' as it should']);
end
%% load data
[~, ~, ~, ~,...
    AUC_N] = extract_grip_force(subBehaviorFolder, sub_nm, run_nm);
nTrialsPerRun = 54;
sumPrevAUC = NaN(1,nTrialsPerRun);
sumPrevAUC(1) = 0;
for iTrial = 2:nTrialsPerRun
    sumPrevAUC(iTrial) = sumPrevAUC(iTrial-1) + AUC_N(iTrial - 1);
end % trial loop
end % function