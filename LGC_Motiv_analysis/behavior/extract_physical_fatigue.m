function[sumPrevAUC_N, sumPrevAUC] = extract_physical_fatigue(subBehaviorFolder, sub_nm, run_nm, task_fullName)
% [sumPrevAUC_N, sumPrevAUC] = extract_physical_fatigue(subBehaviorFolder, sub_nm, run_nm, task_fullName)
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
% sumPrevAUC_N: 1*nTrials vector with information about the cumulated effort
% performed until trial n, based on the force expressed in Newtons.
%
% sumPrevAUC: 1*nTrials vector with information about the cumulated effort
% performed until trial n, based on the force expressed in Voltage (before 
% conversion in Newtons).

%% sanity check
if ~strcmp(task_fullName,'physical')
    error(['task_fullName = ',task_fullName,' instead of ''physical'' as it should']);
end
%% load data
[~, AUC, ~, ~,...
    AUC_N] = extract_grip_force(subBehaviorFolder, sub_nm, run_nm);
AUC_perTrial = AUC.allTrials; % data in Voltage
AUC_N_perTrial = AUC_N.allTrials; % data in Newtons
nTrialsPerRun = 54;
[sumPrevAUC, sumPrevAUC_N] = deal(NaN(1,nTrialsPerRun));
% initialize at 0 (no force performed previously)
sumPrevAUC(1) = 0;
sumPrevAUC_N(1) = 0;
for iTrial = 2:nTrialsPerRun
    sumPrevAUC(iTrial) = sumPrevAUC(iTrial-1) + AUC_perTrial(iTrial - 1);
    sumPrevAUC_N(iTrial) = sumPrevAUC_N(iTrial-1) + AUC_N_perTrial(iTrial - 1);
end % trial loop
end % function