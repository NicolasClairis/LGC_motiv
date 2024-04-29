function[prevEfficacy_with2first,...
    prevEfficacy_pureNback,...
    prevEfficacy_bis_with2first,...
    prevEfficacy_bis_pureNback] = extract_mental_previous_efficacy(subBehaviorFolder, sub_nm, run_nm, task_fullName)
% [prevEfficacy_with2first,...
%     prevEfficacy_pureNback,...
%     prevEfficacy_bis_with2first,...
%     prevEfficacy_bis_pureNback] = extract_mental_previous_efficacy(subBehaviorFolder, sub_nm, run_nm, task_fullName)
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
% prevEfficacy_with2first: 1*nTrials vector with information about the efficacy expressed as
% (n_correct-n_errors)/RT_total_with2firstUseless in the previous trial.
%
% prevEfficacy_pureNback: 1*nTrials vector with information about the efficacy expressed as
% (n_correct-n_errors)/RTtotal_pureNback in the previous trial.
%
% prevEfficacy_bis_with2first: 1*nTrials vector with information about the efficacy expressed as
% (n_correct)/RT_total_with2firstUseless in the previous trial.
%
% prevEfficacy_bis_pureNback: 1*nTrials vector with information about the efficacy expressed as
% (n_correct)/RTtotal_pureNback in the previous trial.

%% sanity check
if ~strcmp(task_fullName,'mental')
    error(['task_fullName = ',task_fullName,' instead of ''mental'' as it should']);
end
%% load data
[~, ~, ~, ~,...
    ~,...
    ~,...
    efficacy_with2first,...
    efficacy_pureNback,...
    efficacy_bis_with2first,...
    efficacy_bis_pureNback] = extract_mental_perf(subBehaviorFolder, sub_nm, run_nm);


nTrialsPerRun = 54;
[prevEfficacy_with2first,...
    prevEfficacy_pureNback,...
    prevEfficacy_bis_with2first,...
    prevEfficacy_bis_pureNback] = deal(NaN(1,nTrialsPerRun));
% initialize at zero for first trial
prevEfficacy_with2first(1) = 0;
prevEfficacy_pureNback(1) = 0;
prevEfficacy_bis_with2first(1) = 0;
prevEfficacy_bis_pureNback(1) = 0;
% fill next trials
for iTrial = 2:nTrialsPerRun
    prevEfficacy_with2first(iTrial) = efficacy_with2first.allTrials(iTrial - 1);
    prevEfficacy_pureNback(iTrial) = efficacy_pureNback.allTrials(iTrial - 1);
    prevEfficacy_bis_with2first(iTrial) = efficacy_bis_with2first.allTrials(iTrial - 1);
    prevEfficacy_bis_pureNback(iTrial) = efficacy_bis_pureNback.allTrials(iTrial - 1);
end % trial loop

end % function