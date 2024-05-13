function[choice_highE] = extract_choice_hE_bis(subBehaviorFolder,...
    sub_nm, run_nm, task_fullName)
% [choice_highE] = extract_choice_hE_bis(subBehaviorFolder,...
%   sub_nm, run_nm, task_fullName)
% extract_choice_hE_bis is like extract_choice_hE, but instead of coding
% choice as a binary variable (0=low E choice; 1=high E choice), it will
% encode it as 4 levels depending on confidence:
% (0)       = low E + high Confidence
% (0.25)    = low E + low Confidence
% (0.75)    = high E + low Confidence
% (1)       = high E + high Confidence
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

% load choice
[choice_highE_binary] = extract_choice_hE(subBehaviorFolder,...
    sub_nm, run_nm, task_fullName);
% load confidence
[conf_rtg] = extract_confidence_rating(subBehaviorFolder,...
    sub_nm, run_nm, task_fullName);

% recompute choice with 4 levels
choice_highE = choice_highE_binary;
high_E_low_Conf_trials = (choice_highE_binary == 1).*(conf_rtg == 0) == 1;
low_E_low_Conf_trials = (choice_highE_binary == 0).*(conf_rtg == 0) == 1;
choice_highE(high_E_low_Conf_trials) = 0.75;
choice_highE(low_E_low_Conf_trials) = 0.25;
end % function