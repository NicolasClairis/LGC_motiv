function[choices] = count_choices_perElvl(subBehaviorFolder, sub_nm, run_nm, task_fullName)
% [choices] = count_choices_perElvl(subBehaviorFolder, sub_nm, run_nm, task_fullName)
% count_choices_perElvl will count the proportion of high vs low effortful
% choices for each effort level (E1/E2/E3).
%
% INPUTS
% subBehaviorFolder: subject behavioral folder
%
% sub_nm: subject name
%
% run_nm: string with run name
%
% task_fullName: string indicating task full name 'mental'/'physical'
%
% OUTPUTS
% choices: structure with information about proportion of high vs low
% effort choices for each level of effort difficulty

%% general parameters
nTrials = 54;

%% count proportion of choices
[choices.E1.lowCh, choices.E2.lowCh, choices.E3.lowCh,...
    choices.E1.highCh, choices.E2.highCh, choices.E3.highCh,...
    choices.E1.highCh_percentage, choices.E2.highCh_percentage, choices.E3.highCh_percentage] = deal(NaN(1,nTrials));
% extract information about effort level
hE_level = extract_hE_level(subBehaviorFolder, sub_nm, run_nm, task_fullName);
E1_trials = hE_level == 1;
E2_trials = hE_level == 2;
E3_trials = hE_level == 3;
% extract choices
choice_highE = extract_choice_hE(subBehaviorFolder, sub_nm, run_nm, task_fullName);
% extract choices for each effort level
% effort 1
choices.E1.lowCh = sum(choice_highE(E1_trials) == 0);
choices.E1.highCh = sum(choice_highE(E1_trials) == 1);
choices.E1.highCh_percentage = choices.E1.highCh./choices.E1.lowCh;
% effort 2
choices.E2.lowCh = sum(choice_highE(E2_trials) == 0);
choices.E2.highCh = sum(choice_highE(E2_trials) == 1);
choices.E2.highCh_percentage = choices.E2.highCh./choices.E2.lowCh;
% effort 3
choices.E3.lowCh = sum(choice_highE(E3_trials) == 0);
choices.E3.highCh = sum(choice_highE(E3_trials) == 1);
choices.E3.highCh_percentage = choices.E3.highCh./choices.E3.lowCh;

%% register if each condition has more than 1 trial
if (choices.E1.lowCh > 1) && (choices.E1.highCh > 1) &&...
        (choices.E2.lowCh > 1) && (choices.E2.highCh > 1) &&...
        (choices.E3.lowCh > 1) && (choices.E3.highCh > 1)
    choices.allTrialsFullFilled = 1;
else
    choices.allTrialsFullFilled = 0;
end % both low and high effort choices present for each effort level (E1/E2/E3)
end % function