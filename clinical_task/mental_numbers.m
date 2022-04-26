function[mental_nbers_per_trial] = mental_numbers(n_trials)
%[mental_nbers_per_trial] = mental_numbers(n_trials)
% mental_numbers will pre-define all the numbers appearing on screen for every participant.
% The script ensures that no sequence is repeated so that they can't expect
% a given number to appear after seeing another.
%
% INPUTS
% n_trials: number of total trials
%
% OUTPUTS
% mental_nbers_per_trial: matrix with all the numbers to display for each
% trial
%
%

%% define numbers to be used
n_to_use_taskSwitch = [1, 2, 3, 4, 6, 7, 8, 9]; % all numbers between 1
% and 9 but 5 to have equal amount of numbers <5 and >5
n_nbers_to_use = length(n_to_use_taskSwitch);

%% define matrix with all the possible numbers
n_sequences_per_trial = 50; % for security prepare a lot just in case the subject keeps committing errors
mental_nbers_per_trial = NaN(n_trials, n_nbers_to_use*n_sequences_per_trial);

%% extract all possible permutations
all_possible_perms = perms(1:n_nbers_to_use);
n_possiblePerms = size(all_possible_perms, 1);

% randomize the order of the permutations
if n_nbers_to_use == 8 % in this case: lots of permutation, you can use them
    rdm_order = randperm(n_possiblePerms);
end
all_possible_perms_rdm_order = all_possible_perms(rdm_order,:);


%% prepare each trial
jIndex = 0;
for iTrial = 1:n_trials
    for iSeq_per_trial = 1:n_sequences_per_trial
        questions_index = (1:n_nbers_to_use) + n_nbers_to_use*(iSeq_per_trial - 1);
        
        % remove possibility that the same number appears twice
        adjacentNbOk = 0;
        while adjacentNbOk == 0
            jIndex = jIndex + 1;
            if iSeq_per_trial == 1 ||...
                    (mental_nbers_per_trial(iTrial, n_nbers_to_use*(iSeq_per_trial-1)) ~= n_to_use_taskSwitch( all_possible_perms_rdm_order(jIndex,1)))
                % last element of previous sequence and first element of next sequence have to be different
                mental_nbers_per_trial(iTrial, questions_index) = n_to_use_taskSwitch( all_possible_perms_rdm_order(jIndex,:));
                adjacentNbOk = 1;
            end
        end
    end
end % trial loop

end % function