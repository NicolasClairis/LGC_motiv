function [trainingChoiceOptions, nTrainingTrials] = training_options(taskTrainingCond, n_R_levels, n_E_levels, R_money)
% [trainingChoiceOptions, nTrainingTrials] = training_options(taskTrainingCond, n_R_levels, n_E_levels)
% design of the reward, effort and punishment options for the learning
% phase
%
% INPUTS
% taskTrainingCond
% 'R': reward
% 'P': punishment
% 'RP': both Reward and Punishment options
%
% n_R_levels: number of reward levels
%
% n_E_levels: number of effort levels
%
% R_money: structure with monetary amount corresponding to each reward and
% punishment level
%
% OUTPUTS
% trainingChoiceOptions: structure with reward and effort level for each
% training trial + reward or punishment trial
%
% nTrainingTrials: number of training trials

%% define options and reward or punishment trials
if n_R_levels == 3 && n_E_levels == 3
    switch taskTrainingCond
        case 'R'
            trainingChoiceOptions.R.left    = [1 3 2];
            trainingChoiceOptions.R.right   = [3 1 1];
            trainingChoiceOptions.E.left    = [1 1 3];
            trainingChoiceOptions.E.right   = [3 2 2];
        case 'P'
            % would be good to map reward and punishment so that there is
            % no discrepancy in the average reward and effort levels
            % between the two training sessions
            trainingChoiceOptions.R.left    = [1 2 2];
            trainingChoiceOptions.R.right   = [3 1 3];
            trainingChoiceOptions.E.left    = [3 1 1];
            trainingChoiceOptions.E.right   = [1 2 3];
        case 'RP'
            trainingChoiceOptions.R.left    = [1 3 2 1 2 2];
            trainingChoiceOptions.R.right   = [3 1 1 3 1 3];
            trainingChoiceOptions.E.left    = [1 1 3 3 1 1];
            trainingChoiceOptions.E.right   = [3 2 2 1 2 3];
    end
else
    error(['Please determine a training sequence for when there are ',num2str(n_R_levels),' reward levels ',...
        ' and ', num2str(n_E_levels),' effort levels']);
end

%% number of trials
nTrainingTrials = size(trainingChoiceOptions.E.left, 2);

%% reward or punishment
switch taskTrainingCond
    case 'R'
        R_or_P = repmat({'R'},1,nTrainingTrials);
    case 'P'
        R_or_P = repmat({'P'},1,nTrainingTrials);
    case 'RP'
        R_or_P = [repmat({'R'},1,nTrainingTrials/2), repmat({'P'},1,nTrainingTrials/2)];
        % use the same sequence for all participants
        if n_R_levels == 3 && n_E_levels == 3 && nTrainingTrials == 12
            rand_RP = [4 9 2 6 5 11 10 12 1 7 8 3];
        elseif n_R_levels == 3 && n_E_levels == 3 && nTrainingTrials == 6
            rand_RP = [2 3 4 5 1 6];
        else
            error(['Please define a pre-determined training sequence for reward and punishment for when ',num2str(n_R_levels),' reward levels ',...
        ' and ', num2str(n_E_levels),' effort levels and ',num2str(nTrainingTrials),' trials']);
        end
        R_or_P = R_or_P(rand_RP);
        
        % re-order options based on RP order randomisation
        trainingChoiceOptions.R.left = trainingChoiceOptions.R.left(rand_RP);
        trainingChoiceOptions.E.left = trainingChoiceOptions.E.left(rand_RP);
        trainingChoiceOptions.R.right = trainingChoiceOptions.R.right(rand_RP);
        trainingChoiceOptions.E.right = trainingChoiceOptions.E.right(rand_RP);
end
% store reward/punishment condition for each trial
trainingChoiceOptions.R_or_P = R_or_P;

%% convert data in monetary amounts
if exist('R_money','var') && ~isempty(R_money) % you also need this function for the preparation of the timings, no need to compute this in this case
    [trainingChoiceOptions.monetary_amount.left] = reward_level_to_moneyAmount_converter(trainingChoiceOptions.R.left, R_money, R_or_P);
    [trainingChoiceOptions.monetary_amount.right] = reward_level_to_moneyAmount_converter(trainingChoiceOptions.R.right, R_money, R_or_P);
end

end % function