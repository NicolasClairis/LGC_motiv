function [trainingChoiceOptions] = training_options(taskTrainingCond, n_R_levels, n_E_levels, R_money, nTrainingTrials)
% [trainingChoiceOptions] = training_options(taskTrainingCond, n_R_levels, n_E_levels, R_money, nTrainingTrials)
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
% nTrainingTrials: number of training trials
%
% OUTPUTS
% trainingChoiceOptions: structure with reward and effort level for each
% training trial + reward or punishment trial

%% define options and reward or punishment trials
if n_R_levels == 4 && n_E_levels == 4
    switch taskTrainingCond
        case 'R'
            R_or_P = repmat({'R'},1,nTrainingTrials);
        case 'P'
            R_or_P = repmat({'P'},1,nTrainingTrials);
        case 'RP'
            switch nTrainingTrials
                case 4
                    R.left    = [0 1 2 0];
                    R.right   = [1 0 0 3];
                    E.left    = [0 3 2 0];
                    E.right   = [1 0 0 2];
                    R_or_P = {'R','P','P','R'};
                    % side of default option
                    default_LR = [-1 1 1 -1];
                case 8
                    R.left    = [0 2 0 0 0 1 2 1];
                    R.right   = [1 0 1 3 3 0 0 0];
                    E.left    = [0 2 0 0 0 3 2 3];
                    E.right   = [1 0 1 2 2 0 0 0];
                    R_or_P = {'R','R','P','R','P','P','P','R'};
                    % side of default option
                    default_LR = [-1 1 -1 -1 -1 1 1 1];
                otherwise
                    error(['case where ',num2str(nTrainingTrials),' training trials not ready yet. Please define design matrix.'])
            end
        otherwise
            error(['Please determine a training sequence for when there are ',num2str(n_R_levels),' reward levels ',...
                ' and ', num2str(n_E_levels),' effort levels']);
    end
else
    switch taskTrainingCond
        case 'R'
            error(['Please determine a training sequence for when there are ',num2str(n_R_levels),' reward levels ',...
                ' and ', num2str(n_E_levels),' effort levels']);
        case 'P'
            error(['Please determine a training sequence for when there are ',num2str(n_R_levels),' reward levels ',...
                ' and ', num2str(n_E_levels),' effort levels']);
        case 'RP'
            error(['Please determine a training sequence for when there are ',num2str(n_R_levels),' reward levels ',...
                ' and ', num2str(n_E_levels),' effort levels']);
    end

end

%% convert data in monetary amounts
if exist('R_money','var') && ~isempty(R_money) % you also need this function for the preparation of the timings, no need to compute this in this case
    [trainingChoiceOptions.monetary_amount.left] = reward_level_to_moneyAmount_converter(R.left, R_money, R_or_P);
    [trainingChoiceOptions.monetary_amount.right] = reward_level_to_moneyAmount_converter(R.right, R_money, R_or_P);
end

%% store everything
trainingChoiceOptions.R = R;
trainingChoiceOptions.E = E;
trainingChoiceOptions.default_LR = default_LR;
trainingChoiceOptions.R_or_P = R_or_P;

end % function