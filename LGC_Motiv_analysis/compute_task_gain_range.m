function[minDefaultNoEffortLoss, minNonDefaultNoEffortLoss, minDefaultGains, minGains, maxGains] = compute_task_gain_range()
% [minDefaultGains, minNonDefaultNoEffortLoss, minGains, maxGains] = compute_task_gains()
%
% OUTPUTS
% minDefaultNoEffortLoss: sum of money lost at the end if the participant never
% performed the effort and always picked the default option
%
% minNonDefaultNoEffortLoss: sum of money lost at the end if the
% participant never performed the effort and always picked the non-default
% option with the highest possible IP (=highest amount of money to lose)
%
% minDefaultGains: sum of money obtained at the end if participant always
% picked the default option
%
% minGains: minimal amount of money obtained if participant always picked
% the non-default option with the lowest possible IP (0.01 chf)
%
% maxGains: maximal amount of money obtained if participant always picked
% the non-default option with the highest possible IP (0.50 chf)

%% define general parameters
IPdata.baselineR = 0.5;
IPdata.baselineP = 0.5;
n_R_levels = 4;
punishment_yn = 'yes';
taskType = 'mental';
nRuns = 4;
nTrialsPerRun = 54;
% load design matrix
designMatrix = getfield(load('DaBestDesignMat.mat'),'bestMatrix');
R_or_P_trial_logic = strcmp(designMatrix.R_or_P,'R');
R_level = (designMatrix.R.left).*(designMatrix.default_LR == 1) + (designMatrix.R.right).*(designMatrix.default_LR == -1);

%% min range (if always taking the default option)
minDefaultGains = (IPdata.baselineR*(nTrialsPerRun/2) - IPdata.baselineP*(nTrialsPerRun/2))*nRuns;
% same but if no effort is performed at all
minDefaultNoEffortLoss = (-2*IPdata.baselineP*(nTrialsPerRun/2))*nRuns;

%% min range with lowest possible IP (if always taking the non-default)
IPdata.mentalDeltaIP = 0.01;
[R_money_min] = R_amounts_IP(n_R_levels, punishment_yn, IPdata, taskType);
minMoneyAmountPerTrialPerRun = NaN(1,nTrialsPerRun);
for iTrial = 1:nTrialsPerRun
    if R_or_P_trial_logic(iTrial) == true
        minMoneyAmountPerTrialPerRun(iTrial) = R_money_min.(['R_',num2str(R_level(iTrial))]);
    elseif R_or_P_trial_logic(iTrial) == false
        minMoneyAmountPerTrialPerRun(iTrial) = -R_money_min.(['P_',num2str(R_level(iTrial))]);
    end
end % trial loop
minGains = sum(minMoneyAmountPerTrialPerRun)*nRuns;

%% max range with highest possible IP (if always taking the non-default)
IPdata.mentalDeltaIP = 0.5;
[R_money_max] = R_amounts_IP(n_R_levels, punishment_yn, IPdata, taskType);
maxMoneyAmountPerTrialPerRun = NaN(1,nTrialsPerRun);
for iTrial = 1:nTrialsPerRun
    if R_or_P_trial_logic(iTrial) == true
        maxMoneyAmountPerTrialPerRun(iTrial) = R_money_max.(['R_',num2str(R_level(iTrial))]);
    elseif R_or_P_trial_logic(iTrial) == false
        maxMoneyAmountPerTrialPerRun(iTrial) = -R_money_max.(['P_',num2str(R_level(iTrial))]);
    end
end % trial loop
maxGains = sum(maxMoneyAmountPerTrialPerRun)*nRuns;

%% minimal possible amount of money to lose: highest possible IP (0.50 chf) but never performing the effort and always picking the default option (ie maximizing losses)
maxLossMoneyAmountPerTrialPerRun = NaN(1,nTrialsPerRun);
for iTrial = 1:nTrialsPerRun
    if R_or_P_trial_logic(iTrial) == true
        maxLossMoneyAmountPerTrialPerRun(iTrial) = 0;
    elseif R_or_P_trial_logic(iTrial) == false
        maxLossMoneyAmountPerTrialPerRun(iTrial) = -2*R_money_max.P_0;
    end
end % trial loop
minNonDefaultNoEffortLoss = sum(maxLossMoneyAmountPerTrialPerRun)*nRuns;

end % function