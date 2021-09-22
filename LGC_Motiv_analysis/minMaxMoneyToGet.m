function[minTotalToObtain, maxTotalToObtain, baselineAmountBasedOnRPDifference] = minMaxMoneyToGet()
%[minTotalToObtain, maxTotalToObtain, baselineAmountBasedOnRPDifference] = minMaxMoneyToGet()
% minMaxMoneyToGet will compute the total money you can obtain by
% performing the LGC motivation task
%
% INPUTS
%
% OUTPUTS
%minTotalToObtain: minimum obtainable
%
% maxTotalToObtain: maximum obtainable
%
% baselineAmountBasedOnRPDifference: baseline difference between rewards
% and punishments which provides with baseline money

%% define number of trials
nTrialsPerSession = 48;
nTrialsPerConditionPerSession = nTrialsPerSession/2;
nSessions = 4;
nTrials = nTrialsPerSession*nSessions;

%% check baseline difference in reward vs punishment amounts
n_R_levels = 3;
punishment_yn = 'yes';
[R_money] = R_amounts(n_R_levels, punishment_yn);
RP_difference = R_money.R_1 - R_money.P_1;
baselineAmountBasedOnRPDifference = RP_difference*nTrials/2;

%% load design matrix
bestMatrix = getfield(load('DaBestDesignMat.mat'),'bestMatrix');
% extract reward vs punishment trials
R_or_P = strcmp(bestMatrix.R_or_P,'R');
Rtrials = R_or_P == 1;
Ptrials = R_or_P == 0;
% extract reward and punishment trials
RmoneyLeft = bestMatrix.R.left(Rtrials);
RmoneyRight = bestMatrix.R.right(Rtrials);
PmoneyLeft = bestMatrix.R.left(Ptrials);
PmoneyRight = bestMatrix.R.right(Ptrials);
% extract max/min reward and punishment level per trial
maxR = nanmax( [RmoneyLeft; RmoneyRight]);
minR = nanmin( [RmoneyLeft; RmoneyRight]);
minP = nanmin( [PmoneyLeft; PmoneyRight]);
maxP = nanmax( [PmoneyLeft; PmoneyRight]);
% convert reward/punishment levels in money
[maxRmoney, minRmoney,...
    maxPmoney, minPmoney] = deal(NaN(nTrialsPerConditionPerSession,1));
for iTrial = 1:nTrialsPerConditionPerSession
    maxRmoney(iTrial) = R_money.(['R_',num2str(maxR(iTrial))]);
    minRmoney(iTrial) = R_money.(['R_',num2str(minR(iTrial))]);
    minPmoney(iTrial) = R_money.(['P_',num2str(minP(iTrial))]);
    maxPmoney(iTrial) = R_money.(['P_',num2str(maxP(iTrial))]);
end
sumMaxRP = (sum(maxRmoney) - sum(minPmoney))*nSessions;
sumMinRP = (sum(minRmoney) - sum(maxPmoney))*nSessions;

%% compute total
maxTotalToObtain = sumMaxRP;
minTotalToObtain = sumMinRP;
end % function