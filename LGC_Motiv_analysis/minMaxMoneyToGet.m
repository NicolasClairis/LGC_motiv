function[minTotalToObtain, maxTotalToObtain] = minMaxMoneyToGet(deltaIP)
%[minTotalToObtain, maxTotalToObtain] = minMaxMoneyToGet(deltaIP)
% minMaxMoneyToGet will compute the total money you can obtain by
% performing the LGC motivation task
%
% INPUTS
% deltaIP: define delta between indifference point and low effort option to
% get an idea about max/min performance depending on the IP
%
% OUTPUTS
% minTotalToObtain: minimum obtainable
%
% maxTotalToObtain: maximum obtainable

%% define number of trials
nTrialsPerSession = 54;
nSessions = 4;

%% check baseline difference in reward vs punishment amounts
n_R_levels = 4;
punishment_yn = 'yes';
IPdata.baselineR = 0.5;
IPdata.baselineP = 0.5;
IPdata.mentalDeltaIP = deltaIP;
[R_money] = R_amounts_IP(n_R_levels, punishment_yn,IPdata,'mental');

%% load design matrix
bestMatrix = getfield(load('DaBestDesignMat.mat'),'bestMatrix');
% extract reward vs punishment trials
R_or_P = bestMatrix.R_or_P;

% convert reward/punishment levels in money
[maxPerformer, minPerformer] = deal(NaN(1,nTrialsPerSession));
for iTrial = 1:nTrialsPerSession
    switch R_or_P{iTrial}
        case 'R'
            R_left = R_money.(['R_',num2str(bestMatrix.R.left(iTrial))]);
            R_right = R_money.(['R_',num2str(bestMatrix.R.right(iTrial))]);
            maxPerformer(iTrial) = max(R_left, R_right);
            minPerformer(iTrial) = min(R_left, R_right);
        case 'P'
            P_left = R_money.(['P_',num2str(bestMatrix.R.left(iTrial))]);
            P_right = R_money.(['P_',num2str(bestMatrix.R.right(iTrial))]);
            maxPerformer(iTrial) = -min(P_left, P_right);
            minPerformer(iTrial) = -max(P_left, P_right);
    end
end
sumMaxRP = sum(maxPerformer);
sumMinRP = sum(minPerformer);

%% compute total across the four sessions
maxTotalToObtain = sumMaxRP*nSessions;
minTotalToObtain = sumMinRP*nSessions;
end % function