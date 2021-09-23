% check monetary amounts that you can win by performing the task

bestMatrix = getfield(load('DaBestDesignMat.mat'),'bestMatrix');

R_or_P = strcmp(bestMatrix.R_or_P,'R');
R_trials = R_or_P == 1;
P_trials = R_or_P == 0;
R_left = bestMatrix.monetary_amount.left(R_trials);
R_right = bestMatrix.monetary_amount.right(R_trials);
P_left = bestMatrix.monetary_amount.left(P_trials);
P_right = bestMatrix.monetary_amount.right(P_trials);

R_max = max([R_left; R_right]);
R_min = min([R_left; R_right]);
P_min = min([P_left; P_right]);
P_max = max([P_left; P_right]);

maxMoneyPerSession = sum(R_max) - sum(P_min);
minMoneyPerSession = sum(R_min) - sum(P_max);
nSessions = 4;

disp(['min money per session = ',num2str(minMoneyPerSession),' up to max money per session = ',num2str(maxMoneyPerSession)]);
disp(['min money = ',num2str(minMoneyPerSession*nSessions),' up to max money = ',num2str(maxMoneyPerSession*nSessions)]);