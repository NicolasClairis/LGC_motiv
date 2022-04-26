function[totalGain] = compute_gains_clinicalTrial()
%[totalGain] = compute_gains_clinicalTrial()
% function to compute the total final gains after the three visits.
% 
% OUTPUTS
% totalGain: amount of money the participant should be paid in the end

%% enter subject ID
% subject = 'CID0XX';
subjectId = [];
while isempty(subjectId)
    info = inputdlg({'Subject CID (XXX)'});
    subjectId = ['CID',info{[1]}];
end

%% go to subject path
cd ..
rootPath = [pwd, filesep];
subPath = [rootPath,'LGC_Motiv_results',filesep,subjectId,filesep,'behavior',filesep];

%% load data of the three visits
totalGain = 0;
% compute money gained during the task
nVisits = 3;
for iVisit = 1:nVisits
    nRuns = 4;
    for iRun = 1:nRuns
        filenm = ls([subPath,subjectId,'_session',num2str(iRun),'_*_task_behavioral_visit',num2str(iVisit),'_tmp.mat']);
        sessionGain = getfield(getfield(load([subPath,filenm],'summary'),'summary'),'totalGain');
        totalGain = totalGain + sessionGain(end);
    end % run loop
end % visit loop

% add money acquired during IP
filenm = ls([subPath,'delta_IP_',subjectId,'.mat']);
totalGain = totalGain + getfield(getfield(load([subPath,filenm]),'IP_variables'),'totalGain');


% add the other resources
questionnairesMoney = 10; % like 1h exp
moneyPerVisit = 20*nVisits;
supplementationMoney = 150*(nVisits - 1);
MVCmoney = 48*nVisits;
totalGain = totalGain + questionnairesMoney + supplementationMoney + moneyPerVisit + MVCmoney;

end % function