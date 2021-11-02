function[totalGain] = compute_gains()
%[totalGain] = compute_gains()
% function to compute the total final gains
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

%% load data
totalGain = 0;
% compute money gained during the task
nRuns = 4;
for iRun = 1:nRuns
    filenm = ls([subPath,subjectId,'_session',num2str(iRun),'_*_task_behavioral_tmp.mat']);
    sessionGain = getfield(getfield(load([subPath,filenm],'summary'),'summary'),'totalGain');
    totalGain = totalGain + sessionGain(end);
end
% add money acquired during IP
filenm = ls([subPath,'delta_IP_',subjectId,'.mat']);
totalGain = totalGain + getfield(getfield(load([subPath,filenm]),'IP_variables'),'totalGain');


% add the other resources
MRStunes = 50;
IRMftunes = 40;
MVCtune = 48;
timeTune = 40;
questionnairesTune = 40; % like 1h exp
% pilots
% MRStunes = 0;
% IRMftunes = 0;
% MVCtune = 12;
% timeTune = 0;
% questionnairesTune = 0; % like 1h exp
totalGain = totalGain + MRStunes + IRMftunes + MVCtune + timeTune + questionnairesTune;

end % function