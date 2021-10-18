function[totalGain] = compute_gains()
%[totalGain] = compute_gains()
% function to compute the total final gains
% 
% OUTPUTS
% totalGain: amount of money the participant should be paid in the end

%% enter subject ID
% subject = 'CID0XX';
subjectid = [];
while isempty(subjectid)
    info = inputdlg({'Subject CID (XXX)'});
    subjectid = ['CID',info{[1]}];
end

%% go to subject path
cd ..
rootPath = [pwd, filesep];
subPath = [rootPath,'LGC_Motiv_results',filesep,subjectid,filesep,'behavior',filesep];

%% load data
totalGain = 0;
nRuns = 4;
for iRun = 1:nRuns
    filenm = ls([subPath,subjectid,'_session',num2str(iRun),'_*_task_behavioral_tmp.mat']);
    sessionGain = getfield(getfield(load([subPath,filenm],'summary'),'summary'),'totalGain');
    totalGain = totalGain + sessionGain(end);
end

% add the other resources
MRStunes = 50;
IRMftunes = 0;
MVCtune = 48;
timeTune = 40;
questionnairesTune = 10; % like 1h exp
totalGain = totalGain + MRStunes + IRMftunes + MVCtune + timeTune + questionnairesTune;

end % function