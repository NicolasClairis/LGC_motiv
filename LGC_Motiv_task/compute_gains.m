function[totalGain] = compute_gains()
%[totalGain] = compute_gains()
% function to compute the total final gains
% 
% OUTPUTS
% totalGain: amount of money the participant should be paid in the end

%% enter subject ID
% subject = 'CID0XX';
subjectid = 'CID0XX';

%% load data
totalGain = 0;
nRuns = 4;
for iRun = 1:nRuns
    filenm = ls([subjectid,'_session',num2str(iRun),'_*_task_behavioral_tmp.mat']);
    sessionGain = getfield(getfield(load(filenm,'summary'),'summary'),'totalGain');
    totalGain = totalGain + sessionGain(end);
end

% add the other resources
MRStune = 50;
MVCtune = 48;
timeTune = 40;
questionnairesTune = 10; % like 1h exp
totalGain = totalGain + MRStune + MVCtune + timeTune + questionnairesTune;

end % function