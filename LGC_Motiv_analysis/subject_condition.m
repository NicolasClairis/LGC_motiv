function [condition] = subject_condition()
%[condition] = subject_condition()
% subject_condition will display a window allowing you to select which
% participants and runs to include in your analysis based on the different
% possibilities for fMRI movement or not and task saturation or not (for
% the choices).
%
% OUTPUT
% condition: string with name of the condition

listPossibleConditions = {'fullList',...
    'fMRI_noSatTaskSub_noSatRun',...
    'fMRI_noSatTaskSub_noSatRun_choiceSplit_Elvl',...
    'fMRI_noSatTaskSub_noMoveSub_noSatRun',...
    'fMRI_noSatTaskSub_noMoveSub_noSatRun_noMoveRun',...
    'fMRI_noSatTaskSub',...
    'behavior',...
    'behavior_noSatTaskSub','behavior_noSatRunSub',...
    'behavior_noSatTask','behavior_noSatRun',...
    'behavior_noSatRun_bayesianMdl',...
    'behavior_noSatTask_bayesianMdl',...
    'fMRI',...
    'fMRI_noMoveSub','fMRI_noMoveSub_bis','fMRI_noMoveSub_ter',...
    'fMRI_noMove_bis','fMRI_noMove_ter',...
    'fMRI_noSatRunSub',...
    'fMRI_noSatTaskSub_noMove_bis_Sub',...
    'fMRI_noSatTask_noMove_bis',...
    'fMRI_noSatRun','fMRI_noSatTask',...
    'fMRI_noSatRun_choiceSplit_Elvl',...
    'fMRI_noSatRun_bayesianMdl',...
    'fMRI_noSatTask_bayesianMdl',...
    'fMRI_noSatTask_noMove_bis_bayesianMdl'};

cond_idx = listdlg('ListString',listPossibleConditions);
condition = listPossibleConditions{cond_idx};
end % function