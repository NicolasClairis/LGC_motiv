function [condition] = subject_condition()
%[condition] = subject_condition()
% subject_condition will display a window allowing you to select which
% participants and runs to include in your analysis based on the different
% possibilities for fMRI movement or not and task saturation or not (for
% the choices).
%
% OUTPUT
% condition: string with name of the condition

listPossibleConditions = {'behavior',...
    'behavior_noSatTask','behavior_noSatRun',...
    'fMRI','fMRI_no_move','fMRI_no_move_bis',...
    'fMRI_noSatRun','fMRI_noSatTask'};

cond_idx = listdlg('ListString',listPossibleConditions);
condition = listPossibleConditions{cond_idx};
end % function