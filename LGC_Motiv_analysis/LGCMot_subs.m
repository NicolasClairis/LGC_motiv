function [subList, NS] = LGCMot_subs(subjectsToConsider)
%[subList, NS] = LGCMot_subs(subjectsToConsider)
% LGCMot_subs will determine the list of subjects to consider.
%
% INPUTS
% subjectsToConsider:
% 'behavioral_pilots': list of behavioral pilots used to confirm the task
% 'study1': list of participants having performed the main task
%
% OUTPUTS
% subList: cell with the list of participants
%
% NS: number of subjects included in subList


switch subjectsToConsider
    case 'behavioral_pilots'
        subList = {'201','202','203','204','205','206','207','208','209'};
    case 'study1'
        subList = {'074','036','095','064'};
end
NS = length(subList);

end % function