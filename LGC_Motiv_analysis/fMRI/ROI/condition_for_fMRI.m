function[condition_for_fMRI_trial_per_trial_extraction] = condition_for_fMRI(condition)
% [condition_for_fMRI_trial_per_trial_extraction] = condition_for_fMRI(condition)
% condition_for_fMRI will create a new condition to extract the ROI data
% extracted trial per trial (because for simplicity we did not filter the
% bad runs there)
%
% INPUTS
% condition: condition as input
%
% OUTPUTS
% condition_for_fMRI_trial_per_trial_extraction: either equal to condition
% or adapted to correspond to the condition without saturation.

switch condition
    case 'fMRI_noSatTaskSub_noSatRun'
        condition_for_fMRI_trial_per_trial_extraction = 'fMRI_noSatTaskSub';
    case {'fMRI_noSatRun',...
            'fMRI_noSatRun_choiceSplit_Elvl',...
            'fMRI_noSatRun_choiceSplit_Elvl_bis',...
            'fMRI_noSatRun_bayesianMdl'}
        condition_for_fMRI_trial_per_trial_extraction = 'fMRI';
    otherwise
        condition_for_fMRI_trial_per_trial_extraction = condition;
end

end % function