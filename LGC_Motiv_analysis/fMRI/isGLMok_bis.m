function[] = isGLMok_bis(GLMprm, condition)
% isGLMok_bis(GLMprm, condition)
% isGLMok_bis will check whether the selected GLM is compatible with the
% condition selected for the first level analysis.
%
% INPUTS
% GLMprm: GLM parameters as extracted with which_GLM.m
%
% condition: condition used to filter subjects and runs as extracted with
% subject_condition.m
%
% See also which_GLM.m, runs_definition.m and subject_condition.m

%% if trials are split, you should be careful which condition is entered
% check for Ep choice period
if ((GLMprm.choice.Ep.splitPerE == 3) &&...
        ((GLMprm.choice.Ep.RP.lEch.E_varOption ~= 0) || (GLMprm.choice.Ep.RP.hEch.E_varOption ~= 0))) &&...
        ~ismember(condition,{'fMRI_noSatRun_choiceSplit_Elvl', 'fMRI_noSatTaskSub_noSatRun_choiceSplit_Elvl'})
    error(['condition ',condition,' is not compatible with splitting trials depending on the effort chosen and including the high effort regressor in Ep task.']);
end

% check for Em choice period
if ((GLMprm.choice.Em.splitPerE == 3) &&...
        ((GLMprm.choice.Em.RP.lEch.E_varOption ~= 0) || (GLMprm.choice.Em.RP.hEch.E_varOption ~= 0))) &&...
        ~ismember(condition,{'fMRI_noSatRun_choiceSplit_Elvl', 'fMRI_noSatTaskSub_noSatRun_choiceSplit_Elvl'})
    error(['condition ',condition,' is not compatible with splitting trials depending on the effort chosen and including the high effort regressor in Em task.']);
end

% check for Ep chosen period
if ((GLMprm.chosen.Ep.splitPerE == 3) &&...
        ((GLMprm.chosen.Ep.RP.lEch.E_varOption ~= 0) || (GLMprm.chosen.Ep.RP.hEch.E_varOption ~= 0))) &&...
        ~ismember(condition,{'fMRI_noSatRun_choiceSplit_Elvl', 'fMRI_noSatTaskSub_noSatRun_choiceSplit_Elvl'})
    error(['condition ',condition,' is not compatible with splitting trials depending on the effort chosen and including the high effort regressor in Ep task.']);
end

% check for Em chosen period
if ((GLMprm.chosen.Em.splitPerE == 3) &&...
        ((GLMprm.chosen.Em.RP.lEch.E_varOption ~= 0) || (GLMprm.chosen.Em.RP.hEch.E_varOption ~= 0))) &&...
        ~ismember(condition,{'fMRI_noSatRun_choiceSplit_Elvl', 'fMRI_noSatTaskSub_noSatRun_choiceSplit_Elvl'})
    error(['condition ',condition,' is not compatible with splitting trials depending on the effort chosen and including the high effort regressor in Em task.']);
end
    
end % function