function [] = checkGLM_and_subjectIncompatibility(study_nm, sub_nm, GLMprm)

%% fixation cross pre-effort not added for the 2 first fMRI pilots
if strcmp(study_nm,'fMRI_pilots') &&...
        ismember(sub_nm,{'pilot_s1','pilot_s2'}) &&...
        (~strcmp(GLMprm.model_onset.Ep.preChoiceCross,'none') ||...
        ~strcmp(GLMprm.model_onset.Em.preChoiceCross,'none'))
    error(['pilot ',sub_nm,' did not had any fixation cross before effort. ',...
        'Please update GLM used for this participant.']);
end

%% should add a check here for subjects who saturated and GLMs including 
% variables such as confidence or value or choice

end % function