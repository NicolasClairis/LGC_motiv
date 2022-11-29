function [] = checkGLM_and_subjectIncompatibility(study_nm, sub_nm, condition, GLMprm)
% checkGLM_and_subjectIncompatibility will check if there are some
% incompatibilities based on the GLM parameters in GLMprm and the study,
% subject and condition selected
%
% INPUTS
% study_nm: study name
%
% sub_nm: subject name
%
% condition: condition
%
% GLMprm: structure with GLM parameters

%% fixation cross pre-effort not added for the 2 first fMRI pilots
if strcmp(study_nm,'fMRI_pilots') &&...
        ismember(sub_nm,{'pilot_s1','pilot_s2'}) &&...
        (~strcmp(GLMprm.model_onset.Ep.preChoiceCross,'none') ||...
        ~strcmp(GLMprm.model_onset.Em.preChoiceCross,'none'))
    error(['pilot ',sub_nm,' did not had any fixation cross before effort. ',...
        'Please update GLM used for this participant.']);
end

%% check for subjects who saturated and GLMs including 
% variables such as confidence or value unchosen or choice (if the subject
% always picked up the same option, there is no room for these variables to
% vary => the GLM will probably not work).
switch study_nm
    case 'study1'
        switch condition
            case {'fMRI',...
                    'fMRI_noMoveSub','fMRI_noMoveSub_bis','fMRI_noMoveSub_ter',...
                    'fMRI_noMove_bis','fMRI_noMove_ter'}
               switch sub_nm
                   case {'027','047','052','095',...
                           '002','005','012','032',...
                           '048','076','100'} %% check subjects with a run fully saturated
                       periodToCheck = {'choice','chosen','Eperf'};
                       taskToCheck = {'Ep','Em'};
                       RPcond = {'R','P','RP'};
                       Econd = {'E','E1','E2','E3','Ech0','Ech1','Ech2','Ech3'};
                       problematicRegsChoice = {'choiceHighE',...
                           'money_unchosen','E_unchosen',...
                           'confidence'};
                       problematicRegsPerf = {'choiceHighE','confidence'};
                       for iP = 1:length(periodToCheck)
                           period_nm = periodToCheck{iP};
                           for iT = 1:length(taskToCheck)
                               task_nm = taskToCheck{iT};
                               for iRP = 1:length(RPcond)
                                   RP_nm = RPcond{iRP};
                                   for iE = 1:length(Econd)
                                       E_nm = Econd{iE};
                                       switch period_nm
                                           case {'choice','chosen'}
                                               for iProbRegC = 1:length(problematicRegsChoice)
                                                   reg_nm = problematicRegsChoice{iProbRegC};
                                                   prm_tmp = GLMprm.(period_nm).(task_nm).(RP_nm).(E_nm).(reg_nm);
                                                   if prm_tmp ~= 0
                                                       error(['subject ',sub_nm,' has been included under the condition ',...
                                                           condition,' while GLMprm.',period_nm,'.',task_nm,'.',...
                                                           RP_nm,'.',E_nm,'.',reg_nm,' = ',num2str(prm_tmp),'.'])
                                                   end
                                               end
                                           case 'Eperf'
                                               for iProbRegP = 1:length(problematicRegsPerf)
                                                   reg_nm = problematicRegsPerf{iProbRegP};
                                                   prm_tmp = GLMprm.(period_nm).(task_nm).(RP_nm).(E_nm).(reg_nm);
                                                   if prm_tmp ~= 0
                                                       error(['subject ',sub_nm,' has been included under the condition ',...
                                                           condition,' while GLMprm.',period_nm,'.',task_nm,'.',...
                                                           RP_nm,'.',E_nm,'.',reg_nm,' = ',num2str(prm_tmp),'.'])
                                                   end
                                               end
                                       end % period
                                   end % effort
                               end % RP
                           end % tasks
                       end % period
               end % subject
        end
end % study

end % function