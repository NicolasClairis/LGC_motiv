function [resultsFolderName, resultsFolderShortName] = fMRI_subFolder_DCM(sm_folderName, GLM, condition, DCM_mode)
% [resultsFolderName, resultsFolderShortName] = fMRI_subFolder_DCM(sm_folderName, GLM, condition, DCM_mode)
% fMRI_subFolder_DCM will determine the result folder name depending on the GLM,
% condition and DCM mode.
%
% INPUTS
% sm_folderName: path to subject results depending on the currently used
% preprocessing smoothing kernel
%
% GLM: number with info about GLM number
%
% condition: string with condition information
%
% DCM_mode:
% (1) all sessions modeled independently like in a classic univariate GLM
% => hard to manipulate for DCM but could be useful for testing
% session-specific effects or comparing sessions
% (2) sessions pooled within each task (ex: session 1 and 3 of physical
% effort will be concatenated into one single regressor) but each task will
% be modeled separately
% (3) all sessions pooled together
% (4) all trial periods are pooled together across sessions except for
% choice and effort which are modeled independently for each task (but
% pooled across sessions of the same task)
% (5) all trial periods are pooled together across sessions except for
% the effort period which is modeled independently for each task (but
% pooled across sessions of the same task)
%
% OUTPUTS
% resultsFolderName: full path to results
%
% resultsFolderShortName: short path to results
%
% See also which_DCM_mode_for_GLM.m for more details on DCM_mode
% and fMRI_subFolder.m for the classic version (not DCM) of this script

%% add a filesep at the end if not present already
if ~strcmp(sm_folderName(end),filesep)
    sm_folderName = [sm_folderName, filesep];
end
%% DCM mode label
DCM_mode_nm = ['_DCM_mode',num2str(DCM_mode)];
%% adapt name depending on condition
GLM_nm = num2str(GLM);
switch condition
    case {'fMRI','fMRI_noSatRunSub','fMRI_noSatTaskSub',...
            'fMRI_noMoveSub','fMRI_noMoveSub_bis','fMRI_noMoveSub_ter',...
            'fMRI_noSatTaskSub_noMove_bis_Sub'}
        resultsFolderShortName = ['GLM',GLM_nm,DCM_mode_nm];
    case {'fMRI_noSatTask','fMRI_noSatTask_bayesianMdl'} % saturation runs removed for the full saturated tasks
        resultsFolderShortName = ['GLM',GLM_nm,'_no_satTask',DCM_mode_nm];
    case {'fMRI_noSatRun','fMRI_noSatRun_bayesianMdl',...
            'fMRI_noSatTaskSub_noSatRun','fMRI_noSatTaskSub_noMoveSub_noSatRun'} % saturation runs removed
        resultsFolderShortName = ['GLM',GLM_nm,'_no_satRun',DCM_mode_nm];
    case {'fMRI_noSatRun_choiceSplit_Elvl',...
            'fMRI_noSatTaskSub_noSatRun_choiceSplit_Elvl'} % saturation runs removed, including those where choice is 100% correlated with effort level
        resultsFolderShortName = ['GLM',GLM_nm,'_no_satRun_Elvl',DCM_mode_nm];
    case 'fMRI_noMove_bis' % any run with movement removed (with some tolerance)
        resultsFolderShortName = ['GLM',GLM_nm,'_noMvmtRun_lenient',DCM_mode_nm];
    case 'fMRI_noMove_ter' % any run with movement removed (even slightest movement removed)
        resultsFolderShortName = ['GLM',GLM_nm,'_noMmvmtRun_stringent',DCM_mode_nm];
    case {'fMRI_noSatTask_noMove_bis','fMRI_noSatTask_noMove_bis_bayesianMdl'} % saturation runs removed for the full saturated tasks
        resultsFolderShortName = ['GLM',GLM_nm,'_no_satTask_noMmvmtRun',DCM_mode_nm];
    case 'fMRI_noSatTaskSub_noMoveSub_noSatRun_noMoveRun'
        resultsFolderShortName = ['GLM',GLM_nm,'_no_satRun_noMmvmtRun',DCM_mode_nm];
    otherwise
        error('problem with folder name');
end
resultsFolderName = [sm_folderName,...
    resultsFolderShortName,filesep];

end % function