function [resultsFolderName, resultsFolderShortName] = fMRI_subFolder(sm_folderName, GLM, condition)
% [resultsFolderName, resultsFolderShortName] = fMRI_subFolder(sm_folderName, GLM, condition)
% fMRI_subFolder will determine the result folder name depending on the GLM
% and condition
%
% INPUTS
% sm_folderName: path to subject results depending on the currently used
% preprocessing smoothing kernel
%
% GLM: number with info about GLM number
%
% condition: string with condition information
%
% OUTPUTS
% resultsFolderName: full path to results
%
% resultsFolderShortName: short path to results

GLM_nm = num2str(GLM);
switch condition
    case {'fMRI','fMRI_noSatRunSub','fMRI_noSatTaskSub',...
            'fMRI_noMoveSub','fMRI_noMoveSub_bis','fMRI_noMoveSub_ter',...
            'fMRI_noSatTaskSub_noMove_bis_Sub'}
        resultsFolderShortName = ['GLM',GLM_nm];
    case 'fMRI_noSatTask' % saturation runs removed for the full saturated tasks
        resultsFolderShortName = ['GLM',GLM_nm,'_no_satTask'];
    case 'fMRI_noSatRun' % saturation runs removed
        resultsFolderShortName = ['GLM',GLM_nm,'_no_satRun'];
    case 'fMRI_noMove_bis' % any run with movement removed (with some tolerance)
        resultsFolderShortName = ['GLM',GLM_nm,'_noMvmtRun_lenient'];
    case 'fMRI_noMove_ter' % any run with movement removed (even slightest movement removed)
        resultsFolderShortName = ['GLM',GLM_nm,'_noMmvmtRun_stringent'];
    case 'fMRI_noSatTask_noMove_bis' % saturation runs removed for the full saturated tasks
        resultsFolderShortName = ['GLM',GLM_nm,'_no_satTask_noMmvmtRun'];
    otherwise
        resultsFolderShortName = [];
end
resultsFolderName = [sm_folderName, resultsFolderShortName,filesep];

end % function