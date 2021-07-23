function[matlabbatch] = First_level_loadEachCondition(matlabbatch, sub_idx, iRun, iCond, cond_nm, cond_onsets, cond_dur, n_mods, mod_nm, mod_vals, orth_vars)
%[matlabbatch] = First_level_loadEachCondition(matlabbatch, sub_idx, iRun, iCond, cond_nm, cond_onsets, cond_dur, n_mods, mod_nm, mod_vals, orth_vars)
% First_level_loadEachCondition will load the data for each condition
%
% INPUTS
% matlabbatch: structure with matlab batch
%
% sub_idx: index for the current subject
%
% iRun: run number
%
% iCond: condition number
%
% cond_nm: name for the current condition
%
% cond_onsets: onsets for the current condition
% 
% cond_dur: duration of the current condition
%
% n_mods: number of modulators
%
% mod_nm: cell with name for each modulator
%
% mod_vals: n_mods*n_trials vector with values of each modulator
%
% orth_vars: orthogonalize the modulators or not (0/1)?
%
% OUTPUTS
% matlabbatch: updated structure with each condition

matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).cond(iCond).name = cond_nm;
matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).cond(iCond).onset = cond_onsets;
matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).cond(iCond).duration = cond_dur;
matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).cond(iCond).tmod = 0;

% add parametric modulators (if there are some)
% pmod need to be initialized, otherwise SPM is not happy
matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).cond(iCond).pmod = struct('name',{''},'param',{},'poly',{});
if n_mods > 0
    for iMod = 1:n_mods
        matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).cond(iCond).pmod(iMod).name = mod_nm{iMod};
        matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).cond(iCond).pmod(iMod).param = mod_vals(iMod,:);
        matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).cond(iCond).pmod(iMod).poly = 1;
    end
end

% orthogonalize regressors
switch orth_vars
    case 0
        matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).cond(iCond).orth = 0;
    case 1
        matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).cond(iCond).orth = 1;
end

end % function