function[matlabbatch] = First_level_loadEachCondition(matlabbatch,...
    sub_idx, iRun, iCond, cond_nm, cond_onsets, cond_dur,...
    n_mods, mod_nm, mod_vals, orth_vars, onsets_only_GLM)
%[matlabbatch] = First_level_loadEachCondition(matlabbatch,...
%   sub_idx, iRun, iCond, cond_nm, cond_onsets, cond_dur,...
%   n_mods, mod_nm, mod_vals, orth_vars, onsets_only_GLM)
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
% onsets_only_GLM: binary variable indicating whether this is an
% onsets-only GLM aimed at extracting one beta per trial or not
%
% OUTPUTS
% matlabbatch: updated structure with each condition
%

switch onsets_only_GLM
    case 0 % regular GLM
        
        matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).cond(iCond).name = cond_nm;
        matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).cond(iCond).onset = cond_onsets;
        matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).cond(iCond).duration = cond_dur;
        matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).cond(iCond).tmod = 0;
        
        %% add parametric modulators (if there are some)
        % pmod need to be initialized, otherwise SPM is not happy
        matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).cond(iCond).pmod = struct('name',{''},'param',{},'poly',{});
        if n_mods > 0
            for iMod = 1:n_mods
                matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).cond(iCond).pmod(iMod).name = mod_nm{iMod};
                matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).cond(iCond).pmod(iMod).param = mod_vals(iMod,:);
                matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).cond(iCond).pmod(iMod).poly = 1;
                
                %% check that regressor is clean before launching 1st level
                % 1) verify that the regressor has some variability (otherwise
                % SPM may crash when performing the contrast)
                % 2) verify that there are no NaNs remaining
                if all(mod_vals(iMod,:) == mod_vals(iMod,1))
                    error([mod_nm{iMod},' condition ',cond_nm,' all values are the same.']);
                elseif any(isnan(mod_vals(iMod,:)))
                    error([mod_nm{iMod},' condition ',cond_nm,' contains NaN values.']);
                end
            end
        end
        
        %% orthogonalize regressors
        switch orth_vars
            case 0
                matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).cond(iCond).orth = 0;
            case 1
                matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).cond(iCond).orth = 1;
        end
        
    case 1 % GLM with onsets_only
        
        nTrials_per_run = length( cond_onsets );
        for iTrial = 1:nTrials_per_run
            jSample = iTrial + nTrials_per_run*(iCond - 1);
            matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).cond(jSample).name = [cond_nm,'_sample_',num2str(iTrial)];
            matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).cond(jSample).onset = cond_onsets(iTrial);
            if length(cond_dur) == 1 % stick function
                matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).cond(jSample).duration = 0;
            elseif length(cond_dur) == nTrials_per_run % boxcar
                matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).cond(jSample).duration = cond_dur(iTrial);
            end
            matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).cond(jSample).tmod = 0;
            
            %% add parametric modulators (if there are some)
            % pmod need to be initialized, otherwise SPM is not happy
            matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).cond(jSample).pmod = struct('name',{''},'param',{},'poly',{});
            if n_mods > 0
                error('You should not add parametric modulators to the onsets-only GLM. Please fix that');
            end
            
            %% no regressors for onsets-only, should not be orthogonalized
            switch orth_vars
                case 0
                    matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).cond(jSample).orth = 0;
                case 1
                    error('no sense to orthogonalize non-existant regressors');
            end
        end % loop through trials included in the model
end

end % function