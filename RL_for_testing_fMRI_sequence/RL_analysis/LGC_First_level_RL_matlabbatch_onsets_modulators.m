function [ matlabbatch ] = LGC_First_level_RL_matlabbatch_onsets_modulators( matlabbatch,...
    cond_nm, onsetCond, durCond,...
    sub_idx, iRun,...
    Nmod, task_id, modName, modulators,...
    gal)
% [ matlabbatch ] = LGC_First_level_RL_matlabbatch_onsets_modulators( matlabbatch,...
%     cond_nm, onsetCond, durCond,...
%     sub_idx, iRun,...
%     Nmod, task_id, modName, modulators,...
%     gal)
% load onsets and modulators of the different conditions inside the
% matlabbatch
%
% INPUTS:
% matlabbatch: batch for SPM before the current condition is added
%
% cond_nm: name of the current condition
% 
% onsetCond: onsets for the current condition
%
% durCond: duration for the current condition
%
% sub_idx: current subject number in the batch
%
% iRun: current run number
%
% Nmod: structure containing number of modulators for each condition
%
% task_id: task identification letter (L/G/S for learning, grip, stroop)
%
% modName: name of the eventual modulators of the current condition
%
% modulators: matrix containing a vector/eventual modulator in the same
% order as in modName
%
% gal: extracted from GLMprm (=which_GLM_MS2.m output) to know if this is a
% GLM with onsets_only (= 1 regressor/trial) or across trials and if
% variables should be orthogonalized or not
%
% OUTPUTS:
% matlabbatch: batch updated with the current condition (onset + duration +
% eventual modulators)
%

%% check if onsets_only GLM
onsets_only = gal.onsets_only;
zscore_vars = gal.zscore;

%% preparing matlab batch
n_cond = length(cond_nm);
for iCond = 1:n_cond
    
    switch onsets_only
        case 0
            curr_cond_nm    = cond_nm{iCond};
            curr_onset      = onsetCond.(curr_cond_nm);
            curr_dur        = durCond.(curr_cond_nm);
            
            % name and onset definition
            matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).cond(iCond).name = curr_cond_nm;
            matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).cond(iCond).onset = curr_onset;
            
            % stick or boxcar
            matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).cond(iCond).duration = curr_dur;
            
            matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).cond(iCond).tmod = 0;
            
            % enter parametric modulator(s)
            matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).cond(iCond).pmod = struct('name',{''},'param',{},'poly',{}); % pmod need to be initialized, otherwise SPM is not happy
            if Nmod.(curr_cond_nm) > 0
                for iMod = 1:Nmod.(curr_cond_nm)
                    curr_mod_nm = [task_id, modName.(curr_cond_nm){iMod}];
                    switch zscore_vars
                        case 0
                            curr_mod = modulators.(curr_cond_nm)(:,iMod);
                        case 1
                            curr_mod = zscore(modulators.(curr_cond_nm)(:,iMod));
                    end
                    
                    matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).cond(iCond).pmod(iMod).name = curr_mod_nm;
                    matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).cond(iCond).pmod(iMod).param = curr_mod;
                    matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).cond(iCond).pmod(iMod).poly = 1;
                end
            end
            
            % orthogonalize data
            switch gal.orth_vars
                case 0
                    matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).cond(iCond).orth = 0;
                case 1
                    matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).cond(iCond).orth = 1;
            end
            
        case 1
            curr_cond_nm    = cond_nm{iCond};
            % trial loop
            nTrials_run = length( onsetCond.(curr_cond_nm) );
            for iTrial = 1:nTrials_run
                jSample = iTrial + nTrials_run*(iCond - 1);
                curr_cond_trial_nm = [curr_cond_nm,'_sample_',num2str(iTrial)];
                
                curr_onset      = onsetCond.(curr_cond_nm)(iTrial);
                % FIR
                if gal.FIR == 1
                    TR = matlabbatch{sub_idx}.spm.stats.fmri_spec.timing.RT;
                    curr_onset = curr_onset - TR;
                end
                % duration
                if length( durCond.(curr_cond_nm) ) == 1 && durCond.(curr_cond_nm) == 0 % stick function
                    curr_dur = 0;
                else
                    curr_dur        = durCond.(curr_cond_nm)(iTrial);
                end
                
                % name and onset definition
                matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).cond(jSample).name = curr_cond_trial_nm;
                matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).cond(jSample).onset = curr_onset;
                
                % stick or boxcar
                matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).cond(jSample).duration = curr_dur;
                
                matlabbatch{sub_idx}.spm.stats.fmri_spec.sess(iRun).cond(jSample).tmod = 0;
                
                % no parametric modulator so no need for orthogonalization
                % either
                
            end % trial loop
            
    end % onsets_only
    
end % condition loop

end % function