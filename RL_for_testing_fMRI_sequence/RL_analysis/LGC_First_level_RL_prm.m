function [matlabbatch] = LGC_First_level_RL_prm(matlabbatch,...
    GLMprm, subid, sub_idx, iRun, subj_analysis_folder)
% [matlabbatch] = LGC_First_level_RL_prm(matlabbatch, GLMprm,...
%     subid, sub_idx, iRun, subj_analysis_folder)
% LGC_First_level_RL_prm loads the first level infos for
% LGC_First_level_RL_NicoC_batch.m
%
% INPUTS
% matlabbatch: batch that needs to be filled with the info of the learning
% task
%
% GLMprm: GLM parameters extracted with which_GLM_MS2 to know which
% conditions (onsets/modulators) are to be used
%
% subid: subject identification number (string)
%
% sub_idx: batch subject number for matlabbatch
%
% iRun: run to extract number
%
% subj_analysis_folder: subject analysis folder path
%
% OUTPUTS
% matlabbatch: matlabbatch for first level fMRI analysis filled with the
% infos of the current task

path_split = strsplit(subj_analysis_folder,filesep);
sub_nm = path_split{end-2};
run_nm = num2str(iRun);
gal = GLMprm.gal;
task_id = 'L_';

%% load onsets and modulators for the GLM
learn       = getfield( load([subj_analysis_folder,...
    'onsets_sub',subid,'_learning_run',run_nm],'learn'), 'learn');
onset       = learn.onset; % onsets
duration    = learn.duration; % durations
mod         = learn.mod; % modulators


%% see how many onsets need to be considered
% n_cond = 0;
cond_nm = {};
modName = struct;
modulators = struct;

%% onset stimuli
o_stim      = GLMprm.o_stim;
dur_stim    = GLMprm.dur_stim;
mod_stim    = GLMprm.mod_stim;
switch o_stim
    case 0
    case {1,2,3,5,6}
        if o_stim == 1 % onsets grouped across conditions
            curr_modStim_names1         = {'RL_stim'};
            curr_modStim_names2.RL_stim = 'allPairs';
        elseif o_stim == 2 % onsets gain pair/neutral pair/loss pair separated
            curr_modStim_names1 = {'RL_stim_gainPair', 'RL_stim_ntalPair', 'RL_stim_lossPair'};
            curr_modStim_names2.RL_stim_gainPair = 'gainPair';
            curr_modStim_names2.RL_stim_ntalPair = 'ntalPair';
            curr_modStim_names2.RL_stim_lossPair = 'lossPair';
        elseif o_stim == 3 % onsets gain+loss/neutral pair
            curr_modStim_names1 = {'RL_stim_GL_Pairs', 'RL_stim_ntalPair'};
            curr_modStim_names2.RL_stim_GL_Pairs    = 'GL_Pairs';
            curr_modStim_names2.RL_stim_ntalPair    = 'ntalPair';
        elseif o_stim == 4 % onsets first/second half
            curr_modStim_names1         = {'RL_stim_first','RL_stim_last'};
            curr_modStim_names2.RL_stim_first   = 'allPairs_first';
            curr_modStim_names2.RL_stim_last    = 'allPairs_last';
        elseif o_stim == 5 % gain/neutral/loss pair & first/second half
            curr_modStim_names1 = {'RL_stim_gainPair_first', 'RL_stim_ntalPair_first', 'RL_stim_lossPair_first',...
                'RL_stim_gainPair_last', 'RL_stim_ntalPair_last', 'RL_stim_lossPair_last'};
            curr_modStim_names2.RL_stim_gainPair_first  = 'gainPair_first';
            curr_modStim_names2.RL_stim_gainPair_last   = 'gainPair_last';
            curr_modStim_names2.RL_stim_ntalPair_first  = 'ntalPair_first';
            curr_modStim_names2.RL_stim_ntalPair_last   = 'ntalPair_last';
            curr_modStim_names2.RL_stim_lossPair_first  = 'lossPair_first';
            curr_modStim_names2.RL_stim_lossPair_last   = 'lossPair_last';
        elseif o_stim == 6 % onsets gain+loss/neutral pair first/second half
            curr_modStim_names1 = {'RL_stim_GL_Pairs_first', 'RL_stim_ntalPair_first',...
                'RL_stim_GL_Pairs_last', 'RL_stim_ntalPair_last'};
            curr_modStim_names2.RL_stim_GL_Pairs_first 	= 'GLPairs_first';
            curr_modStim_names2.RL_stim_ntalPair_first  = 'ntalPair_first';
            curr_modStim_names2.RL_stim_GL_Pairs_last   = 'GLPairs_last';
            curr_modStim_names2.RL_stim_ntalPair_last   = 'ntalPair_last';
        end
        %         n_cond = n_cond + length(curr_modStim_names1);
        cond_nm = [cond_nm, curr_modStim_names1];
        
        % loop through conditions depending on RLprm.o_stim
        for iCond_stim_nm1 = 1:length(curr_modStim_names1)
            curr_modStim_cond1 = curr_modStim_names1{iCond_stim_nm1};
            curr_modStim_cond2 = curr_modStim_names2.(curr_modStim_cond1);
            
            onsetCond.(curr_modStim_cond1) = onset.(curr_modStim_cond2).dispOptions.main;
            
            % duration
            switch dur_stim
                case 0 % stick
                    durCond.(curr_modStim_cond1) = 0;
                case 1 % boxcar from stim until choice in red
                    durCond.(curr_modStim_cond1) = duration.(curr_modStim_cond2).dispOptions.main;
                case 2 % boxcar from stim until RT
                    durCond.(curr_modStim_cond1) = mod.RT_fp.(curr_modStim_cond2).main.raw;
                case 3 % boxcar from stim until feedback
                    durCond.(curr_modStim_cond1) = duration.(curr_modStim_cond2).dispOptions.main +...
                        duration.(curr_modStim_cond2).dispChoice.main;
            end
            
            % modulators
            Nmod.(curr_modStim_cond1) = 0;
            
            %% RT
            if mod_stim.RT ~= 0 && ismember(mod_stim.RT,9:16)
                Nmod.(curr_modStim_cond1) = Nmod.(curr_modStim_cond1) + 1;
                idx_RT = Nmod.(curr_modStim_cond1);
                switch mod_stim.RT
                    case 9
                        modName.(curr_modStim_cond1)(idx_RT) = {'RT_raw'};
                        modulators.(curr_modStim_cond1)(:,idx_RT) = mod.RT_fp.(curr_modStim_cond2).main.raw;
                    case 10
                        modName.(curr_modStim_cond1)(idx_RT) = {'RT_zPerRun'};
                        modulators.(curr_modStim_cond1)(:,idx_RT) = mod.RT_fp.(curr_modStim_cond2).main.zPerRun;
                    case 11
                        error('not ready yet');
                        %                     modName.(curr_cond1)(idx_RT) = 'RT_zAcrossRuns';
                        %                     modulators.(curr_cond1)(:,idx_RT) = mod.RT_fp.(curr_cond2).main.zAcrossRuns;
                    case 12
                        error('not ready yet');
                        %                     modName.(curr_cond1)(idx_RT) = 'RT_zAcrossRuns_zPerRun';
                        %                     modulators.(curr_cond1)(:,idx_RT) = mod.RT_fp.(curr_cond2).main.raw;
                    case {13,14,15,16}
                        error('RT corrected not possible yet');
                        %                 case 5
                        %                     modName.(curr_cond1)(idx_RT) = 'RT_corr';
                        %                     modulators.(curr_cond1)(:,idx_RT) = mod.RT_fp.(curr_cond2).main.raw;
                        %                 case 6
                        %                     modName.(curr_cond1)(idx_RT) = 'RT_corr_zPerRun';
                        %                     modulators.(curr_cond1)(:,idx_RT) = mod.RT_fp.(curr_cond2).main.raw;
                        %                 case 7
                        %                     modName.(curr_cond1)(idx_RT) = 'RT_corr_zAcrossRuns';
                        %                     modulators.(curr_cond1)(:,idx_RT) = mod.RT_fp.(curr_cond2).main.raw;
                        %                 case 8
                        %                     modName.(curr_cond1)(idx_RT) = 'RT_corr_zAcrossRuns_zPerRun';
                        %                     modulators.(curr_cond1)(:,idx_RT) = mod.RT_fp.(curr_cond2).main.raw;
                    otherwise
                        error('not ready yet');
                end
            end
            
            %% trial number
            if ismember(mod_stim.trialN, [1,2,3])
                Nmod.(curr_modStim_cond1) = Nmod.(curr_modStim_cond1) + 1;
                idx_trialN = Nmod.(curr_modStim_cond1);
                switch mod_stim.trialN
                    case 1
                        modName.(curr_modStim_cond1)(idx_trialN)      = {'trialN_raw'};
                        modulators.(curr_modStim_cond1)(:,idx_trialN)   = mod.trialN.raw.(curr_modStim_cond2).main;
                    case 2
                        modName.(curr_modStim_cond1)(idx_trialN)      = {'trialN_zPerRun'};
                        modulators.(curr_modStim_cond1)(:,idx_trialN)   = mod.trialN.zPerRun.(curr_modStim_cond2).main;
                    case 3
                        modName.(curr_modStim_cond1)(idx_trialN)      = {'trialN_zAcrossRuns'};
                        modulators.(curr_modStim_cond1)(:,idx_trialN)   = mod.trialN.zAcrossRuns.(curr_modStim_cond2).main;
                end
            end
            
            %% Subjective Value based on model (before Conf)
            if o_stim == 1 ||...
                    ~ismember(curr_modStim_cond2,{'ntalPair','ntalPair_first','ntalPair_last'}) % no use of SV for neutral pair
                if ismember(mod_stim.SV,1:5)
                    mdl_type_toUse = [mod_stim.mdl_type,'_VBA_models'];
                    n_Qmodel_toUse = mod_stim.mdl_n;
                    Nmod.(curr_modStim_cond1) = Nmod.(curr_modStim_cond1) + 1;
                    idx_SV = Nmod.(curr_modStim_cond1);
                    switch mod_stim.SV
                        case 1 % pA*QA+pB*QB
                            modName.(curr_modStim_cond1)(idx_SV) = {'pAQA_pBQB'};
                            modulators.(curr_modStim_cond1)(:,idx_SV) = mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.SV.(curr_modStim_cond2).([curr_modStim_cond2,'Trials']);
                        case 2 % QA+QB
                            modName.(curr_modStim_cond1)(idx_SV) = {'QA_plus_QB'};
                            modulators.(curr_modStim_cond1)(:,idx_SV) =...
                                mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.(curr_modStim_cond2).([curr_modStim_cond2,'Trials']).bestItem +...
                                mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.(curr_modStim_cond2).([curr_modStim_cond2,'Trials']).worseItem;
                        case 3 % Qch-Qunch
                            modName.(curr_modStim_cond1)(idx_SV) = {'Qch_min_Qunch'};
                            modulators.(curr_modStim_cond1)(:,idx_SV) =...
                                mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.(curr_modStim_cond2).([curr_modStim_cond2,'Trials']).chosenItem -...
                                mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.(curr_modStim_cond2).([curr_modStim_cond2,'Trials']).unchosenItem;
                        case 4 % Qch
                            modName.(curr_modStim_cond1)(idx_SV) = {'Qch'};
                            modulators.(curr_modStim_cond1)(:,idx_SV) =...
                                mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.(curr_modStim_cond2).([curr_modStim_cond2,'Trials']).chosenItem;
                        case 5 % Qch/(QA+QB)
                            modName.(curr_modStim_cond1)(idx_SV) = {'Qch_div_sum_QA_QB'};
                            modulators.(curr_modStim_cond1)(:,idx_SV) =...
                                mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.(curr_modStim_cond2).([curr_modStim_cond2,'Trials']).chosenItem./...
                                (mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.(curr_modStim_cond2).([curr_modStim_cond2,'Trials']).bestItem +...
                                mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.(curr_modStim_cond2).([curr_modStim_cond2,'Trials']).worseItem);
                        otherwise
                            error('not ready yet');
                    end
                end
            end
            
            %% QA-QB
            if o_stim == 1 ||...
                    ~ismember(curr_modStim_cond2,{'ntalPair','ntalPair_first','ntalPair_last'}) % no use of QA-QB for neutral pair
                if mod_stim.dQ ~= 0
                    mdl_type_toUse = [mod_stim.mdl_type,'_VBA_models'];
                    n_Qmodel_toUse = mod_stim.mdl_n;
                    Nmod.(curr_modStim_cond1) = Nmod.(curr_modStim_cond1) + 1;
                    idx_dQ = Nmod.(curr_modStim_cond1);
                    switch mod_stim.dQ
                        case {1,2}
                            modName.(curr_modStim_cond1)(idx_dQ) = {'dQ'};
                            dQ_values = mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.dQ.(curr_modStim_cond2).([curr_modStim_cond2,'Trials']);
                            if ismember(mod_stim.dQ,[2,8])
                                dQ_values = zscore(dQ_values);
                            end
                            modulators.(curr_modStim_cond1)(:,idx_dQ) = dQ_values;
                        case {3,4,5}
                            error('dQ zscore not ready yet (think about how to zscore)');
                        case {6}
                            modName.(curr_modStim_cond1)(idx_dQ) = {'dQsquared'};
                            modulators.(curr_modStim_cond1)(:,idx_dQ) = ( mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.dQ.(curr_modStim_cond2).([curr_modStim_cond2,'Trials']) ).^2;
                        case {7,8}
                            modName.(curr_modStim_cond1)(idx_dQ) = {'|dQ|/sqrt(QbestSquared+QnonBestSquared)'};
                            Qbest = mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.(curr_modStim_cond2).([curr_modStim_cond2,'Trials']).bestItem;
                            Qworse = mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.(curr_modStim_cond2).([curr_modStim_cond2,'Trials']).worseItem;
                            dQ_values = abs( Qbest - Qworse )./sqrt( (Qbest + 0.5).^2 + (Qworse + 0.5).^2 ); % add 0.5 to avoid problem with dividing by zero
                            if ismember(mod_stim.dQ,[15,16]) % zscore/run
                                    dQ_values = zscore(dQ_values);
                            end
                            modulators.(curr_modStim_cond1)(:,idx_dQ) = dQ_values;
                        otherwise
                            error('not ready yet');
                    end
                end
            end
            
            %% p(choice=best option)
            if o_stim == 1 ||...
                    ~ismember(curr_modStim_cond2,{'ntalPair','ntalPair_first','ntalPair_last'}) % no use of QA-QB for neutral pair
                if mod_stim.pBest ~= 0
                    mdl_type_toUse = [mod_stim.mdl_type,'_VBA_models'];
                    n_Qmodel_toUse = mod_stim.mdl_n;
                    Nmod.(curr_modStim_cond1) = Nmod.(curr_modStim_cond1) + 1;
                    idx_pA = Nmod.(curr_modStim_cond1);
                    switch mod_stim.pBest
                        case 1
                            modName.(curr_modStim_cond1)(idx_pA) = {'pBest_centered_squared'};
                            modulators.(curr_modStim_cond1)(:,idx_pA) = ( mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.pChoice.best.(curr_modStim_cond2).([curr_modStim_cond2,'Trials']) - 0.5).^2;
                        case 2
                            modName.(curr_modStim_cond1)(idx_pA) = {'pLeft_centered_squared'};
                            modulators.(curr_modStim_cond1)(:,idx_pA) = ( mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.pChoice.left.(curr_modStim_cond2).([curr_modStim_cond2,'Trials']) - 0.5).^2;
                        case 3
                            modName.(curr_modStim_cond1)(idx_pA) = {'pChosen'};
                            modulators.(curr_modStim_cond1)(:,idx_pA) = mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.pChoice.chosen.(curr_modStim_cond2).([curr_modStim_cond2,'Trials']);
                        otherwise
                            error('not ready yet');
                    end
                end
            end
            
            %% trial number
            if ismember(mod_stim.trialN, [4,5,6])
                Nmod.(curr_modStim_cond1) = Nmod.(curr_modStim_cond1) + 1;
                idx_trialN = Nmod.(curr_modStim_cond1);
                switch mod_stim.trialN
                    case 4
                        modName.(curr_modStim_cond1)(idx_trialN)      = {'trialN_raw'};
                        modulators.(curr_modStim_cond1)(:,idx_trialN)   = mod.trialN.raw.(curr_modStim_cond2).main;
                    case 5
                        modName.(curr_modStim_cond1)(idx_trialN)      = {'trialN_zPerRun'};
                        modulators.(curr_modStim_cond1)(:,idx_trialN)   = mod.trialN.zPerRun.(curr_modStim_cond2).main;
                    case 6
                        modName.(curr_modStim_cond1)(idx_trialN)      = {'trialN_zAcrossRuns'};
                        modulators.(curr_modStim_cond1)(:,idx_trialN)   = mod.trialN.zAcrossRuns.(curr_modStim_cond2).main;
                end
            end
            
            %% ROI activity
            if mod_stim.ROI_activity_yn ~= 0
                Nmod.(curr_modStim_cond1) = Nmod.(curr_modStim_cond1) + 1;
                idx_oStim_ROI = Nmod.(curr_modStim_cond1);
                modName.(curr_modStim_cond1)(idx_oStim_ROI) = {[]};
                % get the data
                trialN_idx = mod.trialN.raw.(curr_modStim_cond2).main; % trial index
                ROI_nm          = mod_stim.ROI_activity_ROI_nm; % ROI name
                GLM_ROI_nm      = mod_stim.ROI_activity_GLM;
                ROI_period_nm   = mod_stim.ROI_activity_period; % trial period when ROI extracted
                jRun_nb = (iRun == 1)*1 + (iRun == 4)*2 + (iRun == 7)*3;
                jRun_nm = ['run',num2str(jRun_nb)];
                ROI_data = getfield(load([subj_analysis_folder,'ROI',filesep,'GLM',GLM_ROI_nm,'_',ROI_nm,'.mat'],...
                    'sub_ROI_trial_b_trial'),'sub_ROI_trial_b_trial');
                modulators.(curr_modStim_cond1)(:,idx_oStim_ROI) = ROI_data.(ROI_nm).L.(ROI_period_nm).(jRun_nm).(sub_nm)(trialN_idx);
            end
            
            %% RT
            if mod_stim.RT ~= 0 && ismember(mod_stim.RT,1:8)
                Nmod.(curr_modStim_cond1) = Nmod.(curr_modStim_cond1) + 1;
                idx_RT = Nmod.(curr_modStim_cond1);
                switch mod_stim.RT
                    case 1
                        modName.(curr_modStim_cond1)(idx_RT) = {'RT_raw'};
                        modulators.(curr_modStim_cond1)(:,idx_RT) = mod.RT_fp.(curr_modStim_cond2).main.raw;
                    case 2
                        modName.(curr_modStim_cond1)(idx_RT) = {'RT_zPerRun'};
                        modulators.(curr_modStim_cond1)(:,idx_RT) = mod.RT_fp.(curr_modStim_cond2).main.zPerRun;
                    case 3
                        error('not ready yet');
                        %                     modName.(curr_cond1)(idx_RT) = 'RT_zAcrossRuns';
                        %                     modulators.(curr_cond1)(:,idx_RT) = mod.RT_fp.(curr_cond2).main.zAcrossRuns;
                    case 4
                        error('not ready yet');
                        %                     modName.(curr_cond1)(idx_RT) = 'RT_zAcrossRuns_zPerRun';
                        %                     modulators.(curr_cond1)(:,idx_RT) = mod.RT_fp.(curr_cond2).main.raw;
                    case {5,6,7,8}
                        error('RT corrected not possible yet');
                        %                 case 5
                        %                     modName.(curr_cond1)(idx_RT) = 'RT_corr';
                        %                     modulators.(curr_cond1)(:,idx_RT) = mod.RT_fp.(curr_cond2).main.raw;
                        %                 case 6
                        %                     modName.(curr_cond1)(idx_RT) = 'RT_corr_zPerRun';
                        %                     modulators.(curr_cond1)(:,idx_RT) = mod.RT_fp.(curr_cond2).main.raw;
                        %                 case 7
                        %                     modName.(curr_cond1)(idx_RT) = 'RT_corr_zAcrossRuns';
                        %                     modulators.(curr_cond1)(:,idx_RT) = mod.RT_fp.(curr_cond2).main.raw;
                        %                 case 8
                        %                     modName.(curr_cond1)(idx_RT) = 'RT_corr_zAcrossRuns_zPerRun';
                        %                     modulators.(curr_cond1)(:,idx_RT) = mod.RT_fp.(curr_cond2).main.raw;
                    otherwise
                        error('not ready yet');
                end
            end
            
            %% Subjective Value based on model (after Conf and DT)
            if o_stim == 1 ||...
                    ~ismember(curr_modStim_cond2,{'ntalPair','ntalPair_first','ntalPair_last'}) % no use of SV for neutral pair
                if ismember(mod_stim.SV,6:10)
                    mdl_type_toUse = [mod_stim.mdl_type,'_VBA_models'];
                    n_Qmodel_toUse = mod_stim.mdl_n;
                    Nmod.(curr_modStim_cond1) = Nmod.(curr_modStim_cond1) + 1;
                    idx_SV = Nmod.(curr_modStim_cond1);
                    switch mod_stim.SV
                        case 6 % pA*QA+pB*QB
                            modName.(curr_modStim_cond1)(idx_SV) = {'pAQA_pBQB'};
                            modulators.(curr_modStim_cond1)(:,idx_SV) = mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.SV.(curr_modStim_cond2).([curr_modStim_cond2,'Trials']);
                        case 7 % QA+QB
                            modName.(curr_modStim_cond1)(idx_SV) = {'QA_plus_QB'};
                            modulators.(curr_modStim_cond1)(:,idx_SV) =...
                                mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.(curr_modStim_cond2).([curr_modStim_cond2,'Trials']).bestItem +...
                                mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.(curr_modStim_cond2).([curr_modStim_cond2,'Trials']).worseItem;
                        case 8 % Qch-Qunch
                            modName.(curr_modStim_cond1)(idx_SV) = {'Qch_min_Qunch'};
                            modulators.(curr_modStim_cond1)(:,idx_SV) =...
                                mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.(curr_modStim_cond2).([curr_modStim_cond2,'Trials']).chosenItem -...
                                mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.(curr_modStim_cond2).([curr_modStim_cond2,'Trials']).unchosenItem;
                        case 9 % Qch
                            modName.(curr_modStim_cond1)(idx_SV) = {'Qch'};
                            modulators.(curr_modStim_cond1)(:,idx_SV) =...
                                mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.(curr_modStim_cond2).([curr_modStim_cond2,'Trials']).chosenItem;
                        case 10 % Qch/(QA+QB)
                            modName.(curr_modStim_cond1)(idx_SV) = {'Qch_div_sum_QA_QB'};
                            modulators.(curr_modStim_cond1)(:,idx_SV) =...
                                mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.(curr_modStim_cond2).([curr_modStim_cond2,'Trials']).chosenItem./...
                                (mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.(curr_modStim_cond2).([curr_modStim_cond2,'Trials']).bestItem +...
                                mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.(curr_modStim_cond2).([curr_modStim_cond2,'Trials']).worseItem);
                        otherwise
                            error('not ready yet');
                    end
                end
            end
            
            %% trial number
            if ismember(mod_stim.trialN, [7,8,9])
                Nmod.(curr_modStim_cond1) = Nmod.(curr_modStim_cond1) + 1;
                idx_trialN = Nmod.(curr_modStim_cond1);
                switch mod_stim.trialN
                    case 7
                        modName.(curr_modStim_cond1)(idx_trialN)      = {'trialN_raw'};
                        modulators.(curr_modStim_cond1)(:,idx_trialN)   = mod.trialN.raw.(curr_modStim_cond2).main;
                    case 8
                        modName.(curr_modStim_cond1)(idx_trialN)      = {'trialN_zPerRun'};
                        modulators.(curr_modStim_cond1)(:,idx_trialN)   = mod.trialN.zPerRun.(curr_modStim_cond2).main;
                    case 9
                        modName.(curr_modStim_cond1)(idx_trialN)      = {'trialN_zAcrossRuns'};
                        modulators.(curr_modStim_cond1)(:,idx_trialN)   = mod.trialN.zAcrossRuns.(curr_modStim_cond2).main;
                end
            end
            
            
        end % condition loop gain/neutral/loss pair
otherwise
    error('condition not ready yet');
end % onset_stim filter

%% onset choice answer (moment subject answers)
o_answer = GLMprm.o_answer;
dur_answer = GLMprm.dur_answer;
% mod_answer = RLprm.mod_answer;
switch o_answer
    case 0
    case {1,2}
        if o_answer == 1 % onsets grouped across conditions
            curr_modChoice_names1 = {'RL_choice'};
            curr_modChoice_names2.RL_choice = 'allPairs';
        elseif o_answer == 2 % onsets gain pair/neutral pair/loss pair separated
            curr_modChoice_names1 = {'RL_choice_gainPair','RL_choice_ntalPair','RL_choice_lossPair'};
            curr_modChoice_names2.RL_choice_gainPair = 'gainPair';
            curr_modChoice_names2.RL_choice_ntalPair = 'ntalPair';
            curr_modChoice_names2.RL_choice_lossPair = 'lossPair';
        end
        %         n_cond = n_cond + length(curr_modChoice_names1);
        cond_nm = [cond_nm, curr_modChoice_names1];
        
        % loop through conditions depending on RLprm.o_stim
        for iCond_answer_nm1 = 1:length(curr_modChoice_names1)
            curr_modAnswer_cond1 = curr_modChoice_names1{iCond_answer_nm1};
            curr_modAnswer_cond2 = curr_modChoice_names2.(curr_modAnswer_cond1);
            
            onsetCond.(curr_modAnswer_cond1) = onset.(curr_modAnswer_cond2).choice.main;
            
            %% duration
            switch dur_answer
                case 0 % stick
                    durCond.(curr_modStim_cond1) = 0;
                case 1 % boxcar from RT until choice in red
                    durCond.(curr_modStim_cond1) = duration.(curr_modStim_cond2).choice.main;
                case 2 % boxcar from RT until feedback
                    durCond.(curr_modStim_cond1) = duration.(curr_modStim_cond2).choice.main +...
                        duration.(curr_modStim_cond2).dispChoice.main;
            end
            
            %             % modulators
            %             Nmod.(curr_modAnswer_cond1) = 0;
            %
            %             % trial number
            %             if mod_answer.trialN ~= 0
            %                 Nmod.(curr_modAnswer_cond1) = Nmod.(curr_modAnswer_cond1) + 1;
            %                 idx_trialN = Nmod.(curr_modAnswer_cond1);
            %                 switch mod_answer.trialN
            %                     case 1
            %                         modName.(curr_modAnswer_cond1)(idx_trialN)      = {'trialN_raw'};
            %                         modulators.(curr_modAnswer_cond1)(:,idx_trialN)   = mod.trialN.raw.(curr_modAnswer_cond2).main;
            %                     case 2
            %                         modName.(curr_modAnswer_cond1)(idx_trialN)      = {'trialN_zPerRun'};
            %                         modulators.(curr_modAnswer_cond1)(:,idx_trialN)   = mod.trialN.zPerRun.(curr_modAnswer_cond2).main;
            %                     case 3
            %                         modName.(curr_modAnswer_cond1)(idx_trialN)      = {'trialN_zAcrossRuns'};
            %                         modulators.(curr_modAnswer_cond1)(:,idx_trialN)   = mod.trialN.zAcrossRuns.(curr_modAnswer_cond2).main;
            %                 end
            %             end
            %
            %             % QA-QB
            %             if o_answer == 1 || ~strcmp(curr_modAnswer_cond2,'ntalPair') % no use of QA-QB for ntal pair
            %                 if mod_answer.dQ ~= 0
            %                     Nmod.(curr_modAnswer_cond1) = Nmod.(curr_modAnswer_cond1) + 1;
            %                     idx_RT = Nmod.(curr_modAnswer_cond1);
            %                     switch mod_answer.dQ
            %                         case 1
            %                             modName.(curr_modAnswer_cond1)(idx_RT) = {'dQ'};
            %                             modulators.(curr_modAnswer_cond1)(:,idx_RT) = mod.(mdl_type_toUse).Q_model(3).raw.dQ.(curr_modAnswer_cond2).([curr_modAnswer_cond2,'Trials']);
            %                         case {2,3,4}
            %                             error('dQ zscore not ready yet (think about how to zscore)');
            %                     end
            %                 end
            %             end
            %
            %
            %             % RT
            %             if mod_answer.RT ~= 0
            %                 Nmod.(curr_modAnswer_cond1) = Nmod.(curr_modAnswer_cond1) + 1;
            %                 idx_RT = Nmod.(curr_modAnswer_cond1);
            %                 switch mod_answer.RT
            %                     case 1
            %                         modName.(curr_modAnswer_cond1)(idx_RT) = {'RT_raw'};
            %                         modulators.(curr_modAnswer_cond1)(:,idx_RT) = mod.RT_fp.(curr_modAnswer_cond2).main.raw;
            %                     case 2
            %                         modName.(curr_modAnswer_cond1)(idx_RT) = {'RT_zPerRun'};
            %                         modulators.(curr_modAnswer_cond1)(:,idx_RT) = mod.RT_fp.(curr_modAnswer_cond2).main.zPerRun;
            %                     case 3
            %                         error('not ready yet');
            %                         %                     modName.(curr_cond1)(idx_RT) = 'RT_zAcrossRuns';
            %                         %                     modulators.(curr_cond1)(idx_RT) = mod.RT_fp.(curr_cond2).main.zAcrossRuns;
            %                     case 4
            %                         error('not ready yet');
            %                         %                     modName.(curr_cond1)(idx_RT) = 'RT_zAcrossRuns_zPerRun';
            %                         %                     modulators.(curr_cond1)(idx_RT) = mod.RT_fp.(curr_cond2).main.raw;
            %                     case {5,6,7,8}
            %                         error('RT corrected not possible yet');
            %                         %                 case 5
            %                         %                     modName.(curr_cond1)(idx_RT) = 'RT_corr';
            %                         %                     modulators.(curr_cond1)(idx_RT) = mod.RT_fp.(curr_cond2).main.raw;
            %                         %                 case 6
            %                         %                     modName.(curr_cond1)(idx_RT) = 'RT_corr_zPerRun';
            %                         %                     modulators.(curr_cond1)(idx_RT) = mod.RT_fp.(curr_cond2).main.raw;
            %                         %                 case 7
            %                         %                     modName.(curr_cond1)(idx_RT) = 'RT_corr_zAcrossRuns';
            %                         %                     modulators.(curr_cond1)(idx_RT) = mod.RT_fp.(curr_cond2).main.raw;
            %                         %                 case 8
            %                         %                     modName.(curr_cond1)(idx_RT) = 'RT_corr_zAcrossRuns_zPerRun';
            %                         %                     modulators.(curr_cond1)(idx_RT) = mod.RT_fp.(curr_cond2).main.raw;
            %                 end
            %             end
            
            
        end % condition loop gain/ntal/loss pair
otherwise
    error('condition not ready yet');
end % onset_answer filter

%% onset chosen option circled in red
o_chosen    = GLMprm.o_chosen;
dur_chosen  = GLMprm.dur_chosen;
mod_chosen  = GLMprm.mod_chosen;
switch o_chosen
    case 0
    case {1,2,3,5}
        if o_chosen == 1 % onsets grouped across conditions
            curr_modChosen_names1         = {'RL_chosen'};
            curr_modChosen_names2.RL_chosen = 'allPairs';
        elseif o_chosen == 2 % onsets gain pair/neutral pair/loss pair separated
            curr_modChosen_names1 = {'RL_chosen_gainPair', 'RL_chosen_ntalPair', 'RL_chosen_lossPair'};
            curr_modChosen_names2.RL_chosen_gainPair = 'gainPair';
            curr_modChosen_names2.RL_chosen_ntalPair = 'ntalPair';
            curr_modChosen_names2.RL_chosen_lossPair = 'lossPair';
        elseif o_chosen == 3 % onsets gain+loss/neutral pair
            curr_modChosen_names1 = {'RL_chosen_GL_Pairs', 'RL_chosen_ntalPair'};
            curr_modChosen_names2.RL_chosen_GL_Pairs    = 'GL_Pairs';
            curr_modChosen_names2.RL_chosen_ntalPair    = 'ntalPair';
        elseif o_chosen == 4 % onsets first/second half
            curr_modChosen_names1         = {'RL_chosen_first','RL_chosen_last'};
            curr_modChosen_names2.RL_chosen_first   = 'allPairs_first';
            curr_modChosen_names2.RL_chosen_last    = 'allPairs_last';
        elseif o_chosen == 5 % gain/neutral/loss pair & first/second half
            curr_modChosen_names1 = {'RL_chosen_gainPair_first', 'RL_chosen_ntalPair_first', 'RL_chosen_lossPair_first',...
                'RL_chosen_gainPair_last', 'RL_chosen_ntalPair_last', 'RL_chosen_lossPair_last'};
            curr_modChosen_names2.RL_chosen_gainPair_first  = 'gainPair_first';
            curr_modChosen_names2.RL_chosen_gainPair_last   = 'gainPair_last';
            curr_modChosen_names2.RL_chosen_ntalPair_first  = 'ntalPair_first';
            curr_modChosen_names2.RL_chosen_ntalPair_last   = 'ntalPair_last';
            curr_modChosen_names2.RL_chosen_lossPair_first  = 'lossPair_first';
            curr_modChosen_names2.RL_chosen_lossPair_last   = 'lossPair_last';
        end
        %         n_cond = n_cond + length(curr_modChosen_names1);
        cond_nm = [cond_nm, curr_modChosen_names1];
        
        % loop through conditions depending on RLprm.o_chosen
        for iCond_chosen_nm1 = 1:length(curr_modChosen_names1)
            curr_modChosen_cond1 = curr_modChosen_names1{iCond_chosen_nm1};
            curr_modChosen_cond2 = curr_modChosen_names2.(curr_modChosen_cond1);
            
            onsetCond.(curr_modChosen_cond1) = onset.(curr_modChosen_cond2).dispChoice.main;
            
            %% duration
            switch dur_chosen
                case 0 % stick
                    durCond.(curr_modChosen_cond1) = 0;
                case 1 % boxcar from choice in red until feedback
                    durCond.(curr_modChosen_cond1) = duration.(curr_modChosen_cond2).dispChoice.main;
            end
            
            % modulators
            Nmod.(curr_modChosen_cond1) = 0;
            
            %% trial number
            if ismember(mod_chosen.trialN, [1,2,3])
                Nmod.(curr_modChosen_cond1) = Nmod.(curr_modChosen_cond1) + 1;
                idx_trialN = Nmod.(curr_modChosen_cond1);
                switch mod_chosen.trialN
                    case 1
                        modName.(curr_modChosen_cond1)(idx_trialN)      = {'trialN_raw'};
                        modulators.(curr_modChosen_cond1)(:,idx_trialN)   = mod.trialN.raw.(curr_modChosen_cond2).main;
                    case 2
                        modName.(curr_modChosen_cond1)(idx_trialN)      = {'trialN_zPerRun'};
                        modulators.(curr_modChosen_cond1)(:,idx_trialN)   = mod.trialN.zPerRun.(curr_modChosen_cond2).main;
                    case 3
                        modName.(curr_modChosen_cond1)(idx_trialN)      = {'trialN_zAcrossRuns'};
                        modulators.(curr_modChosen_cond1)(:,idx_trialN)   = mod.trialN.zAcrossRuns.(curr_modChosen_cond2).main;
                end
            end
            
            %% SV: pA*QA+pB*QB (or QA+QB)
            if o_chosen == 1 ||...
                    ~ismember(curr_modChosen_cond2,{'ntalPair','ntalPair_first','ntalPair_last'}) % no use of SV for neutral pair
                if mod_chosen.SV ~= 0
                    mdl_type_toUse = [mod_chosen.mdl_type,'_VBA_models'];
                    n_Qmodel_toUse = mod_chosen.mdl_n;
                    Nmod.(curr_modChosen_cond1) = Nmod.(curr_modChosen_cond1) + 1;
                    idx_SV = Nmod.(curr_modChosen_cond1);
                    switch mod_chosen.SV
                        case 1 % pA*QA+pB*QB
                            modName.(curr_modChosen_cond1)(idx_SV) = {'pAQA_pBQB'};
                            modulators.(curr_modChosen_cond1)(:,idx_SV) = mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.SV.(curr_modChosen_cond2).([curr_modChosen_cond2,'Trials']);
                        case 2 % QA+QB
                            modName.(curr_modChosen_cond1)(idx_SV) = {'QA_plus_QB'};
                            modulators.(curr_modChosen_cond1)(:,idx_SV) =...
                                mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.(curr_modChosen_cond2).([curr_modChosen_cond2,'Trials']).bestItem +...
                                mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.(curr_modChosen_cond2).([curr_modChosen_cond2,'Trials']).worseItem;
                        case 3 % Qch-Qunch
                            modName.(curr_modChosen_cond1)(idx_SV) = {'Qch_min_Qunch'};
                            modulators.(curr_modChosen_cond1)(:,idx_SV) =...
                                mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.(curr_modChosen_cond2).([curr_modChosen_cond2,'Trials']).chosenItem -...
                                mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.(curr_modChosen_cond2).([curr_modChosen_cond2,'Trials']).unchosenItem;
                        case 4 % Qch
                            modName.(curr_modChosen_cond1)(idx_SV) = {'Qch'};
                            modulators.(curr_modChosen_cond1)(:,idx_SV) =...
                                mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.(curr_modChosen_cond2).([curr_modChosen_cond2,'Trials']).chosenItem;
                        case 5 % Qch/(QA+QB)
                            modName.(curr_modChosen_cond1)(idx_SV) = {'Qch_div_sum_QA_QB'};
                            modulators.(curr_modChosen_cond1)(:,idx_SV) =...
                                mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.(curr_modChosen_cond2).([curr_modChosen_cond2,'Trials']).chosenItem./...
                                (mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.(curr_modChosen_cond2).([curr_modChosen_cond2,'Trials']).bestItem +...
                                mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.(curr_modChosen_cond2).([curr_modChosen_cond2,'Trials']).worseItem);
                        otherwise
                            error('not ready yet');
                    end
                end
            end
            
            %% QA-QB
            if o_chosen == 1 ||...
                    ~ismember(curr_modChosen_cond2,{'ntalPair','ntalPair_first','ntalPair_last'}) % no use of QA-QB for neutral pair
                if mod_chosen.dQ ~= 0
                    mdl_type_toUse = [mod_chosen.mdl_type,'_VBA_models'];
                    n_Qmodel_toUse = mod_chosen.mdl_n;
                    Nmod.(curr_modChosen_cond1) = Nmod.(curr_modChosen_cond1) + 1;
                    idx_dQ = Nmod.(curr_modChosen_cond1);
                    switch mod_chosen.dQ
                        case {1,2}
                            modName.(curr_modChosen_cond1)(idx_dQ) = {'dQ'};
                            dQ_values = mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.dQ.(curr_modChosen_cond2).([curr_modChosen_cond2,'Trials']);
                            if ismember(mod_chosen.dQ,[2,8])
                                dQ_values = zscore(dQ_values);
                            end
                            modulators.(curr_modChosen_cond1)(:,idx_dQ) = dQ_values;
                        case {3,4,5}
                            error('dQ zscore not ready yet (think about how to zscore)');
                        case {6}
                            modName.(curr_modChosen_cond1)(idx_dQ) = {'dQsquared'};
                            modulators.(curr_modChosen_cond1)(:,idx_dQ) = ( mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.dQ.(curr_modChosen_cond2).([curr_modChosen_cond2,'Trials']) ).^2;
                        case {7,8}
                            modName.(curr_modChosen_cond1)(idx_dQ) = {'|dQ|/sqrt(QbestSquared+QnonBestSquared)'};
                            Qbest = mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.(curr_modChosen_cond2).([curr_modChosen_cond2,'Trials']).bestItem;
                            Qworse = mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.(curr_modChosen_cond2).([curr_modChosen_cond2,'Trials']).worseItem;
                            dQ_values = abs( Qbest - Qworse )./sqrt( (Qbest + 0.5).^2 + (Qworse + 0.5).^2 ); % add 0.5 to avoid problem with dividing by zero
                            if ismember(mod_chosen.dQ,[15,16]) % zscore/run
                                    dQ_values = zscore(dQ_values);
                            end
                            modulators.(curr_modChosen_cond1)(:,idx_dQ) = dQ_values;
                        otherwise
                            error('not ready yet');
                    end
                end
            end
            
            %% p(choice=best option)
            if o_chosen == 1 ||...
                    ~ismember(curr_modChosen_cond2,{'ntalPair','ntalPair_first','ntalPair_last'}) % no use of QA-QB for neutral pair
                if mod_chosen.pBest ~= 0
                    mdl_type_toUse = [mod_chosen.mdl_type,'_VBA_models'];
                    n_Qmodel_toUse = mod_chosen.mdl_n;
                    Nmod.(curr_modChosen_cond1) = Nmod.(curr_modChosen_cond1) + 1;
                    idx_pA = Nmod.(curr_modChosen_cond1);
                    switch mod_chosen.pBest
                        case 1
                            modName.(curr_modChosen_cond1)(idx_pA) = {'pBest_centered_squared'};
                            modulators.(curr_modChosen_cond1)(:,idx_pA) = ( mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.pChoice.best.(curr_modChosen_cond2).([curr_modChosen_cond2,'Trials']) - 0.5).^2;
                        case 2
                            modName.(curr_modChosen_cond1)(idx_pA) = {'pLeft_centered_squared'};
                            modulators.(curr_modChosen_cond1)(:,idx_pA) = ( mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.pChoice.left.(curr_modChosen_cond2).([curr_modChosen_cond2,'Trials']) - 0.5).^2;
                        case 3
                            modName.(curr_modChosen_cond1)(idx_pA) = {'pChosen'};
                            modulators.(curr_modChosen_cond1)(:,idx_pA) = mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.pChoice.chosen.(curr_modChosen_cond2).([curr_modChosen_cond2,'Trials']);
                        otherwise
                            error('not ready yet');
                    end
                end
            end
            
            %% trial number
            if ismember(mod_chosen.trialN, [4,5,6])
                Nmod.(curr_modChosen_cond1) = Nmod.(curr_modChosen_cond1) + 1;
                idx_trialN = Nmod.(curr_modChosen_cond1);
                switch mod_chosen.trialN
                    case 4
                        modName.(curr_modChosen_cond1)(idx_trialN)      = {'trialN_raw'};
                        modulators.(curr_modChosen_cond1)(:,idx_trialN)   = mod.trialN.raw.(curr_modChosen_cond2).main;
                    case 5
                        modName.(curr_modChosen_cond1)(idx_trialN)      = {'trialN_zPerRun'};
                        modulators.(curr_modChosen_cond1)(:,idx_trialN)   = mod.trialN.zPerRun.(curr_modChosen_cond2).main;
                    case 6
                        modName.(curr_modChosen_cond1)(idx_trialN)      = {'trialN_zAcrossRuns'};
                        modulators.(curr_modChosen_cond1)(:,idx_trialN)   = mod.trialN.zAcrossRuns.(curr_modChosen_cond2).main;
                end
            end
            
            %% ROI activity
            if mod_chosen.ROI_activity_yn ~= 0
                Nmod.(curr_modChosen_cond1) = Nmod.(curr_modChosen_cond1) + 1;
                idx_oChosen_ROI = Nmod.(curr_modChosen_cond1);
                modName.(curr_modChosen_cond1)(idx_oChosen_ROI) = {[]};
                % get the data
                trialN_idx = mod.trialN.raw.(curr_modChosen_cond2).main; % trial index
                ROI_nm          = mod_chosen.ROI_activity_ROI_nm; % ROI name
                GLM_ROI_nm      = mod_chosen.ROI_activity_GLM;
                ROI_period_nm   = mod_chosen.ROI_activity_period; % trial period when ROI extracted
                jRun_nb = (iRun == 1)*1 + (iRun == 4)*2 + (iRun == 7)*3;
                jRun_nm = ['run',num2str(jRun_nb)];
                ROI_data = getfield(load([subj_analysis_folder,'ROI',filesep,'GLM',GLM_ROI_nm,'_',ROI_nm,'.mat'],...
                    'sub_ROI_trial_b_trial'),'sub_ROI_trial_b_trial');
                modulators.(curr_modChosen_cond1)(:,idx_oChosen_ROI) = ROI_data.(ROI_nm).L.(ROI_period_nm).(jRun_nm).(sub_nm)(trialN_idx);
            end
            
            %% RT
            if mod_chosen.RT ~= 0
                Nmod.(curr_modChosen_cond1) = Nmod.(curr_modChosen_cond1) + 1;
                idx_RT = Nmod.(curr_modChosen_cond1);
                switch mod_chosen.RT
                    case 1
                        modName.(curr_modChosen_cond1)(idx_RT) = {'RT_raw'};
                        modulators.(curr_modChosen_cond1)(:,idx_RT) = mod.RT_fp.(curr_modChosen_cond2).main.raw;
                    case 2
                        modName.(curr_modChosen_cond1)(idx_RT) = {'RT_zPerRun'};
                        modulators.(curr_modChosen_cond1)(:,idx_RT) = mod.RT_fp.(curr_modChosen_cond2).main.zPerRun;
                    case 3
                        error('not ready yet');
                        %                     modName.(curr_cond1)(idx_RT) = 'RT_zAcrossRuns';
                        %                     modulators.(curr_cond1)(:,idx_RT) = mod.RT_fp.(curr_cond2).main.zAcrossRuns;
                    case 4
                        error('not ready yet');
                        %                     modName.(curr_cond1)(idx_RT) = 'RT_zAcrossRuns_zPerRun';
                        %                     modulators.(curr_cond1)(:,idx_RT) = mod.RT_fp.(curr_cond2).main.raw;
                    case {5,6,7,8}
                        error('RT corrected not possible yet');
                        %                 case 5
                        %                     modName.(curr_cond1)(idx_RT) = 'RT_corr';
                        %                     modulators.(curr_cond1)(:,idx_RT) = mod.RT_fp.(curr_cond2).main.raw;
                        %                 case 6
                        %                     modName.(curr_cond1)(idx_RT) = 'RT_corr_zPerRun';
                        %                     modulators.(curr_cond1)(:,idx_RT) = mod.RT_fp.(curr_cond2).main.raw;
                        %                 case 7
                        %                     modName.(curr_cond1)(idx_RT) = 'RT_corr_zAcrossRuns';
                        %                     modulators.(curr_cond1)(:,idx_RT) = mod.RT_fp.(curr_cond2).main.raw;
                        %                 case 8
                        %                     modName.(curr_cond1)(idx_RT) = 'RT_corr_zAcrossRuns_zPerRun';
                        %                     modulators.(curr_cond1)(:,idx_RT) = mod.RT_fp.(curr_cond2).main.raw;
                    otherwise
                        error('not ready yet');
                end
            end
            
            % trial number
            if ismember(mod_chosen.trialN, [7,8,9])
                Nmod.(curr_modChosen_cond1) = Nmod.(curr_modChosen_cond1) + 1;
                idx_trialN = Nmod.(curr_modChosen_cond1);
                switch mod_chosen.trialN
                    case 7
                        modName.(curr_modChosen_cond1)(idx_trialN)      = {'trialN_raw'};
                        modulators.(curr_modChosen_cond1)(:,idx_trialN)   = mod.trialN.raw.(curr_modChosen_cond2).main;
                    case 8
                        modName.(curr_modChosen_cond1)(idx_trialN)      = {'trialN_zPerRun'};
                        modulators.(curr_modChosen_cond1)(:,idx_trialN)   = mod.trialN.zPerRun.(curr_modChosen_cond2).main;
                    case 9
                        modName.(curr_modChosen_cond1)(idx_trialN)      = {'trialN_zAcrossRuns'};
                        modulators.(curr_modChosen_cond1)(:,idx_trialN)   = mod.trialN.zAcrossRuns.(curr_modChosen_cond2).main;
                end
            end
            
            
        end % condition loop gain/neutral/loss pair
otherwise
    error('condition not ready yet');
end % onset chosen option circled in red

%% onset feedback
o_fbk = GLMprm.o_fbk;
dur_fbk = GLMprm.dur_fbk;
mod_fbk = GLMprm.mod_fbk;
switch o_fbk
    case 0
    case {1,2,3,4,5,6}
        
        switch o_fbk
            case 1 % all onsets grouped
                curr_modFbk_names = 'RL_fbk';
                %                 curr_modFbk_names2.RL_fbk = {'allPairs'};
                %                 curr_modFbk_names3.RL_fbk.allPairs = 'main';
                
                onsetCond.(curr_modFbk_names) = onset.allPairs.feedback.main;
                
                % duration
                switch dur_fbk
                    case 0 % stick
                        durCond.(curr_modFbk_names) = 0;
                    case 1 % boxcar
                        durCond.(curr_modFbk_names) = duration.allPairs.feedback.main;
                end
                
                % modulators
                Nmod.(curr_modFbk_names) = 0;
                
                %% trial number
                if mod_fbk.trialN ~= 0
                    Nmod.(curr_modFbk_names) = Nmod.(curr_modFbk_names) + 1;
                    idx_trialN = Nmod.(curr_modFbk_names);
                    switch mod_fbk.trialN
                        case 1
                            modName.(curr_modFbk_names)(idx_trialN)          = {'trialN_raw'};
                            modulators.(curr_modFbk_names)(:,idx_trialN)     = mod.trialN.raw.allPairs.main;
                        case 2
                            modName.(curr_modFbk_names)(idx_trialN)          = {'trialN_zPerRun'};
                            modulators.(curr_modFbk_names)(:,idx_trialN)     = mod.trialN.zPerRun.allPairs.main;
                        case 3
                            modName.(curr_modFbk_names)(idx_trialN)          = {'trialN_zAcrossRuns'};
                            modulators.(curr_modFbk_names)(:,idx_trialN)     = mod.trialN.zAcrossRuns.allPairs.main;
                        otherwise
                            error('not ready yet');
                    end
                end
                
                %% feedback
                if mod_fbk.fbk ~= 0
                    Nmod.(curr_modFbk_names) = Nmod.(curr_modFbk_names) + 1;
                    idx_fbk = Nmod.(curr_modFbk_names);
                    switch mod_fbk.fbk
                        case 1
                            modName.(curr_modFbk_names)(idx_fbk)      = {'feedback'};
                            modulators.(curr_modFbk_names)(:,idx_fbk)   = mod.fbk.raw.allPairs.main;
                        otherwise
                            error('not ready yet');
                    end
                end
                
                %% PE
                if mod_fbk.PE ~= 0
                    mdl_type_toUse = [mod_fbk.mdl_type,'_VBA_models'];
                    n_Qmodel_toUse = mod_fbk.mdl_n;
                    Nmod.(curr_modFbk_names) = Nmod.(curr_modFbk_names) + 1;
                    idx_PE = Nmod.(curr_modFbk_names);
                    switch mod_fbk.PE
                        case 1
                            modName.(curr_modFbk_names)(idx_PE)      = {'PE'};
                            modulators.(curr_modFbk_names)(:,idx_PE)   = mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.PE.allPairs.allPairsTrials;
                            case 2
                            modName.(curr_modFbk_names)(idx_PE)      = {'absPE'};
                            modulators.(curr_modFbk_names)(:,idx_PE)   = abs(mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.PE.allPairs.allPairsTrials);
                        otherwise
                            error('not ready yet');
                    end
                end
                %% PE bis
                if mod_fbk.PE_bis ~= 0
                    mdl_type_toUse = [mod_fbk.mdl_type,'_VBA_models'];
                    n_Qmodel_toUse = mod_fbk.mdl_n;
                    Nmod.(curr_modFbk_names) = Nmod.(curr_modFbk_names) + 1;
                    idx_PE_bis = Nmod.(curr_modFbk_names);
                    switch mod_fbk.PE
                        case 1
                            modName.(curr_modFbk_names)(idx_PE_bis)      = {'PE'};
                            modulators.(curr_modFbk_names)(:,idx_PE_bis)   = mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.PE.allPairs.allPairsTrials;
                            case 2
                            modName.(curr_modFbk_names)(idx_PE_bis)      = {'absPE'};
                            modulators.(curr_modFbk_names)(:,idx_PE_bis)   = abs(mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.PE.allPairs.allPairsTrials);
                        otherwise
                            error('not ready yet');
                    end
                end
                %% total gain
                if mod_fbk.totalGain ~= 0
                    Nmod.(curr_modFbk_names) = Nmod.(curr_modFbk_names) + 1;
                    idx_totalGain = Nmod.(curr_modFbk_names);
                    switch mod_fbk.totalGain
                        case 1
                            modName.(curr_modFbk_names)(idx_totalGain)      = {'totalGain'};
                            modulators.(curr_modFbk_names)(:,idx_totalGain)   = mod.totalGain.raw.allPairs.main;
                        otherwise
                            error('not ready yet');
                    end
                end
                
                %% ROI activity
                if mod_fbk.ROI_activity_yn ~= 0
                    Nmod.(curr_modFbk_names) = Nmod.(curr_modFbk_names) + 1;
                    idx_oFbk_ROI = Nmod.(curr_modFbk_names);
                    modName.(curr_modFbk_names)(idx_oFbk_ROI) = {[]};
                    % get the data
                    trialN_idx = mod.trialN.raw.(curr_modFbk_names).main; % trial index
                    ROI_nm          = mod_fbk.ROI_activity_ROI_nm; % ROI name
                    GLM_ROI_nm      = mod_fbk.ROI_activity_GLM;
                    ROI_period_nm   = mod_fbk.ROI_activity_period; % trial period when ROI extracted
                    jRun_nb = (iRun == 1)*1 + (iRun == 4)*2 + (iRun == 7)*3;
                    jRun_nm = ['run',num2str(jRun_nb)];
                    ROI_data = getfield(load([subj_analysis_folder,'ROI',filesep,'GLM',GLM_ROI_nm,'_',ROI_nm,'.mat'],...
                        'sub_ROI_trial_b_trial'),'sub_ROI_trial_b_trial');
                    modulators.(curr_modFbk_names)(:,idx_oFbk_ROI) = ROI_data.(ROI_nm).L.(ROI_period_nm).(jRun_nm).(sub_nm)(trialN_idx);
                end
                
            case 2 % onsets separated per pair: gain pair/neutral pair/loss pair
                curr_modFbk_names = {'RL_fbk_gainPair','RL_fbk_ntalPair', 'RL_fbk_lossPair'};
                curr_modFbk_names2.RL_fbk_gainPair  = 'gainPair';
                curr_modFbk_names2.RL_fbk_ntalPair  = 'ntalPair';
                curr_modFbk_names2.RL_fbk_lossPair  = 'lossPair';
                %                 curr_modFbk_names3.RL_fbk_gainPair.gainPair     = 'main';
                %                 curr_modFbk_names3.RL_fbk_ntalPair.ntalPair  = 'main';
                %                 curr_modFbk_names3.RL_fbk_lossPair.lossPair     = 'main';
                
                for iCond_modFbk = 1:length(curr_modFbk_names)
                    curr_modFbk_name1 = curr_modFbk_names{iCond_modFbk};
                    curr_modFbk_name2 = curr_modFbk_names2.(curr_modFbk_name1);
                    
                    onsetCond.(curr_modFbk_name1) = onset.(curr_modFbk_name2).feedback.main;
                    
                    % duration
                    switch dur_fbk
                        case 0 % stick
                            durCond.(curr_modFbk_name1) = 0;
                        case 1 % boxcar
                            durCond.(curr_modFbk_name1) = duration.(curr_modFbk_name2).feedback.main;
                    end
                    
                    % modulators
                    Nmod.(curr_modFbk_name1) = 0;
                    
                    %% trial number
                    if mod_fbk.trialN ~= 0
                        Nmod.(curr_modFbk_name1) = Nmod.(curr_modFbk_name1) + 1;
                        idx_trialN = Nmod.(curr_modFbk_name1);
                        switch mod_fbk.trialN
                            case 1
                                modName.(curr_modFbk_name1)(idx_trialN)          = {'trialN_raw'};
                                modulators.(curr_modFbk_name1)(:,idx_trialN)     = mod.trialN.raw.(curr_modFbk_name2).main;
                            case 2
                                modName.(curr_modFbk_name1)(idx_trialN)          = {'trialN_zPerRun'};
                                modulators.(curr_modFbk_name1)(:,idx_trialN)     = mod.trialN.zPerRun.(curr_modFbk_name2).main;
                            case 3
                                modName.(curr_modFbk_name1)(idx_trialN)          = {'trialN_zAcrossRuns'};
                                modulators.(curr_modFbk_name1)(:,idx_trialN)     = mod.trialN.zAcrossRuns.(curr_modFbk_name2).main;
                            otherwise
                                error('not ready yet');
                        end
                    end % trial number
                    
                    % feedback and total gain only make sense for gain or
                    % loss condition, not for neutral case
                    if ~strcmp(curr_modFbk_name1,'RL_fbk_ntalPair')
                        
                        %% feedback
                        if mod_fbk.fbk ~= 0
                            Nmod.(curr_modFbk_name1) = Nmod.(curr_modFbk_name1) + 1;
                            idx_fbk = Nmod.(curr_modFbk_name1);
                            switch mod_fbk.fbk
                                case 1
                                    modName.(curr_modFbk_name1)(idx_fbk)      = {'feedback'};
                                    modulators.(curr_modFbk_name1)(:,idx_fbk)   = mod.fbk.raw.(curr_modFbk_name2).main;
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        %% PE
                        if mod_fbk.PE ~= 0
                            mdl_type_toUse = [mod_fbk.mdl_type,'_VBA_models'];
                            n_Qmodel_toUse = mod_fbk.mdl_n;
                            Nmod.(curr_modFbk_name1) = Nmod.(curr_modFbk_name1) + 1;
                            idx_PE = Nmod.(curr_modFbk_name1);
                            switch mod_fbk.PE
                                case 1
                                    modName.(curr_modFbk_name1)(idx_PE)      = {'PE'};
                                    modulators.(curr_modFbk_name1)(:,idx_PE) = mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.PE.(curr_modFbk_name2).([curr_modFbk_name2,'Trials']);
                                case 2
                                    modName.(curr_modFbk_name1)(idx_PE)      = {'PE'};
                                    modulators.(curr_modFbk_name1)(:,idx_PE) = abs(mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.PE.(curr_modFbk_name2).([curr_modFbk_name2,'Trials']));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        %% PE bis
                        if mod_fbk.PE_bis ~= 0
                            mdl_type_toUse = [mod_fbk.mdl_type,'_VBA_models'];
                            n_Qmodel_toUse = mod_fbk.mdl_n;
                            Nmod.(curr_modFbk_name1) = Nmod.(curr_modFbk_name1) + 1;
                            idx_PE_bis = Nmod.(curr_modFbk_name1);
                            switch mod_fbk.PE_bis
                                case 1
                                    modName.(curr_modFbk_name1)(idx_PE_bis)      = {'PE'};
                                    modulators.(curr_modFbk_name1)(:,idx_PE_bis) = mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.PE.(curr_modFbk_name2).([curr_modFbk_name2,'Trials']);
                                case 2
                                    modName.(curr_modFbk_name1)(idx_PE_bis)      = {'PE'};
                                    modulators.(curr_modFbk_name1)(:,idx_PE_bis) = abs(mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.PE.(curr_modFbk_name2).([curr_modFbk_name2,'Trials']));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                    end % filter for neutral pairs
                    
                    %% total gain
                    if mod_fbk.totalGain ~= 0
                        Nmod.(curr_modFbk_name1) = Nmod.(curr_modFbk_name1) + 1;
                        idx_totalGain = Nmod.(curr_modFbk_name1);
                        switch mod_fbk.totalGain
                            case 1
                                modName.(curr_modFbk_name1)(idx_totalGain)      = {'totalGain'};
                                modulators.(curr_modFbk_name1)(:,idx_totalGain)   = mod.totalGain.raw.(curr_modFbk_name2).main;
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    %% ROI activity
                    if mod_fbk.ROI_activity_yn ~= 0
                        Nmod.(curr_modFbk_name1) = Nmod.(curr_modFbk_name1) + 1;
                        idx_oFbk_ROI = Nmod.(curr_modFbk_name1);
                        modName.(curr_modFbk_name1)(idx_oFbk_ROI) = {[]};
                        % get the data
                        trialN_idx = mod.trialN.raw.(curr_modFbk_name2).main; % trial index
                        ROI_nm          = mod_fbk.ROI_activity_ROI_nm; % ROI name
                        GLM_ROI_nm      = mod_fbk.ROI_activity_GLM;
                        ROI_period_nm   = mod_fbk.ROI_activity_period; % trial period when ROI extracted
                        jRun_nb = (iRun == 1)*1 + (iRun == 4)*2 + (iRun == 7)*3;
                        jRun_nm = ['run',num2str(jRun_nb)];
                        ROI_data = getfield(load([subj_analysis_folder,'ROI',filesep,'GLM',GLM_ROI_nm,'_',ROI_nm,'.mat'],...
                            'sub_ROI_trial_b_trial'),'sub_ROI_trial_b_trial');
                        modulators.(curr_modFbk_name1)(:,idx_oFbk_ROI) = ROI_data.(ROI_nm).L.(ROI_period_nm).(jRun_nm).(sub_nm)(trialN_idx);
                    end
                end % condition loop gain pair/neutral pair/loss pair
                
                
            case 3 % onsets separated per feedback: gain feedback/neutral feedback/loss feedback
                curr_modFbk_names = {'RL_fbk_gainFbk','RL_fbk_ntalFbk', 'RL_fbk_lossFbk'};
                %                 curr_modFbk_names2.RL_fbk_gainFbk   = {'gainPair'};
                %                 curr_modFbk_names2.RL_fbk_ntalFbk   = {'gainPair','ntalPair','lossPair'};
                %                 curr_modFbk_names2.RL_fbk_lossFbk   = {'lossPair'};
                %                 curr_modFbk_names3.RL_fbk_gainFbk.gainPair      = 'GainTrial';
                %                 curr_modFbk_names3.RL_fbk_ntalFbk.gainPair      = 'NeutralTrial';
                %                 curr_modFbk_names3.RL_fbk_ntalFbk.ntalPair   = 'NeutralTrial';
                %                 curr_modFbk_names3.RL_fbk_ntalFbk.lossPair      = 'NeutralTrial';
                %                 curr_modFbk_names3.RL_fbk_lossFbk.lossPair      = 'LossTrial';
                
                for iCond_modFbk = 1:length(curr_modFbk_names)
                    curr_modFbk_name1 = curr_modFbk_names{iCond_modFbk};
                    
                    switch curr_modFbk_name1
                        case {'RL_fbk_gainFbk','RL_fbk_lossFbk'}
                            if strcmp(curr_modFbk_name1,'RL_fbk_gainFbk')
                                fbk_pairType = 'gainPair';
                                fbk_trialType = 'GainTrial';
                            elseif strcmp(curr_modFbk_name1,'RL_fbk_lossFbk')
                                fbk_pairType = 'lossPair';
                                fbk_trialType = 'LossTrial';
                            end
                            onsetCond.(curr_modFbk_name1) = onset.(fbk_pairType).feedback.(fbk_trialType);
                            
                            % duration
                            switch dur_fbk
                                case 0 % stick
                                    durCond.(curr_modFbk_name1) = 0;
                                case 1 % boxcar
                                    durCond.(curr_modFbk_name1) = duration.(fbk_pairType).feedback.(fbk_trialType);
                            end
                            
                            % modulators
                            Nmod.(curr_modFbk_name1) = 0;
                            
                            %% trial number
                            if mod_fbk.trialN ~= 0
                                Nmod.(curr_modFbk_name1) = Nmod.(curr_modFbk_name1) + 1;
                                idx_trialN = Nmod.(curr_modFbk_name1);
                                switch mod_fbk.trialN
                                    case 1
                                        modName.(curr_modFbk_name1)(idx_trialN)          = {'trialN_raw'};
                                        modulators.(curr_modFbk_name1)(:,idx_trialN)     = mod.trialN.raw.(fbk_pairType).(fbk_trialType);
                                    case 2
                                        modName.(curr_modFbk_name1)(idx_trialN)          = {'trialN_zPerRun'};
                                        modulators.(curr_modFbk_name1)(:,idx_trialN)     = mod.trialN.zPerRun.(fbk_pairType).(fbk_trialType);
                                    case 3
                                        modName.(curr_modFbk_name1)(idx_trialN)          = {'trialN_zAcrossRuns'};
                                        modulators.(curr_modFbk_name1)(:,idx_trialN)     = mod.trialN.zAcrossRuns.(fbk_pairType).(fbk_trialType);
                                    otherwise
                                        error('not ready yet');
                                end
                            end % trial number
                            
                            %% feedback
                            if mod_fbk.fbk ~= 0
                                Nmod.(curr_modFbk_name1) = Nmod.(curr_modFbk_name1) + 1;
                                idx_fbk = Nmod.(curr_modFbk_name1);
                                switch mod_fbk.fbk
                                    case 1
                                        modName.(curr_modFbk_name1)(idx_fbk)      = {'feedback'};
                                        modulators.(curr_modFbk_name1)(:,idx_fbk)   = mod.fbk.raw.(fbk_pairType).(fbk_trialType);
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                            
                            %% PE
                            if mod_fbk.PE ~= 0
                                mdl_type_toUse = [mod_fbk.mdl_type,'_VBA_models'];
                                n_Qmodel_toUse = mod_fbk.mdl_n;
                                Nmod.(curr_modFbk_name1) = Nmod.(curr_modFbk_name1) + 1;
                                idx_PE = Nmod.(curr_modFbk_name1);
                                switch mod_fbk.PE
                                    case 1
                                        modName.(curr_modFbk_name1)(idx_PE)      = {'PE'};
                                        modulators.(curr_modFbk_name1)(:,idx_PE) = mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.PE.(fbk_pairType).(fbk_trialType);
                                    case 2
                                        modName.(curr_modFbk_name1)(idx_PE)      = {'absPE'};
                                        modulators.(curr_modFbk_name1)(:,idx_PE) = abs(mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.PE.(fbk_pairType).(fbk_trialType));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                            
                            %% PE bis
                            if mod_fbk.PE_bis ~= 0
                                mdl_type_toUse = [mod_fbk.mdl_type,'_VBA_models'];
                                n_Qmodel_toUse = mod_fbk.mdl_n;
                                Nmod.(curr_modFbk_name1) = Nmod.(curr_modFbk_name1) + 1;
                                idx_PE_bis = Nmod.(curr_modFbk_name1);
                                switch mod_fbk.PE
                                    case 1
                                        modName.(curr_modFbk_name1)(idx_PE_bis)      = {'PE'};
                                        modulators.(curr_modFbk_name1)(:,idx_PE_bis) = mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.PE.(fbk_pairType).(fbk_trialType);
                                    case 2
                                        modName.(curr_modFbk_name1)(idx_PE_bis)      = {'absPE'};
                                        modulators.(curr_modFbk_name1)(:,idx_PE_bis) = abs(mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.PE.(fbk_pairType).(fbk_trialType));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                            
                            %% total gain
                            if mod_fbk.totalGain ~= 0
                                Nmod.(curr_modFbk_name1) = Nmod.(curr_modFbk_name1) + 1;
                                idx_totalGain = Nmod.(curr_modFbk_name1);
                                switch mod_fbk.totalGain
                                    case 1
                                        modName.(curr_modFbk_name1)(idx_totalGain)      = {'totalGain'};
                                        modulators.(curr_modFbk_name1)(:,idx_totalGain)   = mod.totalGain.raw.(fbk_pairType).(fbk_trialType);
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                            
                            %% ROI activity
                            if mod_fbk.ROI_activity_yn ~= 0
                                Nmod.(curr_modFbk_name1) = Nmod.(curr_modFbk_name1) + 1;
                                idx_oFbk_ROI = Nmod.(curr_modFbk_name1);
                                modName.(curr_modFbk_name1)(idx_oFbk_ROI) = {[]};
                                % get the data
                                trialN_idx = mod.trialN.raw.(fbk_pairType).(fbk_trialType); % trial index
                                ROI_nm          = mod_fbk.ROI_activity_ROI_nm; % ROI name
                                GLM_ROI_nm      = mod_fbk.ROI_activity_GLM;
                                ROI_period_nm   = mod_fbk.ROI_activity_period; % trial period when ROI extracted
                                jRun_nb = (iRun == 1)*1 + (iRun == 4)*2 + (iRun == 7)*3;
                                jRun_nm = ['run',num2str(jRun_nb)];
                                ROI_data = getfield(load([subj_analysis_folder,'ROI',filesep,'GLM',GLM_ROI_nm,'_',ROI_nm,'.mat'],...
                                    'sub_ROI_trial_b_trial'),'sub_ROI_trial_b_trial');
                                modulators.(curr_modFbk_name1)(:,idx_oFbk_ROI) = ROI_data.(ROI_nm).L.(ROI_period_nm).(jRun_nm).(sub_nm)(trialN_idx);
                            end
                            
                        case 'RL_fbk_ntalFbk'
                            curr_onset = [onset.gainPair.feedback.NeutralTrial;...
                                onset.ntalPair.feedback.NeutralTrial;...
                                onset.lossPair.feedback.NeutralTrial];
                            curr_dur = [duration.gainPair.feedback.NeutralTrial;...
                                duration.ntalPair.feedback.NeutralTrial;...
                                duration.lossPair.feedback.NeutralTrial];
                            
                            % sort by ascending order and get index to fix
                            % modulators order accordingly
                            [onsetCond.(curr_modFbk_name1), ntalTrials_idx] = sort(curr_onset);
                            
                            % duration
                            switch dur_fbk
                                case 0 % stick
                                    durCond.(curr_modFbk_name1) = 0;
                                case 1 % boxcar
                                    durCond.(curr_modFbk_name1) = curr_dur(ntalTrials_idx);
                            end
                            
                            % modulators
                            Nmod.(curr_modFbk_name1) = 0;
                            
                            %% trial number
                            if mod_fbk.trialN ~= 0
                                Nmod.(curr_modFbk_name1) = Nmod.(curr_modFbk_name1) + 1;
                                idx_trialN = Nmod.(curr_modFbk_name1);
                                switch mod_fbk.trialN
                                    case 1
                                        modName.(curr_modFbk_name1)(idx_trialN)          = {'trialN_raw'};
                                        curr_trial_nb = [mod.trialN.raw.gainPair.NeutralTrial,...
                                            mod.trialN.raw.ntalPair.NeutralTrial,...
                                            mod.trialN.raw.lossPair.NeutralTrial];
                                    case 2
                                        modName.(curr_modFbk_name1)(idx_trialN)          = {'trialN_zPerRun'};
                                        curr_trial_nb = [mod.trialN.zPerRun.gainPair.NeutralTrial,...
                                            mod.trialN.zPerRun.ntalPair.NeutralTrial,...
                                            mod.trialN.zPerRun.lossPair.NeutralTrial];
                                    case 3
                                        modName.(curr_modFbk_name1)(idx_trialN)          = {'trialN_zAcrossRuns'};
                                        curr_trial_nb = [mod.trialN.zAcrossRuns.gainPair.NeutralTrial,...
                                            mod.trialN.zAcrossRuns.ntalPair.NeutralTrial,...
                                            mod.trialN.zAcrossRuns.lossPair.NeutralTrial];
                                    otherwise
                                        error('not ready yet');
                                end
                                modulators.(curr_modFbk_name1)(:,idx_trialN)     = curr_trial_nb(ntalTrials_idx);
                            end % trial number
                            
                            %% feedback = pointless, always zero in this case...
                            
                            %% PE (may be different from zero in gain and loss pair conditions)
                            if mod_fbk.PE ~= 0
                                mdl_type_toUse = [mod_fbk.mdl_type,'_VBA_models'];
                                n_Qmodel_toUse = mod_fbk.mdl_n;
                                Nmod.(curr_modFbk_name1) = Nmod.(curr_modFbk_name1) + 1;
                                idx_PE = Nmod.(curr_modFbk_name1);
                                switch mod_fbk.PE
                                    case 1
                                        modName.(curr_modFbk_name1)(idx_PE)      = {'PE'};
                                        curr_PE = [mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.PE.gainPair.NeutralTrial,...
                                            mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.PE.ntalPair.NeutralTrial,...
                                            mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.PE.lossPair.NeutralTrial];
                                        modulators.(curr_modFbk_name1)(:,idx_PE)     = curr_PE(ntalTrials_idx);
                                    case 2
                                        modName.(curr_modFbk_name1)(idx_PE)      = {'absPE'};
                                        curr_PE = [mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.PE.gainPair.NeutralTrial,...
                                            mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.PE.ntalPair.NeutralTrial,...
                                            mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.PE.lossPair.NeutralTrial];
                                        modulators.(curr_modFbk_name1)(:,idx_PE)     = abs(curr_PE(ntalTrials_idx));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                            
                            %% PE bis (may be different from zero in gain and loss pair conditions)
                            if mod_fbk.PE_bis ~= 0
                                mdl_type_toUse = [mod_fbk.mdl_type,'_VBA_models'];
                                n_Qmodel_toUse = mod_fbk.mdl_n;
                                Nmod.(curr_modFbk_name1) = Nmod.(curr_modFbk_name1) + 1;
                                idx_PE_bis = Nmod.(curr_modFbk_name1);
                                switch mod_fbk.PE_bis
                                    case 1
                                        modName.(curr_modFbk_name1)(idx_PE_bis)      = {'PE'};
                                        curr_PE = [mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.PE.gainPair.NeutralTrial,...
                                            mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.PE.ntalPair.NeutralTrial,...
                                            mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.PE.lossPair.NeutralTrial];
                                        modulators.(curr_modFbk_name1)(:,idx_PE_bis)     = curr_PE(ntalTrials_idx);
                                    case 2
                                        modName.(curr_modFbk_name1)(idx_PE_bis)      = {'absPE'};
                                        curr_PE = [mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.PE.gainPair.NeutralTrial,...
                                            mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.PE.ntalPair.NeutralTrial,...
                                            mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.PE.lossPair.NeutralTrial];
                                        modulators.(curr_modFbk_name1)(:,idx_PE_bis)     = abs(curr_PE(ntalTrials_idx));
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                            
                            %% total gain
                            if mod_fbk.totalGain ~= 0
                                Nmod.(curr_modFbk_name1) = Nmod.(curr_modFbk_name1) + 1;
                                idx_totalGain = Nmod.(curr_modFbk_name1);
                                switch mod_fbk.totalGain
                                    case 1
                                        modName.(curr_modFbk_name1)(idx_totalGain)      = {'totalGain'};
                                        curr_totalGain = [mod.totalGain.raw.gainPair.NeutralTrial,...
                                            mod.totalGain.raw.ntalPair.NeutralTrial,...
                                            mod.totalGain.raw.lossPair.NeutralTrial];
                                        modulators.(curr_modFbk_name1)(:,idx_totalGain)   = curr_totalGain;
                                    otherwise
                                        error('not ready yet');
                                end
                            end
                            
                        %% ROI activity
                        if mod_fbk.ROI_activity_yn ~= 0
                            Nmod.(curr_modFbk_name1) = Nmod.(curr_modFbk_name1) + 1;
                            idx_oFbk_ROI = Nmod.(curr_modFbk_name1);
                            modName.(curr_modFbk_name1)(idx_oFbk_ROI) = {[]};
                            % get the data
                            trialN_idx = [mod.trialN.raw.gainPair.NeutralTrial,...
                                mod.trialN.raw.ntalPair.NeutralTrial,...
                                mod.trialN.raw.lossPair.NeutralTrial]; % trial index
                            ROI_nm          = mod_fbk.ROI_activity_ROI_nm; % ROI name
                            GLM_ROI_nm      = mod_fbk.ROI_activity_GLM;
                            ROI_period_nm   = mod_fbk.ROI_activity_period; % trial period when ROI extracted
                            jRun_nb = (iRun == 1)*1 + (iRun == 4)*2 + (iRun == 7)*3;
                            jRun_nm = ['run',num2str(jRun_nb)];
                            ROI_data = getfield(load([subj_analysis_folder,'ROI',filesep,'GLM',GLM_ROI_nm,'_',ROI_nm,'.mat'],...
                                'sub_ROI_trial_b_trial'),'sub_ROI_trial_b_trial');
                            modulators.(curr_modFbk_name1)(:,idx_oFbk_ROI) = ROI_data.(ROI_nm).L.(ROI_period_nm).(jRun_nm).(sub_nm)(trialN_idx);
                        end
                    end % condition type
                    
                    
                end % condition loop gain feedback/neutral feedback/loss feedback
                
            case 4 % onsets separated per feedback per pair: gain feedback gain pair/ neutral feedback gain pair/ neutral feedback neutral pair/ neutral feedback loss pair/ loss feedback loss pair
                curr_modFbk_names = {'RL_fbk_gainPair_gainFbk','RL_fbk_gainPair_ntalFbk',...
                    'RL_fbk_ntalPair',...
                    'RL_fbk_lossPair_ntalFbk','RL_fbk_lossPair_lossFbk'};
                curr_modFbk_names2.RL_fbk_gainPair_gainFbk  = {'gainPair',      'GainTrial'};
                curr_modFbk_names2.RL_fbk_gainPair_ntalFbk  = {'gainPair',      'NeutralTrial'};
                curr_modFbk_names2.RL_fbk_ntalPair          = {'ntalPair',      'NeutralTrial'};
                curr_modFbk_names2.RL_fbk_lossPair_ntalFbk  = {'LossPair',      'NeutralTrial'};
                curr_modFbk_names2.RL_fbk_lossPair_lossFbk  = {'LossPair',      'LossTrial'};
                
                for iCond_modFbk = 1:length(curr_modFbk_names)
                    curr_modFbk_name1 = curr_modFbk_names{iCond_modFbk};
                    fbk_pairType = curr_modFbk_names2.(curr_modFbk_name1){1};
                    fbk_fbkType = curr_modFbk_names2.(curr_modFbk_name1){2};
                    
                    onsetCond.(curr_modFbk_name1) = onset.(fbk_pairType).feedback.(fbk_fbkType);
                    
                    % duration
                    switch dur_fbk
                        case 0 % stick
                            durCond.(curr_modFbk_name1) = 0;
                        case 1 % boxcar
                            durCond.(curr_modFbk_name1) = duration.(fbk_pairType).feedback.(fbk_fbkType);
                    end
                    
                    % modulators
                    Nmod.(curr_modFbk_name1) = 0;
                    
                    %% trial number
                    if mod_fbk.trialN ~= 0
                        Nmod.(curr_modFbk_name1) = Nmod.(curr_modFbk_name1) + 1;
                        idx_trialN = Nmod.(curr_modFbk_name1);
                        switch mod_fbk.trialN
                            case 1
                                modName.(curr_modFbk_name1)(idx_trialN)          = {'trialN_raw'};
                                modulators.(curr_modFbk_name1)(:,idx_trialN)     = mod.trialN.raw.(fbk_pairType).(fbk_fbkType);
                            case 2
                                modName.(curr_modFbk_name1)(idx_trialN)          = {'trialN_zPerRun'};
                                modulators.(curr_modFbk_name1)(:,idx_trialN)     = mod.trialN.zPerRun.(fbk_pairType).(fbk_fbkType);
                            case 3
                                modName.(curr_modFbk_name1)(idx_trialN)          = {'trialN_zAcrossRuns'};
                                modulators.(curr_modFbk_name1)(:,idx_trialN)     = mod.trialN.zAcrossRuns.(fbk_pairType).(fbk_fbkType);
                            otherwise
                                error('not ready yet');
                        end
                    end % trial number
                    
                    %% feedback
                    if mod_fbk.fbk ~= 0
                        error('Pointless: feedback = stable for a given feedback type.');
                    end
                    
                    %% PE
                    if ~strcmp(fbk_pairType,'ntalPair') % no change for neutral pair
                        if mod_fbk.PE ~= 0
                            mdl_type_toUse = [mod_fbk.mdl_type,'_VBA_models'];
                            n_Qmodel_toUse = mod_fbk.mdl_n;
                            Nmod.(curr_modFbk_name1) = Nmod.(curr_modFbk_name1) + 1;
                            idx_PE = Nmod.(curr_modFbk_name1);
                            switch mod_fbk.PE
                                case 1
                                    modName.(curr_modFbk_name1)(idx_PE)      = {'PE'};
                                    modulators.(curr_modFbk_name1)(:,idx_PE)   = mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.PE.(fbk_pairType).(fbk_fbkType);
                                case 2
                                    modName.(curr_modFbk_name1)(idx_PE)      = {'absPE'};
                                    modulators.(curr_modFbk_name1)(:,idx_PE)   = abs(mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.PE.(fbk_pairType).(fbk_fbkType));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                    end
                    
                    %% PE bis
                    if ~strcmp(fbk_pairType,'ntalPair') % no change for neutral pair
                        if mod_fbk.PE_bis ~= 0
                            mdl_type_toUse = [mod_fbk.mdl_type,'_VBA_models'];
                            n_Qmodel_toUse = mod_fbk.mdl_n;
                            Nmod.(curr_modFbk_name1) = Nmod.(curr_modFbk_name1) + 1;
                            idx_PE_bis = Nmod.(curr_modFbk_name1);
                            switch mod_fbk.PE_bis
                                case 1
                                    modName.(curr_modFbk_name1)(idx_PE_bis)      = {'PE'};
                                    modulators.(curr_modFbk_name1)(:,idx_PE_bis)   = mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.PE.(fbk_pairType).(fbk_fbkType);
                                case 2
                                    modName.(curr_modFbk_name1)(idx_PE_bis)      = {'absPE'};
                                    modulators.(curr_modFbk_name1)(:,idx_PE_bis)   = abs(mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.PE.(fbk_pairType).(fbk_fbkType));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                    end
                    
                    %% total gain
                    if mod_fbk.totalGain ~= 0
                        Nmod.(curr_modFbk_name1) = Nmod.(curr_modFbk_name1) + 1;
                        idx_totalGain = Nmod.(curr_modFbk_name1);
                        switch mod_fbk.totalGain
                            case 1
                                modName.(curr_modFbk_name1)(idx_totalGain)      = {'totalGain'};
                                modulators.(curr_modFbk_name1)(:,idx_totalGain)   = mod.totalGain.raw.(fbk_pairType).(fbk_fbkType);
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    %% ROI activity
                    if mod_fbk.ROI_activity_yn ~= 0
                        Nmod.(curr_modFbk_name1) = Nmod.(curr_modFbk_name1) + 1;
                        idx_oFbk_ROI = Nmod.(curr_modFbk_name1);
                        modName.(curr_modFbk_name1)(idx_oFbk_ROI) = {[]};
                        % get the data
                        trialN_idx = mod.trialN.raw.(fbk_pairType).(fbk_fbkType); % trial index
                        ROI_nm          = mod_fbk.ROI_activity_ROI_nm; % ROI name
                        GLM_ROI_nm      = mod_fbk.ROI_activity_GLM;
                        ROI_period_nm   = mod_fbk.ROI_activity_period; % trial period when ROI extracted
                        jRun_nb = (iRun == 1)*1 + (iRun == 4)*2 + (iRun == 7)*3;
                        jRun_nm = ['run',num2str(jRun_nb)];
                        ROI_data = getfield(load([subj_analysis_folder,'ROI',filesep,'GLM',GLM_ROI_nm,'_',ROI_nm,'.mat'],...
                            'sub_ROI_trial_b_trial'),'sub_ROI_trial_b_trial');
                        modulators.(curr_modFbk_name1)(:,idx_oFbk_ROI) = ROI_data.(ROI_nm).L.(ROI_period_nm).(jRun_nm).(sub_nm)(trialN_idx);
                    end
                    
                end % condition loop
                
            case 5 % pool gain and loss pair together but leave neutral pair apart
                curr_modFbk_names = {'RL_fbk_GL_Pairs', 'RL_fbk_ntalPair'};
                curr_modFbk_names2.RL_fbk_GL_Pairs  = 'GL_Pairs';
                curr_modFbk_names2.RL_fbk_ntalPair  = 'ntalPair';
                
                for iCond_modFbk = 1:length(curr_modFbk_names)
                    curr_modFbk_name1 = curr_modFbk_names{iCond_modFbk};
                    curr_modFbk_name2 = curr_modFbk_names2.(curr_modFbk_name1);
                    
                    onsetCond.(curr_modFbk_name1) = onset.(curr_modFbk_name2).feedback.main;
                    
                    % duration
                    switch dur_fbk
                        case 0 % stick
                            durCond.(curr_modFbk_name1) = 0;
                        case 1 % boxcar
                            durCond.(curr_modFbk_name1) = duration.(curr_modFbk_name2).feedback.main;
                    end
                    
                    % modulators
                    Nmod.(curr_modFbk_name1) = 0;
                    
                    %% trial number
                    if mod_fbk.trialN ~= 0
                        Nmod.(curr_modFbk_name1) = Nmod.(curr_modFbk_name1) + 1;
                        idx_trialN = Nmod.(curr_modFbk_name1);
                        switch mod_fbk.trialN
                            case 1
                                modName.(curr_modFbk_name1)(idx_trialN)          = {'trialN_raw'};
                                modulators.(curr_modFbk_name1)(:,idx_trialN)     = mod.trialN.raw.(curr_modFbk_name2).main;
                            case 2
                                modName.(curr_modFbk_name1)(idx_trialN)          = {'trialN_zPerRun'};
                                modulators.(curr_modFbk_name1)(:,idx_trialN)     = mod.trialN.zPerRun.(curr_modFbk_name2).main;
                            case 3
                                modName.(curr_modFbk_name1)(idx_trialN)          = {'trialN_zAcrossRuns'};
                                modulators.(curr_modFbk_name1)(:,idx_trialN)     = mod.trialN.zAcrossRuns.(curr_modFbk_name2).main;
                            otherwise
                                error('not ready yet');
                        end
                    end % trial number
                    
                    % feedback and total gain only make sense for gain or
                    % loss condition, not for neutral case
                    if ~strcmp(curr_modFbk_name1,'RL_fbk_ntalPair')
                        
                        %% feedback
                        if mod_fbk.fbk ~= 0
                            Nmod.(curr_modFbk_name1) = Nmod.(curr_modFbk_name1) + 1;
                            idx_fbk = Nmod.(curr_modFbk_name1);
                            switch mod_fbk.fbk
                                case 1
                                    modName.(curr_modFbk_name1)(idx_fbk)      = {'feedback'};
                                    modulators.(curr_modFbk_name1)(:,idx_fbk)   = mod.fbk.raw.(curr_modFbk_name2).main;
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        %% PE
                        if mod_fbk.PE ~= 0
                            mdl_type_toUse = [mod_fbk.mdl_type,'_VBA_models'];
                            n_Qmodel_toUse = mod_fbk.mdl_n;
                            Nmod.(curr_modFbk_name1) = Nmod.(curr_modFbk_name1) + 1;
                            idx_PE = Nmod.(curr_modFbk_name1);
                            switch mod_fbk.PE
                                case 1
                                    modName.(curr_modFbk_name1)(idx_PE)      = {'PE'};
                                    modulators.(curr_modFbk_name1)(:,idx_PE) = mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.PE.(curr_modFbk_name2).([curr_modFbk_name2,'Trials']);
                                case 2
                                    modName.(curr_modFbk_name1)(idx_PE)      = {'PE'};
                                    modulators.(curr_modFbk_name1)(:,idx_PE) = abs(mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.PE.(curr_modFbk_name2).([curr_modFbk_name2,'Trials']));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        %% PE bis
                        if mod_fbk.PE_bis ~= 0
                            mdl_type_toUse = [mod_fbk.mdl_type,'_VBA_models'];
                            n_Qmodel_toUse = mod_fbk.mdl_n;
                            Nmod.(curr_modFbk_name1) = Nmod.(curr_modFbk_name1) + 1;
                            idx_PE_bis = Nmod.(curr_modFbk_name1);
                            switch mod_fbk.PE_bis
                                case 1
                                    modName.(curr_modFbk_name1)(idx_PE_bis)      = {'PE'};
                                    modulators.(curr_modFbk_name1)(:,idx_PE_bis) = mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.PE.(curr_modFbk_name2).([curr_modFbk_name2,'Trials']);
                                case 2
                                    modName.(curr_modFbk_name1)(idx_PE_bis)      = {'PE'};
                                    modulators.(curr_modFbk_name1)(:,idx_PE_bis) = abs(mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.PE.(curr_modFbk_name2).([curr_modFbk_name2,'Trials']));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                    end % filter for neutral pairs
                    
                    %% total gain
                    if mod_fbk.totalGain ~= 0
                        Nmod.(curr_modFbk_name1) = Nmod.(curr_modFbk_name1) + 1;
                        idx_totalGain = Nmod.(curr_modFbk_name1);
                        switch mod_fbk.totalGain
                            case 1
                                modName.(curr_modFbk_name1)(idx_totalGain)      = {'totalGain'};
                                modulators.(curr_modFbk_name1)(:,idx_totalGain)   = mod.totalGain.raw.(curr_modFbk_name2).main;
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    %% ROI activity
                    if mod_fbk.ROI_activity_yn ~= 0
                        Nmod.(curr_modFbk_name1) = Nmod.(curr_modFbk_name1) + 1;
                        idx_oFbk_ROI = Nmod.(curr_modFbk_name1);
                        modName.(curr_modFbk_name1)(idx_oFbk_ROI) = {[]};
                        % get the data
                        trialN_idx = mod.trialN.raw.(curr_modFbk_name2).main; % trial index
                        ROI_nm          = mod_fbk.ROI_activity_ROI_nm; % ROI name
                        GLM_ROI_nm      = mod_fbk.ROI_activity_GLM;
                        ROI_period_nm   = mod_fbk.ROI_activity_period; % trial period when ROI extracted
                        jRun_nb = (iRun == 1)*1 + (iRun == 4)*2 + (iRun == 7)*3;
                        jRun_nm = ['run',num2str(jRun_nb)];
                        ROI_data = getfield(load([subj_analysis_folder,'ROI',filesep,'GLM',GLM_ROI_nm,'_',ROI_nm,'.mat'],...
                            'sub_ROI_trial_b_trial'),'sub_ROI_trial_b_trial');
                        modulators.(curr_modFbk_name1)(:,idx_oFbk_ROI) = ROI_data.(ROI_nm).L.(ROI_period_nm).(jRun_nm).(sub_nm)(trialN_idx);
                    end
                end % condition loop gain+loss pair/neutral pair
                
            case 6 % pool gain and loss pair together but leave neutral pair apart and also split first/second half trials
                curr_modFbk_names = {'RL_fbk_GL_Pairs_first',...
                    'RL_fbk_GL_Pairs_last',...
                    'RL_fbk_ntalPair'};
                curr_modFbk_names2.RL_fbk_GL_Pairs_first  = 'GLPairs_first';
                curr_modFbk_names2.RL_fbk_GL_Pairs_last  = 'GLPairs_last';
                curr_modFbk_names2.RL_fbk_ntalPair  = 'ntalPair';
                
                for iCond_modFbk = 1:length(curr_modFbk_names)
                    curr_modFbk_name1 = curr_modFbk_names{iCond_modFbk};
                    curr_modFbk_name2 = curr_modFbk_names2.(curr_modFbk_name1);
                    
                    onsetCond.(curr_modFbk_name1) = onset.(curr_modFbk_name2).feedback.main;
                    
                    % duration
                    switch dur_fbk
                        case 0 % stick
                            durCond.(curr_modFbk_name1) = 0;
                        case 1 % boxcar
                            durCond.(curr_modFbk_name1) = duration.(curr_modFbk_name2).feedback.main;
                    end
                    
                    % modulators
                    Nmod.(curr_modFbk_name1) = 0;
                    
                    %% trial number
                    if mod_fbk.trialN ~= 0
                        Nmod.(curr_modFbk_name1) = Nmod.(curr_modFbk_name1) + 1;
                        idx_trialN = Nmod.(curr_modFbk_name1);
                        switch mod_fbk.trialN
                            case 1
                                modName.(curr_modFbk_name1)(idx_trialN)          = {'trialN_raw'};
                                modulators.(curr_modFbk_name1)(:,idx_trialN)     = mod.trialN.raw.(curr_modFbk_name2).main;
                            case 2
                                modName.(curr_modFbk_name1)(idx_trialN)          = {'trialN_zPerRun'};
                                modulators.(curr_modFbk_name1)(:,idx_trialN)     = mod.trialN.zPerRun.(curr_modFbk_name2).main;
                            case 3
                                modName.(curr_modFbk_name1)(idx_trialN)          = {'trialN_zAcrossRuns'};
                                modulators.(curr_modFbk_name1)(:,idx_trialN)     = mod.trialN.zAcrossRuns.(curr_modFbk_name2).main;
                            otherwise
                                error('not ready yet');
                        end
                    end % trial number
                    
                    % feedback and total gain only make sense for gain or
                    % loss condition, not for neutral case
                    if ~strcmp(curr_modFbk_name1,'RL_fbk_ntalPair')
                        
                        %% feedback
                        if mod_fbk.fbk ~= 0
                            Nmod.(curr_modFbk_name1) = Nmod.(curr_modFbk_name1) + 1;
                            idx_fbk = Nmod.(curr_modFbk_name1);
                            switch mod_fbk.fbk
                                case 1
                                    modName.(curr_modFbk_name1)(idx_fbk)      = {'feedback'};
                                    modulators.(curr_modFbk_name1)(:,idx_fbk)   = mod.fbk.raw.(curr_modFbk_name2).main;
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        %% PE
                        if mod_fbk.PE ~= 0
                            mdl_type_toUse = [mod_fbk.mdl_type,'_VBA_models'];
                            n_Qmodel_toUse = mod_fbk.mdl_n;
                            Nmod.(curr_modFbk_name1) = Nmod.(curr_modFbk_name1) + 1;
                            idx_PE = Nmod.(curr_modFbk_name1);
                            switch mod_fbk.PE
                                case 1
                                    modName.(curr_modFbk_name1)(idx_PE)      = {'PE'};
                                    modulators.(curr_modFbk_name1)(:,idx_PE) = mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.PE.(curr_modFbk_name2).([curr_modFbk_name2,'Trials']);
                                case 2
                                    modName.(curr_modFbk_name1)(idx_PE)      = {'PE'};
                                    modulators.(curr_modFbk_name1)(:,idx_PE) = abs(mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.PE.(curr_modFbk_name2).([curr_modFbk_name2,'Trials']));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                        %% PE bis
                        if mod_fbk.PE_bis ~= 0
                            mdl_type_toUse = [mod_fbk.mdl_type,'_VBA_models'];
                            n_Qmodel_toUse = mod_fbk.mdl_n;
                            Nmod.(curr_modFbk_name1) = Nmod.(curr_modFbk_name1) + 1;
                            idx_PE_bis = Nmod.(curr_modFbk_name1);
                            switch mod_fbk.PE_bis
                                case 1
                                    modName.(curr_modFbk_name1)(idx_PE_bis)      = {'PE'};
                                    modulators.(curr_modFbk_name1)(:,idx_PE_bis) = mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.PE.(curr_modFbk_name2).([curr_modFbk_name2,'Trials']);
                                case 2
                                    modName.(curr_modFbk_name1)(idx_PE_bis)      = {'PE'};
                                    modulators.(curr_modFbk_name1)(:,idx_PE_bis) = abs(mod.(mdl_type_toUse).Q_model(n_Qmodel_toUse).raw.PE.(curr_modFbk_name2).([curr_modFbk_name2,'Trials']));
                                otherwise
                                    error('not ready yet');
                            end
                        end
                        
                    end % filter for neutral pairs
                    
                    %% total gain
                    if mod_fbk.totalGain ~= 0
                        Nmod.(curr_modFbk_name1) = Nmod.(curr_modFbk_name1) + 1;
                        idx_totalGain = Nmod.(curr_modFbk_name1);
                        switch mod_fbk.totalGain
                            case 1
                                modName.(curr_modFbk_name1)(idx_totalGain)      = {'totalGain'};
                                modulators.(curr_modFbk_name1)(:,idx_totalGain)   = mod.totalGain.raw.(curr_modFbk_name2).main;
                            otherwise
                                error('not ready yet');
                        end
                    end
                    
                    %% ROI activity
                    if mod_fbk.ROI_activity_yn ~= 0
                        Nmod.(curr_modFbk_name1) = Nmod.(curr_modFbk_name1) + 1;
                        idx_oFbk_ROI = Nmod.(curr_modFbk_name1);
                        modName.(curr_modFbk_name1)(idx_oFbk_ROI) = {[]};
                        % get the data
                        trialN_idx = mod.trialN.raw.(curr_modFbk_name2).main; % trial index
                        ROI_nm          = mod_fbk.ROI_activity_ROI_nm; % ROI name
                        GLM_ROI_nm      = mod_fbk.ROI_activity_GLM;
                        ROI_period_nm   = mod_fbk.ROI_activity_period; % trial period when ROI extracted
                        jRun_nb = (iRun == 1)*1 + (iRun == 4)*2 + (iRun == 7)*3;
                        jRun_nm = ['run',num2str(jRun_nb)];
                        ROI_data = getfield(load([subj_analysis_folder,'ROI',filesep,'GLM',GLM_ROI_nm,'_',ROI_nm,'.mat'],...
                            'sub_ROI_trial_b_trial'),'sub_ROI_trial_b_trial');
                        modulators.(curr_modFbk_name1)(:,idx_oFbk_ROI) = ROI_data.(ROI_nm).L.(ROI_period_nm).(jRun_nm).(sub_nm)(trialN_idx);
                    end
                end % condition loop gain+loss pair/neutral pair
                
        end % o_fbk
        
        %         n_cond = n_cond + length(curr_modFbk_names);
        cond_nm = [cond_nm, curr_modFbk_names];
        
otherwise
    error('condition not ready yet');
end % onset feedback

%% onset cross
o_cross = GLMprm.o_cross;
switch o_cross
case 0
    case 1 % onset for fixation cross across all trials and conditions
        %         n_cond = n_cond + 1;
        cond_nm = [cond_nm, 'RL_cross'];
        onsetCond.RL_cross = onset.ITI;
        
        % duration
        switch dur_fbk
            case 0 % stick
                durCond.RL_cross = 0;
            case 1 % boxcar
                durCond.RL_cross = duration.ITI;
        end
                    
        Nmod.RL_cross = 0;
otherwise
    error('condition not ready yet');
end

%% missed trials
o_missed_trials_stim = GLMprm.o_missed_trials_stim;
switch o_missed_trials_stim
    case 0
    case 1
        %         n_cond = n_cond + 1;
        cond_nm = [cond_nm, 'RL_missed'];
        onsetCond.RL_missed = onset.missedTrials.allPairs.dispOptions.main;
        
        % duration
        durCond.RL_missed = duration.missedTrials.allPairs.dispOptions.main;
        
        % simplest version: no modulator
        Nmod.RL_missed = 0;
    otherwise
        error('condition not ready yet');
end


%% preparing matlab batch
[ matlabbatch ] = LGC_First_level_RL_matlabbatch_onsets_modulators( matlabbatch,...
    cond_nm, onsetCond, durCond,...
    sub_idx, iRun, Nmod,...
    task_id, modName, modulators,...
    gal);

end % function end