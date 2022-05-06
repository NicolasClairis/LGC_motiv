function[con_names, con_vector] = LGCM_contrasts(study_nm, sub_nm, GLM,...
    computer_root, preproc_sm_kernel, condition)
% [con_names, con_vector] = LGCM_contrasts(study_nm, sub_nm, GLM,...
%   computer_root, preproc_sm_kernel, condition)
% LGCM_contrasts will define the contrast names and contrast vector for the
% subject and study entered in input.
%
% INPUTS
% computer_root: path where data is stored
%
% study_nm: definition of the study on which you want to analyze the data
% 'fMRI_pilots': pilots
% 'study1': first study (dmPFC + AI)
% 'study2': second study (clinical trial)
%
% sub_nm: subject name
%
% GLM: GLM number
%
% preproc_sm_kernel: kernel used in preprocessing for smoothing the data
%
% condition: define subjects and runs to include
% 'fMRI': all subjects where fMRI ok
% 'fMRI_no_move': remove runs with too much movement
%
% OUTPUTS
% con_names: list of contrast names
%
% con_vector: matrix with all contrasts organized in nb_contrasts (rows)*nb_regressors (columns)
%
% See also LGCM_contrasts_spm.m


%% working directories
switch study_nm
    case 'fMRI_pilots'
        root = fullfile(computer_root,'fMRI_pilots');
    case 'study1'
        root = fullfile(computer_root,'study1');
    case 'study2'
        root = fullfile(computer_root,'study2');
end
subj_folder             = [root, filesep, 'CID',sub_nm];
subj_analysis_folder    = [subj_folder, filesep, 'fMRI_analysis' filesep];
switch condition
    case {'fMRI','fMRI_no_move_bis'}
        resultsFolderName = [subj_analysis_folder 'functional', filesep,...
            'preproc_sm_',num2str(preproc_sm_kernel),'mm', filesep,...
            'GLM',num2str(GLM),filesep];
    case 'fMRI_no_move'
        resultsFolderName = [subj_analysis_folder 'functional', filesep,...
            'preproc_sm_',num2str(preproc_sm_kernel),'mm', filesep,...
            'GLM',num2str(GLM),'_no_movementRun',filesep];
end

%% extract GLM informations
[reg_names, n_regsPerTask] = GLM_details(GLM);
GLMprm = which_GLM(GLM);

%% define runs based on current subject
runs = runs_definition(study_nm, sub_nm, condition);

%% initialize variables of interest
con_names = {};
con_vector = [];
% define total number of regressors
% get number of regressors/task*number of runs/task + 1 constant/run
n_totalRegs = n_regsPerTask.Ep*(runs.nb_runs.Ep) + n_regsPerTask.Em*(runs.nb_runs.Em) + runs.nb_runs.Ep + runs.nb_runs.Em;

% then compare that n_regs is equal to the total number of betas output
betaFiles = ls([resultsFolderName,'beta_*.nii']);
n_1stLevelBetas = size(betaFiles,1);
if n_1stLevelBetas ~= n_totalRegs
    error(['Problem: number of regressors predicted by contrasts = ',num2str(n_totalRegs),...
        ' while number of betas produced by 1st level = ',num2str(n_1stLevelBetas)]);
else
    disp(['CID',sub_nm,' number of regressors and number of betas produced in 1st level - ok'])
end

%% extract run informations (to know which run corresponds to which task)
Ep_runs = strcmp(runs.tasks,'Ep');
Em_runs = strcmp(runs.tasks,'Em');

%% extract the contrasts of interest
jReg = 0;
%% physical effort
if runs.nb_runs.Ep > 0
    for iReg_Ep = 1:n_regsPerTask.Ep
        reg_nm = reg_names.Ep{iReg_Ep};
        
        if ~strcmp(reg_nm,'movement') && ~isempty(reg_nm) % ignore movement regressors (labelled as 'movement') and temporal derivative regressors (empty name)
            
            % run 1
            run1_Ep_vec = (strcmp(reg_names.Ep, reg_nm).*(Ep_runs(1) == 1)) +...
                (zeros(1,n_regsPerTask.Em).*(Em_runs(1) == 1));
            
            % run 2
            if runs.nb_runs.Ep + runs.nb_runs.Em >= 2
                run2_Ep_vec = (strcmp(reg_names.Ep, reg_nm).*(Ep_runs(2) == 1)) +...
                    (zeros(1,n_regsPerTask.Em).*(Em_runs(2) == 1));
            else
                run2_Ep_vec = [];
            end
            
            % run 3
            if runs.nb_runs.Ep + runs.nb_runs.Em >= 3
                run3_Ep_vec =  (strcmp(reg_names.Ep, reg_nm).*(Ep_runs(3) == 1)) +...
                    (zeros(1,n_regsPerTask.Em).*(Em_runs(3) == 1));
            else
                run3_Ep_vec = [];
            end
            
            % run 4
            if runs.nb_runs.Ep + runs.nb_runs.Em == 4
                run4_Ep_vec =  (strcmp(reg_names.Ep, reg_nm).*(Ep_runs(4) == 1)) +...
                    (zeros(1,n_regsPerTask.Em).*(Em_runs(4) == 1));
            else
                run4_Ep_vec = [];
            end
            
            con_vec_Ep_tmp = [run1_Ep_vec,...
                        run2_Ep_vec,...
                        run3_Ep_vec,...
                        run4_Ep_vec,...
                        zeros(1, runs.nb_runs.Ep), zeros(1, runs.nb_runs.Em)];
            
            % positive contrast
            jReg = jReg + 1;
            con_names{jReg} = ['Ep ',reg_nm];
            con_vector(jReg, 1:n_totalRegs) = con_vec_Ep_tmp;
            
            % negative contrast
            jReg = jReg + 1;
            con_names{jReg} = ['Ep -',reg_nm];
            con_vector(jReg, 1:n_totalRegs) = -con_vec_Ep_tmp;
            
        end % movement and derivative filter
    end % loop through physical effort regressors
end % physical effort performed at least in one full run?

%% mental effort
if runs.nb_runs.Em > 0
    for iReg_Em = 1:n_regsPerTask.Em
        reg_nm = reg_names.Em{iReg_Em};
        
        if ~strcmp(reg_nm,'movement') && ~isempty(reg_nm) % ignore movement regressors (labelled as 'movement') and temporal derivative regressors (empty name)
            
            % run 1
            run1_Em_vec = (strcmp(reg_names.Em, reg_nm).*(Em_runs(1) == 1)) +...
                (zeros(1,n_regsPerTask.Ep).*(Ep_runs(1) == 1));
            
            % run 2
            if runs.nb_runs.Ep + runs.nb_runs.Em >= 2
                run2_Em_vec = (strcmp(reg_names.Em, reg_nm).*(Em_runs(2) == 1)) +...
                    (zeros(1,n_regsPerTask.Ep).*(Ep_runs(2) == 1));
            else
                run2_Em_vec = [];
            end
            
            % run 3
            if runs.nb_runs.Ep + runs.nb_runs.Em >= 3
                run3_Em_vec =  (strcmp(reg_names.Em, reg_nm).*(Em_runs(3) == 1)) +...
                    (zeros(1,n_regsPerTask.Ep).*(Ep_runs(3) == 1));
            else
                run3_Em_vec = [];
            end
            
            % run 4
            if runs.nb_runs.Ep + runs.nb_runs.Em == 4
                run4_Em_vec =  (strcmp(reg_names.Em, reg_nm).*(Em_runs(4) == 1)) +...
                    (zeros(1,n_regsPerTask.Ep).*(Ep_runs(4) == 1));
            else
                run4_Em_vec = [];
            end
            
            con_vec_Em_tmp = [run1_Em_vec,...
                        run2_Em_vec,...
                        run3_Em_vec,...
                        run4_Em_vec,...
                        zeros(1, runs.nb_runs.Ep), zeros(1, runs.nb_runs.Em)];
            
            % positive contrast
            jReg = jReg + 1;
            con_names{jReg} = ['Em ',reg_nm];
            con_vector(jReg, 1:n_totalRegs) = con_vec_Em_tmp;
            
            % negative contrast
            jReg = jReg + 1;
            con_names{jReg} = ['Em -',reg_nm];
            con_vector(jReg, 1:n_totalRegs) = -con_vec_Em_tmp;
            
        end % movement and derivative filter
    end % loop through mental effort regressors
end % mental effort performed at least in one full run?

%% pool contrasts across both tasks
if runs.nb_runs.Em > 0 && runs.nb_runs.Ep > 0
    
    for iReg_Ep_bis = 1:n_regsPerTask.Ep
        reg_nm = reg_names.Ep{iReg_Ep_bis};
        if (sum(strcmp(reg_names.Em, reg_nm)) > 0) &&...
                ~strcmp(reg_nm,'movement') && ~isempty(reg_nm) % ignore movement regressors (labelled as 'movement') and temporal derivative regressors (empty name)
            % extract previously defined contrast for each task
            jReg_Ep = strcmp(con_names,['Ep ',reg_nm]);
            jReg_Em = strcmp(con_names,['Em ',reg_nm]);
            % pool the two contrasts
            con_vec_EpEm_tmp = con_vector(jReg_Ep, :) + con_vector(jReg_Em, :);
            
            % positive contrast
            jReg = jReg + 1;
            con_names{jReg} = ['Ep+Em ',reg_nm];
            con_vector(jReg, 1:n_totalRegs) = con_vec_EpEm_tmp;
            
            % negative contrast
            jReg = jReg + 1;
            con_names{jReg} = ['Ep+Em -',reg_nm];
            con_vector(jReg, 1:n_totalRegs) = -con_vec_EpEm_tmp;
        end % check if contrast exists for both tasks (+ ignore movement and derivative regressors)
    end % loop through Ep regressors
    
end % at least one run of each task has been performed?

%% pool contrasts across R and P trials if have been modelled separately
timePeriod = {'choice','chosen','Eperf','fbk'};
n_timePeriods = length(timePeriod);
timePeriod_names = {'choice','chosen','effort','feedback'};
regTypes = {'ONSET','REG'};
nRegTypes = length(regTypes);
tasks = {'Ep','Em','Ep+Em'};
nTasks = length(tasks);
% loop through time periods
for iTimePeriod = 1:n_timePeriods
    timePhase_nm = timePeriod{iTimePeriod};
    timePhase_nm_bis = timePeriod_names{iTimePeriod};
    
    % loop through tasks (including tasks pooled)
    for iTask = 1:nTasks
        task_nm = tasks{iTask};
        
        % check if there is a pool during this phase, otherwise, it is
        % worthless
        if (ismember(task_nm,{'Ep','Em'}) && GLMprm.(timePhase_nm).(task_nm).RPpool == 0) ||...
                (GLMprm.(timePhase_nm).Ep.RPpool == 0 && GLMprm.(timePhase_nm).Em.RPpool == 0)
            
            % pool through onsets and regressors
            for iRegType = 1:nRegTypes
                regType = regTypes{iRegType};
                
                % loop through positive and negative contrasts
                for iPosNeg = 1:2
                    switch iPosNeg
                        case 1
                            posNeg = '';
                        case 2
                            posNeg = '-';
                    end
                    % to avoid pooling R+P if already done, add ':' for REG
                    % specifically (no need for ONSET as full expression is
                    % already ok)
                    switch regType
                        case 'ONSET'
                            regRExpression = [task_nm,' ',posNeg,'ONSET ',timePhase_nm_bis,' R'];
                            regPExpression = [task_nm,' ',posNeg,'ONSET ',timePhase_nm_bis,' P'];
                            areThereRregs = strcmp(con_names, regRExpression);
                            areTherePregs = strcmp(con_names, regPExpression);
                        case 'REG'
                            regRExpression = [task_nm,' ',posNeg,regType,' ',timePhase_nm_bis,' R:'];
                            regPExpression = [task_nm,' ',posNeg,regType,' ',timePhase_nm_bis,' P:'];
                            n_regNameSize = length(regRExpression);
                            areThereRregs = strncmp(con_names, regRExpression, n_regNameSize);
                            areTherePregs = strncmp(con_names, regPExpression, n_regNameSize);
                    end
                    
                    n_RPregs = sum(areThereRregs);
                    if n_RPregs > 0
                        jRregList = find(areThereRregs);
                        jPregList = find(areTherePregs);
                        
                        % loop through regressors to pool
                        for iRPreg = 1:n_RPregs
                            Rreg_full_nm = con_names{jRregList(iRPreg)};
                            Preg_full_nm = con_names{jPregList(iRPreg)};
                            if strcmp(regType,'REG')
                                conR_nm = Rreg_full_nm(n_regNameSize+1:end);
                                conP_nm = Preg_full_nm(n_regNameSize+1:end);
                                if ~strcmp(conR_nm, conP_nm)
                                    error('problem mismatch R and P regressor for pooling. Please fix it.');
                                end % comparison R/P regressor
                            end
                            
                            % pool the two contrasts
                            con_vec_RP_tmp = con_vector(jRregList(iRPreg), :) +...
                                con_vector(jPregList(iRPreg), :);
                            
                            % positive contrast (no need for negative contrast
                            % for once, because will automatically be created
                            % by pooling negative R and P together)
                            jReg = jReg + 1;
                            switch regType
                                case 'ONSET'
                                    con_names{jReg} = [regRExpression,'P'];
                                case 'REG'
                                    con_names{jReg} = [regRExpression(1:(end-1)),'P:',conR_nm];
                            end
                            con_vector(jReg, 1:n_totalRegs) = con_vec_RP_tmp;
                            
                        end % loop through regressors to pool
                    end % filter: more than 0 regressors then ok
                end % positive/negative contrast
            end % ONSET/REG
        end % RPpool
    end % task loop
end % time loop

%% compare Vchosen - Vunchosen if Vch and Vunch have been modelled
Epm = {'Ep','Em'};
RP_cond = {'R','P','RP'};
task_phases = {'choice','chosen','effort'};
for iTaskPhase = 1:length(task_phases)
    taskPhase_nm = task_phases{iTaskPhase};
    for iRP = 1:length(RP_cond)
        RP_nm = RP_cond{iRP};
        
        % money chosen - money unchosen
        money_chosen_ConNm = ['REG ',taskPhase_nm,' ',RP_nm,': money chosen'];
        money_unchosen_ConNm = ['REG ',taskPhase_nm,' ',RP_nm,': money unchosen'];
        % effort chosen - effort unchosen
        effort_chosen_ConNm = ['REG ',taskPhase_nm,' ',RP_nm,': effort chosen'];
        effort_unchosen_ConNm = ['REG ',taskPhase_nm,' ',RP_nm,': effort unchosen'];
        
        for iTask = 1:length(Epm)
            task_nm = Epm{iTask};
            
            % money chosen - money unchosen
            if (sum(strcmp(reg_names.(task_nm), money_chosen_ConNm)) > 0 &&...
                    sum(strcmp(reg_names.(task_nm), money_unchosen_ConNm)) > 0)
                
                % extract contrast
                moneyCh_idx = strcmp(con_names,[task_nm,' ',money_chosen_ConNm]);
                moneyUnch_idx = strcmp(con_names,[task_nm,' ',money_unchosen_ConNm]);
                con_moneyCh_min_moneyUnch = con_vector(moneyCh_idx, :) - con_vector(moneyUnch_idx, :);
                % positive contrast
                jReg = jReg + 1;
                con_names{jReg} = [task_nm,' REG ',taskPhase_nm,' ',RP_nm,': money chosen - unchosen'];
                con_vector(jReg, 1:n_totalRegs) = con_moneyCh_min_moneyUnch;
                
                % negative contrast
                jReg = jReg + 1;
                con_names{jReg} = [task_nm,' REG ',taskPhase_nm,' ',RP_nm,': money unchosen - chosen'];
                con_vector(jReg, 1:n_totalRegs) = -con_moneyCh_min_moneyUnch;
            end % money chosen - money unchosen
            
            % effort chosen - effort unchosen
            if (sum(strcmp(reg_names.(task_nm), effort_chosen_ConNm)) > 0 &&...
                    sum(strcmp(reg_names.(task_nm), effort_unchosen_ConNm)) > 0)
                % extract contrast
                effortCh_idx = strcmp(con_names,[task_nm,' ',effort_chosen_ConNm]);
                effortUnch_idx = strcmp(con_names,[task_nm,' ',effort_unchosen_ConNm]);
                con_effortCh_min_effortUnch = con_vector(effortCh_idx, :) - con_vector(effortUnch_idx, :);
                % positive contrast
                jReg = jReg + 1;
                con_names{jReg} = [task_nm,' REG ',taskPhase_nm,' ',RP_nm,': effort chosen - unchosen'];
                con_vector(jReg, 1:n_totalRegs) = con_effortCh_min_effortUnch;
                
                % negative contrast
                jReg = jReg + 1;
                con_names{jReg} = [task_nm,' REG ',taskPhase_nm,' ',RP_nm,': effort unchosen - chosen'];
                con_vector(jReg, 1:n_totalRegs) = -con_effortCh_min_effortUnch;
            end % effort chosen - effort unchosen
            
            % Note: possible to compute net value (R-E) and do the same for it
            
        end % task
        
        % pool both tasks together
        % money chosen - money unchosen
            if (sum(strcmp(reg_names.Ep, money_chosen_ConNm)) > 0 &&...
                    sum(strcmp(reg_names.Ep, money_unchosen_ConNm)) > 0) &&...
                    (sum(strcmp(reg_names.Em, money_chosen_ConNm)) > 0 &&...
                    sum(strcmp(reg_names.Em, money_unchosen_ConNm)) > 0)
                
                % extract contrast
                moneyCh_Ep_idx = strcmp(con_names,['Ep ',money_chosen_ConNm]);
                moneyUnch_Ep_idx = strcmp(con_names,['Ep ',money_unchosen_ConNm]);
                moneyCh_Em_idx = strcmp(con_names,['Em ',money_chosen_ConNm]);
                moneyUnch_Em_idx = strcmp(con_names,['Em ',money_unchosen_ConNm]);
                con_moneyCh_min_moneyUnch_EpEm = con_vector(moneyCh_Ep_idx, :) - con_vector(moneyUnch_Ep_idx, :) +...
                    con_vector(moneyCh_Em_idx, :) - con_vector(moneyUnch_Em_idx, :);
                % positive contrast
                jReg = jReg + 1;
                con_names{jReg} = ['Ep+Em REG ',taskPhase_nm,' ',RP_nm,': money chosen - unchosen'];
                con_vector(jReg, 1:n_totalRegs) = con_moneyCh_min_moneyUnch_EpEm;
                
                % negative contrast
                jReg = jReg + 1;
                con_names{jReg} = ['Ep+Em REG ',taskPhase_nm,' ',RP_nm,': money unchosen - chosen'];
                con_vector(jReg, 1:n_totalRegs) = -con_moneyCh_min_moneyUnch_EpEm;
            end % money chosen - money unchosen
            
            % effort chosen - effort unchosen
            if (sum(strcmp(reg_names.Ep, effort_chosen_ConNm)) > 0 &&...
                    sum(strcmp(reg_names.Ep, effort_unchosen_ConNm)) > 0) &&...
                    (sum(strcmp(reg_names.Em, effort_chosen_ConNm)) > 0 &&...
                    sum(strcmp(reg_names.Em, effort_unchosen_ConNm)) > 0)
                % extract contrast
                effortCh_Ep_idx = strcmp(con_names,['Ep ',effort_chosen_ConNm]);
                effortUnch_Ep_idx = strcmp(con_names,['Ep ',effort_unchosen_ConNm]);
                effortCh_Em_idx = strcmp(con_names,['Em ',effort_chosen_ConNm]);
                effortUnch_Em_idx = strcmp(con_names,['Em ',effort_unchosen_ConNm]);
                con_effortCh_min_effortUnch_EpEm = con_vector(effortCh_Ep_idx, :) - con_vector(effortUnch_Ep_idx, :) +...
                    con_vector(effortCh_Em_idx, :) - con_vector(effortUnch_Em_idx, :);
                % positive contrast
                jReg = jReg + 1;
                con_names{jReg} = ['Ep+Em REG ',taskPhase_nm,' ',RP_nm,': effort chosen - unchosen'];
                con_vector(jReg, 1:n_totalRegs) = con_effortCh_min_effortUnch_EpEm;
                
                % negative contrast
                jReg = jReg + 1;
                con_names{jReg} = ['Ep+Em REG ',taskPhase_nm,' ',RP_nm,': effort unchosen - chosen'];
                con_vector(jReg, 1:n_totalRegs) = -con_effortCh_min_effortUnch_EpEm;
            end % effort chosen - effort unchosen
        
        % Note: possible to compute net value (R-E) and do the same for it
    end % RP
end % task phase

end % function