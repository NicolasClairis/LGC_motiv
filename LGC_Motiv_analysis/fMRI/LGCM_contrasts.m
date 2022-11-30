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
subj_analysis_folder    = [subj_folder, filesep, 'fMRI_analysis' filesep,...
    'functional',filesep,'preproc_sm_',num2str(preproc_sm_kernel),'mm',filesep];
[resultsFolderName] = fMRI_subFolder(subj_analysis_folder, GLM, condition);

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
n_totalRegs = n_regsPerTask.Ep*(runs.nb_runs.Ep) +...
    n_regsPerTask.Em*(runs.nb_runs.Em) +...
    runs.nb_runs.Ep + runs.nb_runs.Em;

% then compare that n_regs is equal to the total number of betas output
betaFiles = ls([resultsFolderName,'beta_*.nii']);
n_1stLevelBetas = size(betaFiles,1);
if n_1stLevelBetas ~= n_totalRegs
    error(['Problem: number of regressors predicted by contrasts = ',num2str(n_totalRegs),...
        ' while number of betas produced by 1st level = ',num2str(n_1stLevelBetas),...
        ' for subject ',sub_nm]);
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
% same for Effort trials if have been modelled separately
timePeriod = {'choice','chosen','Eperf','fbk'};
n_timePeriods = length(timePeriod);
timePeriod_names = {'choice','chosen','effort','feedback'};
regTypes = {'ONSET','REG'};
nRegTypes = length(regTypes);
tasks = {'Ep','Em'};
nTasks = length(tasks);
% loop through time periods
for iTimePeriod = 1:n_timePeriods
    timePhase_nm = timePeriod{iTimePeriod};
    timePhase_nm_bis = timePeriod_names{iTimePeriod};
    
    % loop through tasks (including tasks pooled)
    for iTask = 1:nTasks
        task_nm = tasks{iTask};
        
        % check effort conditions
        if (ismember(task_nm,{'Ep','Em'})) ||...
                (GLMprm.(timePhase_nm).Ep.splitPerE == GLMprm.(timePhase_nm).Em.splitPerE)
            switch GLMprm.(timePhase_nm).Ep.splitPerE
                case 0
                    Econd = {'E'};
                case 1
                    Econd = {'E1','E2','E3'};
                case 2
                    Econd = {'Ech0','Ech1','Ech2','Ech3'};
                case 3
                    Econd = {'lEch','hEch'};
            end
        else
            error(['case not ready yet where split of effort levels ',...
                'differs between mental and physical effort']);
        end
        
        for iEcond = 1:length(Econd)
            Econd_nm = Econd{iEcond};
            
            % check if R and P trials have been split, if not no need to pool
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
                                regRExpression_short = [task_nm,' ',posNeg,'ONSET ',timePhase_nm_bis,' R'];
                                regPExpression_short = [task_nm,' ',posNeg,'ONSET ',timePhase_nm_bis,' P'];
                                regRExpression = [regRExpression_short,' ',Econd_nm];
                                regPExpression = [regPExpression_short,' ',Econd_nm];
                                areThereRregs = strcmp(con_names, regRExpression);
                                areTherePregs = strcmp(con_names, regPExpression);
                            case 'REG'
                                regRExpression_short = [task_nm,' ',posNeg,regType,' ',timePhase_nm_bis,' R'];
                                regPExpression_short = [task_nm,' ',posNeg,regType,' ',timePhase_nm_bis,' P'];
                                regRExpression = [regRExpression_short,' ',Econd_nm,':'];
                                regPExpression = [regPExpression_short,' ',Econd_nm,':'];
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
                                
                                %% pool the two contrasts
                                con_vec_RP_tmp = con_vector(jRregList(iRPreg), :) +...
                                    con_vector(jPregList(iRPreg), :);
                                
                                % positive contrast (no need for negative contrast
                                % for once, because will automatically be created
                                % by pooling negative R and P together)
                                jReg = jReg + 1;
                                switch regType
                                    case 'ONSET'
                                        con_names{jReg} = [regRExpression_short,'P ',Econd_nm];
                                    case 'REG'
                                        con_names{jReg} = [regRExpression_short,'P ',Econd_nm,':',conR_nm];
                                end
                                con_vector(jReg, 1:n_totalRegs) = con_vec_RP_tmp;
                                
                                %% compare R and P also
                                con_vec_RminP_tmp = con_vector(jRregList(iRPreg), :) -...
                                    con_vector(jPregList(iRPreg), :);
                                jReg = jReg + 1;
                                switch regType
                                    case 'ONSET'
                                        con_names{jReg} = [regRExpression_short,'-P ',Econd_nm];
                                    case 'REG'
                                        con_names{jReg} = [regRExpression_short,'-P ',Econd_nm,':',conR_nm];
                                end
                                con_vector(jReg, 1:n_totalRegs) = con_vec_RminP_tmp;
                                
                            end % loop through regressors to pool
                        end % filter: more than 0 regressors then ok
                    end % positive/negative contrast
                end % ONSET/REG
            end % RPpool
        end % effort condition loop
    end % task loop
end % time loop

%% pool regressors depending on effort condition
tasks = {'Ep','Em','Ep+Em'};
nTasks = length(tasks);
timePeriod = {'choice','chosen','Eperf','fbk'};
n_timePeriods = length(timePeriod);
timePeriod_names = {'choice','chosen','effort','feedback'};
regTypes = {'ONSET','REG'};
nRegTypes = length(regTypes);
RP_cond = {'R','P','RP'};
% loop through time periods
for iTimePeriod = 1:n_timePeriods
    timePhase_nm = timePeriod{iTimePeriod};
    timePhase_nm_bis = timePeriod_names{iTimePeriod};
    
    % loop through tasks (including tasks pooled)
    for iTask = 1:nTasks
        task_nm = tasks{iTask};
        
        for iRP = 1:length(RP_cond) % loop through R/P/RP conditions
            RP_nm = RP_cond{iRP};
            % check if E trials have been split, if yes, perform the pooling
            if (ismember(task_nm,{'Ep','Em'}) && GLMprm.(timePhase_nm).(task_nm).splitPerE > 0) ||...
                    (GLMprm.(timePhase_nm).Ep.splitPerE > 0 && GLMprm.(timePhase_nm).Em.splitPerE > 0) % check if more than 1 effort condition modeled
                if (ismember(task_nm,{'Ep','Em'}) && GLMprm.(timePhase_nm).(task_nm).splitPerE == 1) ||...
                        (GLMprm.(timePhase_nm).Ep.splitPerE == 1 && GLMprm.(timePhase_nm).Em.splitPerE == 1) % split per effort proposed (3 levels)
                    
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
                            % to avoid pooling E conditions if already done, add ':' for REG
                            % specifically (no need for ONSET as full expression is
                            % already ok)
                            switch regType
                                case 'ONSET'
                                    regE1xpression = [task_nm,' ',posNeg,'ONSET ',timePhase_nm_bis,' ',RP_nm,' E1'];
                                    regE2xpression = [task_nm,' ',posNeg,'ONSET ',timePhase_nm_bis,' ',RP_nm,' E2'];
                                    regE3xpression = [task_nm,' ',posNeg,'ONSET ',timePhase_nm_bis,' ',RP_nm,' E3'];
                                    areThereE1regs = strcmp(con_names, regE1xpression);
                                    areThereE2regs = strcmp(con_names, regE2xpression);
                                    areThereE3regs = strcmp(con_names, regE3xpression);
                                case 'REG'
                                    regE1Expression = [task_nm,' ',posNeg,regType,' ',timePhase_nm_bis,' ',RP_nm,' E1:'];
                                    regE2Expression = [task_nm,' ',posNeg,regType,' ',timePhase_nm_bis,' ',RP_nm,' E2:'];
                                    regE3Expression = [task_nm,' ',posNeg,regType,' ',timePhase_nm_bis,' ',RP_nm,' E3:'];
                                    n_regNameSize = length(regE1Expression);
                                    areThereE1regs = strncmp(con_names, regE1xpression, n_regNameSize);
                                    areThereE2regs = strncmp(con_names, regE2Expression, n_regNameSize);
                                    areThereE3regs = strncmp(con_names, regE3Expression, n_regNameSize);
                            end
                            
                            n_Eregs = sum(areThereE1regs);
                            if n_Eregs > 0
                                jE1regList = find(areThereE1regs);
                                jE2regList = find(areThereE2regs);
                                jE3regList = find(areThereE3regs);
                                
                                % loop through regressors to pool
                                for iEreg = 1:n_Eregs
                                    E1reg_full_nm = con_names{jE1regList(iEreg)};
                                    E2reg_full_nm = con_names{jE2regList(iEreg)};
                                    E3reg_full_nm = con_names{jE3regList(iEreg)};
                                    if strcmp(regType,'REG')
                                        conE1_nm = E1reg_full_nm(n_regNameSize+1:end);
                                        conE2_nm = E2reg_full_nm(n_regNameSize+1:end);
                                        conE3_nm = E3reg_full_nm(n_regNameSize+1:end);
                                        if ~strcmp(conE1_nm, conE2_nm) || ~strcmp(conE1_nm, conE3_nm) || ~strcmp(conE2_nm, conE3_nm)
                                            error('problem mismatch E1/E2/E3 regressor for pooling. Please fix it.');
                                        end % comparison R/P regressor
                                    end
                                    
                                    % pool the contrasts
                                    con_vec_E_tmp = con_vector(jE1regList(iEreg), :) +...
                                        con_vector(jE2regList(iEreg), :) +...
                                        con_vector(jE3regList(iEreg), :);
                                    
                                    % positive contrast (no need for negative contrast
                                    % for once, because will automatically be created
                                    % by pooling negative R and P together)
                                    jReg = jReg + 1;
                                    switch regType
                                        case 'ONSET'
                                            con_names{jReg} = regE1xpression(1:(end-1));
                                        case 'REG'
                                            con_names{jReg} = [regE1xpression(1:(end-1)),':',conE1_nm];
                                    end
                                    con_vector(jReg, 1:n_totalRegs) = con_vec_E_tmp;
                                    
                                end % loop through regressors to pool
                            end % filter: more than 0 regressors then ok
                        end % positive/negative contrast
                    end % ONSET/REG
                    
                elseif (ismember(task_nm,{'Ep','Em'}) && GLMprm.(timePhase_nm).(task_nm).splitPerE == 2) ||...
                        (GLMprm.(timePhase_nm).Ep.splitPerE == 2 && GLMprm.(timePhase_nm).Em.splitPerE == 2) % split per effort chosen (4 levels)
                    
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
                            % to avoid pooling E conditions if already done, add ':' for REG
                            % specifically (no need for ONSET as full expression is
                            % already ok)
                            switch regType
                                case 'ONSET'
                                    regEch0Expression = [task_nm,' ',posNeg,'ONSET ',timePhase_nm_bis,' ',RP_nm,' Ech0'];
                                    regEch1Expression = [task_nm,' ',posNeg,'ONSET ',timePhase_nm_bis,' ',RP_nm,' Ech1'];
                                    regEch2Expression = [task_nm,' ',posNeg,'ONSET ',timePhase_nm_bis,' ',RP_nm,' Ech2'];
                                    regEch3Expression = [task_nm,' ',posNeg,'ONSET ',timePhase_nm_bis,' ',RP_nm,' Ech3'];
                                    areThereEch0regs = strcmp(con_names, regEch0Expression);
                                    areThereEch1regs = strcmp(con_names, regEch1Expression);
                                    areThereEch2regs = strcmp(con_names, regEch2Expression);
                                    areThereEch3regs = strcmp(con_names, regEch3Expression);
                                case 'REG'
                                    regEch0Expression = [task_nm,' ',posNeg,regType,' ',timePhase_nm_bis,' ',RP_nm,' Ech0:'];
                                    regEch1Expression = [task_nm,' ',posNeg,regType,' ',timePhase_nm_bis,' ',RP_nm,' Ech1:'];
                                    regEch2Expression = [task_nm,' ',posNeg,regType,' ',timePhase_nm_bis,' ',RP_nm,' Ech2:'];
                                    regEch3Expression = [task_nm,' ',posNeg,regType,' ',timePhase_nm_bis,' ',RP_nm,' Ech3:'];
                                    n_regNameSize = length(regEch0Expression);
                                    areThereEch0regs = strncmp(con_names, regEch0Expression, n_regNameSize);
                                    areThereEch1regs = strncmp(con_names, regEch1Expression, n_regNameSize);
                                    areThereEch2regs = strncmp(con_names, regEch2Expression, n_regNameSize);
                                    areThereEch3regs = strncmp(con_names, regEch3Expression, n_regNameSize);
                            end
                            
                            n_Eregs = sum(areThereE0regs);
                            if n_Eregs > 0
                                jEch0regList = find(areThereEch0regs);
                                jEch1regList = find(areThereEch1regs);
                                jEch2regList = find(areThereEch2regs);
                                jEch3regList = find(areThereEch3regs);
                                
                                % loop through regressors to pool
                                for iEreg = 1:n_Eregs
                                    Ech0reg_full_nm = con_names{jEch0regList(iEreg)};
                                    Ech1reg_full_nm = con_names{jEch1regList(iEreg)};
                                    Ech2reg_full_nm = con_names{jEch2regList(iEreg)};
                                    Ech3reg_full_nm = con_names{jEch3regList(iEreg)};
                                    if strcmp(regType,'REG')
                                        conEch0_nm = Ech0reg_full_nm(n_regNameSize+1:end);
                                        conEch1_nm = Ech1reg_full_nm(n_regNameSize+1:end);
                                        conEch2_nm = Ech2reg_full_nm(n_regNameSize+1:end);
                                        conEch3_nm = Ech3reg_full_nm(n_regNameSize+1:end);
                                        if ~strcmp(conEch0_nm, conEch1_nm) || ~strcmp(conEch0_nm, conEch2_nm) ||...
                                                ~strcmp(conEch0_nm, conEch3_nm) || ~strcmp(conEch1_nm, conEch2_nm) ||...
                                                ~strcmp(conEch1_nm, conEch3_nm) || ~strcmp(conEch2_nm, conEch3_nm)
                                            error('problem mismatch Ech0/Ech1/Ech2/Ech3 regressor for pooling. Please fix it.');
                                        end % comparison R/P regressor
                                    end
                                    
                                    % pool the contrasts
                                    con_vec_E_tmp = con_vector(jEch0regList(iEreg), :) +...
                                        con_vector(jEch1regList(iEreg), :) +...
                                        con_vector(jEch2regList(iEreg), :) +...
                                        con_vector(jEch3regList(iEreg), :);
                                    
                                    % positive contrast (no need for negative contrast
                                    % for once, because will automatically be created
                                    % by pooling negative R and P together)
                                    jReg = jReg + 1;
                                    switch regType
                                        case 'ONSET'
                                            con_names{jReg} = regEch0Expression(1:(end-1));
                                        case 'REG'
                                            con_names{jReg} = [regEch0Expression(1:(end-1)),':',conEch0_nm];
                                    end
                                    con_vector(jReg, 1:n_totalRegs) = con_vec_E_tmp;
                                    
                                end % loop through regressors to pool
                            end % filter: more than 0 regressors then ok
                        end % positive/negative contrast
                    end % ONSET/REG
                    
                elseif (ismember(task_nm,{'Ep','Em'}) && GLMprm.(timePhase_nm).(task_nm).splitPerE == 3) ||...
                        (GLMprm.(timePhase_nm).Ep.splitPerE == 3 && GLMprm.(timePhase_nm).Em.splitPerE == 3) % split per low vs high effort chosen
                    
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
                            % to avoid pooling E conditions if already done, add ':' for REG
                            % specifically (no need for ONSET as full expression is
                            % already ok)
                            switch regType
                                case 'ONSET'
                                    regLowEchxpression = [task_nm,' ',posNeg,'ONSET ',timePhase_nm_bis,' ',RP_nm,' lEch'];
                                    regHighEchxpression = [task_nm,' ',posNeg,'ONSET ',timePhase_nm_bis,' ',RP_nm,' hEch'];
                                    areThereLowEchregs = strcmp(con_names, regLowEchxpression);
                                    areThereHighEchregs = strcmp(con_names, regHighEchxpression);
                                case 'REG'
                                    regLowEchxpression = [task_nm,' ',posNeg,regType,' ',timePhase_nm_bis,' ',RP_nm,' lEch:'];
                                    regHighEchxpression = [task_nm,' ',posNeg,regType,' ',timePhase_nm_bis,' ',RP_nm,' hEch:'];
                                    n_regNameSize = length(regLowEchxpression);
                                    areThereLowEchregs = strncmp(con_names, regLowEchxpression, n_regNameSize);
                                    areThereHighEchregs = strncmp(con_names, regHighEchxpression, n_regNameSize);
                            end
                            
                            n_Eregs = sum(areThereLowEchregs);
                            if n_Eregs > 0
                                jLowEchRegList = find(areThereLowEchregs);
                                jHighEchRegList = find(areThereHighEchregs);
                                
                                % loop through regressors to pool
                                for iEreg = 1:n_Eregs
                                    lowEchReg_full_nm = con_names{jLowEchRegList(iEreg)};
                                    highEchReg_full_nm = con_names{jHighEchRegList(iEreg)};
                                    if strcmp(regType,'REG')
                                        conLowEch_nm = lowEchReg_full_nm(n_regNameSize+1:end);
                                        conHighEch_nm = highEchReg_full_nm(n_regNameSize+1:end);
                                        if ~strcmp(conLowEch_nm, conHighEch_nm)
                                            error('problem mismatch low E vs high E chosen regressor for pooling. Please fix it.');
                                        end % comparison R/P regressor
                                    end
                                    
                                    %% pool the contrasts
                                    con_vec_E_tmp = con_vector(jLowEchRegList(iEreg), :) +...
                                        con_vector(jHighEchRegList(iEreg), :);
                                    
                                    % positive contrast (no need for negative contrast
                                    % for once, because will automatically be created
                                    % by pooling negative R and P together)
                                    jReg = jReg + 1;
                                    switch regType
                                        case 'ONSET'
                                            con_names{jReg} = [regLowEchxpression(1:(end-4)),'E'];
                                        case 'REG'
                                            con_names{jReg} = [regLowEchxpression(1:(end-5)),'E:',conLowEch_nm];
                                    end
                                    con_vector(jReg, 1:n_totalRegs) = con_vec_E_tmp;
                                    
                                    %% compare low vs high effort choice also
                                    con_vec_HighEchMinLowEch_tmp = -con_vector(jLowEchRegList(iEreg), :) +...
                                        con_vector(jHighEchRegList(iEreg), :);
                                    jReg = jReg + 1;
                                    switch regType
                                        case 'ONSET'
                                            con_names{jReg} = [regHighEchxpression,' - lEch'];
                                        case 'REG'
                                            con_names{jReg} = [regHighEchxpression(1:(end-1)),' - lEch',':',conLowEch_nm];
                                    end
                                    con_vector(jReg, 1:n_totalRegs) = con_vec_HighEchMinLowEch_tmp;
                                    
                                end % loop through regressors to pool
                            end % filter: more than 0 regressors then ok
                        end % positive/negative contrast
                    end % ONSET/REG
                    
                end % effort pool
            end % effort to be pooled
        end % RP loop
    end % task loop
end % time period

%% compare Vchosen - Vunchosen if Vch and Vunch have been modelled separately
Epm = {'Ep','Em'};
RP_cond = {'R','P','RP'};
Esplit_cond = {'E','E1','E2','E3',...
    'Ech0','Ech1','Ech2','Ech3',...
    'lEch','hEch'};
task_phases = {'choice','chosen','effort'};
for iTaskPhase = 1:length(task_phases)
    taskPhase_nm = task_phases{iTaskPhase};
    for iRP = 1:length(RP_cond)
        RP_nm = RP_cond{iRP};
        
        for iE = length(Esplit_cond)
            Esplit_nm = Esplit_cond{iE};
            
            % money chosen - money unchosen
            money_chosen_ConNm = ['REG ',taskPhase_nm,' ',RP_nm,' ',Esplit_nm,': money chosen'];
            money_unchosen_ConNm = ['REG ',taskPhase_nm,' ',RP_nm,' ',Esplit_nm,': money unchosen'];
            % effort chosen - effort unchosen
            effort_chosen_ConNm = ['REG ',taskPhase_nm,' ',RP_nm,' ',Esplit_nm,': effort chosen'];
            effort_unchosen_ConNm = ['REG ',taskPhase_nm,' ',RP_nm,' ',Esplit_nm,': effort unchosen'];
            
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
                    con_names{jReg} = [task_nm,' REG ',taskPhase_nm,' ',RP_nm,' ',Esplit_nm,': money chosen - unchosen'];
                    con_vector(jReg, 1:n_totalRegs) = con_moneyCh_min_moneyUnch;
                    
                    % negative contrast
                    jReg = jReg + 1;
                    con_names{jReg} = [task_nm,' REG ',taskPhase_nm,' ',RP_nm,' ',Esplit_nm,': money unchosen - chosen'];
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
                    con_names{jReg} = [task_nm,' REG ',taskPhase_nm,' ',RP_nm,' ',Esplit_nm,': effort chosen - unchosen'];
                    con_vector(jReg, 1:n_totalRegs) = con_effortCh_min_effortUnch;
                    
                    % negative contrast
                    jReg = jReg + 1;
                    con_names{jReg} = [task_nm,' REG ',taskPhase_nm,' ',RP_nm,' ',Esplit_nm,': effort unchosen - chosen'];
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
                con_names{jReg} = ['Ep+Em REG ',taskPhase_nm,' ',RP_nm,' ',Esplit_nm,': money chosen - unchosen'];
                con_vector(jReg, 1:n_totalRegs) = con_moneyCh_min_moneyUnch_EpEm;
                
                % negative contrast
                jReg = jReg + 1;
                con_names{jReg} = ['Ep+Em REG ',taskPhase_nm,' ',RP_nm,' ',Esplit_nm,': money unchosen - chosen'];
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
                con_names{jReg} = ['Ep+Em REG ',taskPhase_nm,' ',RP_nm,' ',Esplit_nm,': effort chosen - unchosen'];
                con_vector(jReg, 1:n_totalRegs) = con_effortCh_min_effortUnch_EpEm;
                
                % negative contrast
                jReg = jReg + 1;
                con_names{jReg} = ['Ep+Em REG ',taskPhase_nm,' ',RP_nm,' ',Esplit_nm,': effort unchosen - chosen'];
                con_vector(jReg, 1:n_totalRegs) = -con_effortCh_min_effortUnch_EpEm;
            end % effort chosen - effort unchosen
            
            % Note: possible to compute net value (R-E) and do the same for it
        end % Effort loop
    end % RP
end % task phase

end % function