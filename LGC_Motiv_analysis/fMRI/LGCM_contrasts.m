function[con_names, con_vector] = LGCM_contrasts(study_nm, sub_nm, GLM,...
    computer_root, preproc_sm_kernel, condition, biasFieldCorr)
% [con_names, con_vector] = LGCM_contrasts(study_nm, sub_nm, GLM,...
%   computer_root, preproc_sm_kernel, condition, biasFieldCorr)
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
% biasFieldCorr: binary variable indicating whether results are based on
% bias-field corrected files (1) or not (0)
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
switch biasFieldCorr
    case 0
        subj_analysis_folder    = [subj_folder, filesep, 'fMRI_analysis' filesep,...
            'functional',filesep,'preproc_sm_',num2str(preproc_sm_kernel),'mm',filesep];
    case 1
        subj_analysis_folder    = [subj_folder, filesep, 'fMRI_analysis' filesep,...
            'functional',filesep,'preproc_sm_',num2str(preproc_sm_kernel),'mm_with_BiasFieldCorrection',filesep];
    otherwise
        error('biasFieldCorr not defined in the input. Please enter a value');
end
[resultsFolderName] = fMRI_subFolder(subj_analysis_folder, GLM, condition);

%% extract GLM informations
dispGLM = 0;
% disp(['GLM',num2str(GLM),' launching contrasts']);
[reg_names, n_regsPerTask] = GLM_details(GLM, dispGLM);
GLMprm = which_GLM(GLM);

%% define runs based on current subject
runs = runs_definition(study_nm, sub_nm, condition);

%% general parameters
timePeriod = {'choice','chosen','Eperf','fbk'};
timePeriod_names = {'choice','chosen','effort','feedback'};
n_timePeriods = length(timePeriod);
regTypes = {'ONSET','REG'};
nRegTypes = length(regTypes);
tasks = {'Ep','Em'};
nTasks = length(tasks);
tasks_bis = {'Ep','Em','Ep+Em'};
nTasks_bis = length(tasks_bis);
RP_cond = {'R','P','RP'};
nRP = length(RP_cond);
Esplit_cond = {'E','E1','E2','E3',...
    'Ech0','Ech1','Ech2','Ech3',...
    'lEch','hEch'};
nE_conds = length(Esplit_cond);

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
% else
%     disp(['CID',sub_nm,' number of regressors and number of betas produced in 1st level - ok'])
end

%% extract run informations (to know which run corresponds to which task)
Ep_runs = strcmp(runs.tasks,'Ep');
Em_runs = strcmp(runs.tasks,'Em');

% extract number of runs per task to normalize the data by the number of
% runs:
n_Ep_runs = sum(Ep_runs);
n_Em_runs = sum(Em_runs);

%% extract the contrasts of interest
jReg = 0;
%% physical effort
if runs.nb_runs.Ep > 0
    for iReg_Ep = 1:n_regsPerTask.Ep
        reg_nm = reg_names.Ep{iReg_Ep};
        
        if ~strcmp(reg_nm,'movement') && ~isempty(reg_nm) % ignore movement regressors (labelled as 'movement') and temporal derivative regressors (empty name)
            
            % run 1
            if Ep_runs(1) == 1 && Em_runs(1) == 0
                run1_Ep_vec = strcmp(reg_names.Ep, reg_nm);
            elseif Em_runs(1) == 1 && Ep_runs(1) == 0
                run1_Ep_vec = zeros(1,n_regsPerTask.Em);
            else
                error('problem with run 1');
            end
            
            % run 2
            if runs.nb_runs.Ep + runs.nb_runs.Em >= 2
                if Ep_runs(2) == 1 && Em_runs(2) == 0
                    run2_Ep_vec = strcmp(reg_names.Ep, reg_nm);
                elseif Em_runs(2) == 1 && Ep_runs(2) == 0
                    run2_Ep_vec = zeros(1,n_regsPerTask.Em);
                else
                    error('problem with run 2');
                end
            else
                run2_Ep_vec = [];
            end
            
            % run 3
            if runs.nb_runs.Ep + runs.nb_runs.Em >= 3
                if Ep_runs(3) == 1 && Em_runs(3) == 0
                    run3_Ep_vec =  strcmp(reg_names.Ep, reg_nm);
                elseif Em_runs(3) == 1 && Ep_runs(3) == 0
                    run3_Ep_vec = zeros(1,n_regsPerTask.Em);
                else
                    error('problem with run 3');
                end
            else
                run3_Ep_vec = [];
            end
            
            % run 4
            if runs.nb_runs.Ep + runs.nb_runs.Em == 4
                if Ep_runs(4) == 1 && Em_runs(4) == 0
                    run4_Ep_vec =  strcmp(reg_names.Ep, reg_nm);
                elseif Em_runs(4) == 1 && Ep_runs(4) == 0
                    run4_Ep_vec = zeros(1,n_regsPerTask.Em);
                else
                    error('problem with run 4');
                end
            else
                run4_Ep_vec = [];
            end
            
            %% global variable
            con_vec_Ep_tmp1 = [run1_Ep_vec,...
                run2_Ep_vec,...
                run3_Ep_vec,...
                run4_Ep_vec,...
                zeros(1, runs.nb_runs.Ep), zeros(1, runs.nb_runs.Em)];
            %% normalize by number of Ep runs
            con_vec_Ep_tmp1 = con_vec_Ep_tmp1./n_Ep_runs;
            
            %% positive contrast
            jReg = jReg + 1;
            con_names{jReg} = ['Ep ',reg_nm];
            con_vector(jReg, 1:n_totalRegs) = con_vec_Ep_tmp1;
            
            %% negative contrast
            jReg = jReg + 1;
            con_names{jReg} = ['Ep -',reg_nm];
            con_vector(jReg, 1:n_totalRegs) = -con_vec_Ep_tmp1;
            
        end % movement and derivative filter
    end % loop through physical effort regressors
end % physical effort performed at least in one full run?

%% mental effort
if runs.nb_runs.Em > 0
    for iReg_Em = 1:n_regsPerTask.Em
        reg_nm = reg_names.Em{iReg_Em};
        
        if ~strcmp(reg_nm,'movement') && ~isempty(reg_nm) % ignore movement regressors (labelled as 'movement') and temporal derivative regressors (empty name)
            
            % run 1
            if Em_runs(1) == 1 && Ep_runs(1) == 0
                run1_Em_vec = strcmp(reg_names.Em, reg_nm);
            elseif Ep_runs(1) == 1 && Em_runs(1) == 0
                run1_Em_vec = zeros(1,n_regsPerTask.Ep);
            else
                error('problem with run 1');
            end
            
            % run 2
            if runs.nb_runs.Ep + runs.nb_runs.Em >= 2
                if Em_runs(2) == 1 && Ep_runs(2) == 0
                    run2_Em_vec = strcmp(reg_names.Em, reg_nm);
                elseif Ep_runs(2) == 1 && Em_runs(2) == 0
                    run2_Em_vec = zeros(1,n_regsPerTask.Ep);
                else
                    error('problem with run 2');
                end
            else
                run2_Em_vec = [];
            end
            
            % run 3
            if runs.nb_runs.Ep + runs.nb_runs.Em >= 3
                if Em_runs(3) == 1 && Ep_runs(3) == 0
                    run3_Em_vec =  strcmp(reg_names.Em, reg_nm);
                elseif Ep_runs(3) == 1 && Em_runs(3) == 0
                    run3_Em_vec = zeros(1,n_regsPerTask.Ep);
                else
                    error('problem with run 3');
                end
            else
                run3_Em_vec = [];
            end
            
            % run 4
            if runs.nb_runs.Ep + runs.nb_runs.Em == 4
                if Em_runs(4) == 1 && Ep_runs(4) == 0
                    run4_Em_vec =  strcmp(reg_names.Em, reg_nm);
                elseif Ep_runs(4) == 1 && Em_runs(4) == 0
                    run4_Em_vec = zeros(1,n_regsPerTask.Ep);
                else
                    error('problem with run 4');
                end
            else
                run4_Em_vec = [];
            end
            
            %% global variable
            con_vec_Em_tmp1 = [run1_Em_vec,...
                        run2_Em_vec,...
                        run3_Em_vec,...
                        run4_Em_vec,...
                        zeros(1, runs.nb_runs.Ep), zeros(1, runs.nb_runs.Em)];
            %% normalize by number of Em runs
            con_vec_Em_tmp1 = con_vec_Em_tmp1./n_Em_runs;
            
            %% positive contrast
            jReg = jReg + 1;
            con_names{jReg} = ['Em ',reg_nm];
            con_vector(jReg, 1:n_totalRegs) = con_vec_Em_tmp1;
            
            %% negative contrast
            jReg = jReg + 1;
            con_names{jReg} = ['Em -',reg_nm];
            con_vector(jReg, 1:n_totalRegs) = -con_vec_Em_tmp1;
            
        end % movement and derivative filter
    end % loop through mental effort regressors
end % mental effort performed at least in one full run?

%% pool Reward and Punishment regressors if modeled as separate regressors (but on the same onset)
for iTP = 1:n_timePeriods
    period_nm1 = timePeriod{iTP};
    period_nm1_bis = timePeriod_names{iTP};
    if ismember(period_nm1, {'choice','chosen','Eperf'})
        for iTask = 1:nTasks
            task_nm1 = tasks{iTask}; % Ep/Em
            for iRP = 1:nRP
                RP_nm1 = RP_cond{iRP};
                for iE = 1:nE_conds
                    Econd_nm1 = Esplit_cond{iE};
                    
                    % R/P variable option (high effort)
                    if GLMprm.(period_nm1).(task_nm1).(RP_nm1).(Econd_nm1).R_varOption ~= 0 &&...
                            GLMprm.(period_nm1).(task_nm1).(RP_nm1).(Econd_nm1).P_varOption ~= 0
                        if GLMprm.(period_nm1).(task_nm1).(RP_nm1).(Econd_nm1).R_varOption == 1 &&...
                                GLMprm.(period_nm1).(task_nm1).(RP_nm1).(Econd_nm1).P_varOption == 1
                            jReg_Rvar = strcmp(con_names,[task_nm1,' REG ',period_nm1_bis,' ',RP_nm1,' ',Econd_nm1,': R amount highE']);
                            jReg_Pvar = strcmp(con_names,[task_nm1,' REG ',period_nm1_bis,' ',RP_nm1,' ',Econd_nm1,': P amount highE']);
                            RPvar_name = ['REG ',period_nm1_bis,' ',RP_nm1,' ',Econd_nm1,': R+P amount highE'];
                        elseif GLMprm.(period_nm1).(task_nm1).(RP_nm1).(Econd_nm1).R_varOption == 2 &&...
                                GLMprm.(period_nm1).(task_nm1).(RP_nm1).(Econd_nm1).P_varOption == 2
                            jReg_Rvar = strcmp(con_names,[task_nm1,' REG ',period_nm1_bis,' ',RP_nm1,' ',Econd_nm1,': R level highE']);
                            jReg_Pvar = strcmp(con_names,[task_nm1,' REG ',period_nm1_bis,' ',RP_nm1,' ',Econd_nm1,': P level highE']);
                            RPvar_name = ['REG ',period_nm1_bis,' ',RP_nm1,' ',Econd_nm1,': R+P level highE'];
                        end % incentive variable option
                        % extract each task
                        con_vec_Rvar_tmp = con_vector(jReg_Rvar, :);
                        con_vec_Pvar_tmp = con_vector(jReg_Pvar, :);
                        % pool the two contrasts
                        con_vec_RP_tmp = (con_vec_Rvar_tmp + con_vec_Pvar_tmp)./2;
                        
                        
                        % positive contrast
                        jReg = jReg + 1;
                        con_names{jReg} = [task_nm1,' ',RPvar_name];
                        con_vector(jReg, 1:n_totalRegs) = con_vec_RP_tmp;
                        
                        % negative contrast
                        jReg = jReg + 1;
                        con_names{jReg} = [task_nm1,' -',RPvar_name];
                        con_vector(jReg, 1:n_totalRegs) = -con_vec_RP_tmp;
                    end
                    
                    % R/P chosen
                    if GLMprm.(period_nm1).(task_nm1).(RP_nm1).(Econd_nm1).R_chosen ~= 0 &&...
                            GLMprm.(period_nm1).(task_nm1).(RP_nm1).(Econd_nm1).P_chosen ~= 0
                        if (GLMprm.(period_nm1).(task_nm1).(RP_nm1).(Econd_nm1).R_chosen == 1 &&...
                                GLMprm.(period_nm1).(task_nm1).(RP_nm1).(Econd_nm1).P_chosen == 1) ||...
                                (GLMprm.(period_nm1).(task_nm1).(RP_nm1).(Econd_nm1).R_chosen == 3 &&...
                                GLMprm.(period_nm1).(task_nm1).(RP_nm1).(Econd_nm1).P_chosen == 3)
                            jReg_Rch = strcmp(con_names,[task_nm1,' REG ',period_nm1_bis,' ',RP_nm1,' ',Econd_nm1,': R amount chosen']);
                            jReg_Pch = strcmp(con_names,[task_nm1,' REG ',period_nm1_bis,' ',RP_nm1,' ',Econd_nm1,': P amount chosen']);
                            RPvar_name = ['REG ',period_nm1_bis,' ',RP_nm1,' ',Econd_nm1,': R+P amount chosen'];
                        elseif (GLMprm.(period_nm1).(task_nm1).(RP_nm1).(Econd_nm1).R_chosen == 2 &&...
                                GLMprm.(period_nm1).(task_nm1).(RP_nm1).(Econd_nm1).P_chosen == 2) ||...
                                (GLMprm.(period_nm1).(task_nm1).(RP_nm1).(Econd_nm1).R_chosen == 4 &&...
                                GLMprm.(period_nm1).(task_nm1).(RP_nm1).(Econd_nm1).P_chosen == 4)
                            jReg_Rch = strcmp(con_names,[task_nm1,' REG ',period_nm1_bis,' ',RP_nm1,' ',Econd_nm1,': R level chosen']);
                            jReg_Pch = strcmp(con_names,[task_nm1,' REG ',period_nm1_bis,' ',RP_nm1,' ',Econd_nm1,': P level chosen']);
                            RPvar_name = ['REG ',period_nm1_bis,' ',RP_nm1,' ',Econd_nm1,': R+P level chosen'];
                        else
                            error(['R_chosen = ',num2str(GLMprm.(period_nm1).(task_nm1).(RP_nm1).(Econd_nm1).R_chosen),...
                                ' and P_chosen = ',num2str(GLMprm.(period_nm1).(task_nm1).(RP_nm1).(Econd_nm1).P_chosen),...
                                ' not ready yet in LGCM_contrasts.m. Please fix.']);
                        end % incentive variable option
                        % extract each task
                        con_vec_Rch_tmp = con_vector(jReg_Rch, :);
                        con_vec_Pch_tmp = con_vector(jReg_Pch, :);
                        % pool the two contrasts
                        con_vec_RPch_tmp = (con_vec_Rch_tmp + con_vec_Pch_tmp)./2;
                        
                        % positive contrast
                        jReg = jReg + 1;
                        con_names{jReg} = [task_nm1,' ',RPvar_name];
                        con_vector(jReg, 1:n_totalRegs) = con_vec_RPch_tmp;
                        
                        % negative contrast
                        jReg = jReg + 1;
                        con_names{jReg} = [task_nm1,' -',RPvar_name];
                        con_vector(jReg, 1:n_totalRegs) = -con_vec_RPch_tmp;
                    end % incentive chosen option
                    
                    % incentive*effort variable option
                    if GLMprm.(period_nm1).(task_nm1).(RP_nm1).(Econd_nm1).R_level_x_E_varOption ~= 0 &&...
                            GLMprm.(period_nm1).(task_nm1).(RP_nm1).(Econd_nm1).P_level_x_E_varOption ~= 0
                        if GLMprm.(period_nm1).(task_nm1).(RP_nm1).(Econd_nm1).R_level_x_E_varOption == 1 &&...
                                GLMprm.(period_nm1).(task_nm1).(RP_nm1).(Econd_nm1).P_level_x_E_varOption == 1
                            jReg_REvar = strcmp(con_names,[task_nm1,' REG ',period_nm1_bis,' ',RP_nm1,' ',Econd_nm1,': R x E non-default']);
                            jReg_PEvar = strcmp(con_names,[task_nm1,' REG ',period_nm1_bis,' ',RP_nm1,' ',Econd_nm1,': P x E non-default']);
                            REPEvar_name = ['REG ',period_nm1_bis,' ',RP_nm1,' ',Econd_nm1,': (R x E)+(P x E) non-default'];
                        end % incentive variable option
                        % extract each task
                        con_vec_REvar_tmp = con_vector(jReg_REvar, :);
                        con_vec_PEvar_tmp = con_vector(jReg_PEvar, :);
                        % pool the two contrasts
                        con_vec_REPEvar_tmp = (con_vec_REvar_tmp + con_vec_PEvar_tmp)./2;
                        
                        % positive contrast
                        jReg = jReg + 1;
                        con_names{jReg} = [task_nm1,' ',REPEvar_name];
                        con_vector(jReg, 1:n_totalRegs) = con_vec_REPEvar_tmp;
                        
                        % negative contrast
                        jReg = jReg + 1;
                        con_names{jReg} = [task_nm1,' -',REPEvar_name];
                        con_vector(jReg, 1:n_totalRegs) = -con_vec_REPEvar_tmp;
                    end % incentive*effort variable option
                    
                    % incentive*effort chosen option
                    if GLMprm.(period_nm1).(task_nm1).(RP_nm1).(Econd_nm1).R_level_x_E_chosen ~= 0 &&...
                            GLMprm.(period_nm1).(task_nm1).(RP_nm1).(Econd_nm1).P_level_x_E_chosen ~= 0
                        if GLMprm.(period_nm1).(task_nm1).(RP_nm1).(Econd_nm1).R_level_x_E_chosen == 1 &&...
                                GLMprm.(period_nm1).(task_nm1).(RP_nm1).(Econd_nm1).P_level_x_E_chosen == 1
                            jReg_REch = strcmp(con_names,[task_nm1,' REG ',period_nm1_bis,' ',RP_nm1,' ',Econd_nm1,': R x E chosen']);
                            jReg_PEch = strcmp(con_names,[task_nm1,' REG ',period_nm1_bis,' ',RP_nm1,' ',Econd_nm1,': P x E chosen']);
                            REPEch_name = ['REG ',period_nm1_bis,' ',RP_nm1,' ',Econd_nm1,': (R x E)+(P x E) chosen'];
                        end % incentive variable option
                        % extract each task
                        con_vec_REch_tmp = con_vector(jReg_REch, :);
                        con_vec_PEch_tmp = con_vector(jReg_PEch, :);
                        % pool the two contrasts
                        con_vec_REPEch_tmp = (con_vec_REch_tmp + con_vec_PEch_tmp)./2;
                        
                        % positive contrast
                        jReg = jReg + 1;
                        con_names{jReg} = [task_nm1,' ',REPEch_name];
                        con_vector(jReg, 1:n_totalRegs) = con_vec_REPEch_tmp;
                        
                        % negative contrast
                        jReg = jReg + 1;
                        con_names{jReg} = [task_nm1,' -',REPEch_name];
                        con_vector(jReg, 1:n_totalRegs) = -con_vec_REPEch_tmp;
                    end % incentive*effort chosen option
                end % loop over effort conditions
            end % loop over R/P/RP
        end % filter relevant periods
    end % loop over task
end % loop over trial period
%% pool contrasts across both tasks
con_names_v0 = con_names; % need to create a second variable otherwise the script will use the contrasts created within the loop as well
if runs.nb_runs.Em > 0 && runs.nb_runs.Ep > 0
    
    for iReg = 1:size(con_names_v0,2) % loop through all regressors
        if strcmp(con_names_v0{iReg}(1:3),'Ep ') && ~strcmp(con_names_v0{iReg}(4),'-') % avoid negative contrasts cause that would be redundant with creating the negative contrast
            reg_nm = con_names_v0{iReg}(4:end);
            if (sum(strcmp(con_names_v0,['Em ',reg_nm])) > 0) % check that regressor also exists for mental
                    
                % extract previously defined contrast for each task
                jReg_Ep = strcmp(con_names_v0,['Ep ',reg_nm]);
                jReg_Em = strcmp(con_names_v0,['Em ',reg_nm]);
                % extract each task
                con_vec_Ep_tmp2 = con_vector(jReg_Ep, :);
                con_vec_Em_tmp2 = con_vector(jReg_Em, :);
                % pool the two contrasts
                con_vec_EpEm_tmp = (con_vec_Ep_tmp2 + con_vec_Em_tmp2)./2;
                
                % positive contrast
                jReg = jReg + 1;
                con_names{jReg} = ['Ep+Em ',reg_nm];
                con_vector(jReg, 1:n_totalRegs) = con_vec_EpEm_tmp;
                
                % negative contrast
                jReg = jReg + 1;
                con_names{jReg} = ['Ep+Em -',reg_nm];
                con_vector(jReg, 1:n_totalRegs) = -con_vec_EpEm_tmp;
            elseif ((length(reg_nm) > 14) && strcmp(reg_nm((end-14):end),'effort integral')) ||...
                    ((length(reg_nm) > 24) && strcmp(reg_nm((end-24):end),'effort integral overshoot')) % pool effort between two tasks
                
                % extract physical effort perf index
                jReg_Ep_perf = strcmp(con_names_v0,['Ep ',reg_nm]);
                
                % extract all infos to use for Em (RP/R/P, E/E1/E2/E3/etc.)
                reg_infos = reg_nm(1:(find(reg_nm==':')));
                
                % extract mental effort performance now
                jReg_Em_perf = strcmp(con_names_v0,['Em ',reg_infos,' efficacy']);
                
                % extract each task
                con_vec_Ep_tmp2 = con_vector(jReg_Ep_perf, :);
                con_vec_Em_tmp2 = con_vector(jReg_Em_perf, :);
                
                % pool the two contrasts
                con_vec_EpEm_tmp = (con_vec_Ep_tmp2 + con_vec_Em_tmp2)./2;
                % create new regressor name
                reg_nm_bis = [reg_infos,' F_integral x efficacy'];
                
                % positive contrast
                jReg = jReg + 1;
                con_names{jReg} = ['Ep+Em ',reg_nm_bis];
                con_vector(jReg, 1:n_totalRegs) = con_vec_EpEm_tmp;
                
                % negative contrast
                jReg = jReg + 1;
                con_names{jReg} = ['Ep+Em -',reg_nm_bis];
                con_vector(jReg, 1:n_totalRegs) = -con_vec_EpEm_tmp;
            end % check if contrast exists for both tasks (+ ignore movement and derivative regressors)
        end % Ep filter in order to do that only once (on the Ep regressors) as Em regressors should be redundant
    end % loop through previous contrasts
end % at least one run of each task has been performed?

%% pool contrasts across R and P trials if have been modelled separately
% loop through time periods
for iTimePeriod = 1:n_timePeriods
    period_nm3 = timePeriod{iTimePeriod};
    period_nm3_bis = timePeriod_names{iTimePeriod};
    
    % loop through tasks (including tasks pooled)
    for iTask = 1:nTasks
        task_nm3 = tasks{iTask}; % Ep/Em
        
        % check effort conditions
        if (ismember(task_nm3,{'Ep','Em'})) ||...
                (GLMprm.(period_nm3).Ep.splitPerE == GLMprm.(period_nm3).Em.splitPerE)
            switch GLMprm.(period_nm3).Ep.splitPerE
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
        
        % loop through E conditions
        for iEcond = 1:length(Econd)
            Econd_nm3 = Econd{iEcond};
            
            %% check if R and P trials have been split in both tasks, if not no need to pool
            if (ismember(task_nm3,{'Ep','Em'}) && GLMprm.(period_nm3).(task_nm3).RPpool == 0) ||...
                    (GLMprm.(period_nm3).Ep.RPpool == 0 && GLMprm.(period_nm3).Em.RPpool == 0)
                % pool through onsets and regressors
                for iRegType = 1:nRegTypes
                    regType1 = regTypes{iRegType};
                    
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
                        switch regType1
                            case 'ONSET'
                                regRExpression_short = [task_nm3,' ',posNeg,'ONSET ',period_nm3_bis,' R'];
                                regPExpression_short = [task_nm3,' ',posNeg,'ONSET ',period_nm3_bis,' P'];
                                regRExpression = [regRExpression_short,' ',Econd_nm3];
                                regPExpression = [regPExpression_short,' ',Econd_nm3];
                                areThereRregs = strcmp(con_names, regRExpression);
                                areTherePregs = strcmp(con_names, regPExpression);
                            case 'REG'
                                regRExpression_short = [task_nm3,' ',posNeg,'REG ',period_nm3_bis,' R'];
                                regPExpression_short = [task_nm3,' ',posNeg,'REG ',period_nm3_bis,' P'];
                                regRExpression = [regRExpression_short,' ',Econd_nm3,':'];
                                regPExpression = [regPExpression_short,' ',Econd_nm3,':'];
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
                                if strcmp(regType1,'REG')
                                    conR_nm = Rreg_full_nm(n_regNameSize+1:end);
                                    conP_nm = Preg_full_nm(n_regNameSize+1:end);
                                    if ~strcmp(conR_nm, conP_nm) &&...
                                            ~strcmp(conR_nm,' R level chosen') &&...
                                            ~strcmp(conR_nm,' P level chosen')
                                        error('problem mismatch R and P regressor for pooling. Please fix it.');
                                    end % comparison R/P regressor
                                end
                                
                                %% extract R and P
                                con_R_tmp = con_vector(jRregList(iRPreg), :);
                                con_P_tmp = con_vector(jPregList(iRPreg), :);
                                
                                %% pool the two contrasts
                                con_vec_RP_tmp = (con_R_tmp + con_P_tmp)./2;
                                
                                % positive contrast (no need for negative contrast
                                % for once, because will automatically be created
                                % by pooling negative R and P together)
                                jReg = jReg + 1;
                                switch regType1
                                    case 'ONSET'
                                        con_names{jReg} = [regRExpression_short,'P ',Econd_nm3];
                                    case 'REG'
                                        con_names{jReg} = [regRExpression_short,'P ',Econd_nm3,':',conR_nm];
                                end
                                con_vector(jReg, 1:n_totalRegs) = con_vec_RP_tmp;
                                
                                %% compare R and P also
                                con_vec_RminP_tmp = con_R_tmp - con_P_tmp;
                                jReg = jReg + 1;
                                switch regType1
                                    case 'ONSET'
                                        con_names{jReg} = [regRExpression_short,'-P ',Econd_nm3];
                                    case 'REG'
                                        con_names{jReg} = [regRExpression_short,'-P ',Econd_nm3,':',conR_nm];
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
% loop through time periods
for iTimePeriod = 1:n_timePeriods
    period_nm4 = timePeriod{iTimePeriod};
    period_nm4_bis = timePeriod_names{iTimePeriod};
    
    % loop through tasks (including tasks pooled)
    for iTask = 1:nTasks_bis
        task_nm_bis = tasks_bis{iTask}; % includes Ep/Em/Ep+Em
        
        % loop through R/P/RP conditions
        for iRP = 1:length(RP_cond)
            RP_nm3 = RP_cond{iRP};
            
            % check if E trials have been split, if yes, perform the pooling
            if (ismember(task_nm_bis,{'Ep','Em'}) && GLMprm.(period_nm4).(task_nm_bis).splitPerE > 0) ||...
                    (GLMprm.(period_nm4).Ep.splitPerE > 0 && GLMprm.(period_nm4).Em.splitPerE > 0) % check if more than 1 effort condition modeled
                
                %% split per effort proposed (3 levels)
                if (ismember(task_nm_bis,{'Ep','Em'}) && GLMprm.(period_nm4).(task_nm_bis).splitPerE == 1) ||...
                        (GLMprm.(period_nm4).Ep.splitPerE == 1 && GLMprm.(period_nm4).Em.splitPerE == 1)
                    
                    % pool through onsets and regressors
                    for iRegType = 1:nRegTypes
                        regType1 = regTypes{iRegType};
                        
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
                            switch regType1
                                case 'ONSET'
                                    regE1Expression = [task_nm_bis,' ',posNeg,'ONSET ',period_nm4_bis,' ',RP_nm3,' E1'];
                                    regE2Expression = [task_nm_bis,' ',posNeg,'ONSET ',period_nm4_bis,' ',RP_nm3,' E2'];
                                    regE3Expression = [task_nm_bis,' ',posNeg,'ONSET ',period_nm4_bis,' ',RP_nm3,' E3'];
                                    areThereE1regs = strcmp(con_names, regE1Expression);
                                    areThereE2regs = strcmp(con_names, regE2Expression);
                                    areThereE3regs = strcmp(con_names, regE3Expression);
                                case 'REG'
                                    regE1Expression = [task_nm_bis,' ',posNeg,'REG ',period_nm4_bis,' ',RP_nm3,' E1'];
                                    regE2Expression = [task_nm_bis,' ',posNeg,'REG ',period_nm4_bis,' ',RP_nm3,' E2'];
                                    regE3Expression = [task_nm_bis,' ',posNeg,'REG ',period_nm4_bis,' ',RP_nm3,' E3'];
                                    n_regNameSize = length(regE1Expression);
                                    areThereE1regs = strncmp(con_names, regE1Expression, n_regNameSize);
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
                                    if strcmp(regType1,'REG')
                                        conE1_nm = E1reg_full_nm(n_regNameSize+1:end);
                                        conE2_nm = E2reg_full_nm(n_regNameSize+1:end);
                                        conE3_nm = E3reg_full_nm(n_regNameSize+1:end);
                                        if ~strcmp(conE1_nm, conE2_nm) || ~strcmp(conE1_nm, conE3_nm) || ~strcmp(conE2_nm, conE3_nm)
                                            error('problem mismatch E1/E2/E3 regressor for pooling. Please fix it.');
                                        end % comparison R/P regressor
                                    end
                                    
                                    %% extract individual contrasts
                                    con_E1_tmp = con_vector(jE1regList(iEreg), :);
                                    con_E2_tmp = con_vector(jE2regList(iEreg), :);
                                    con_E3_tmp = con_vector(jE3regList(iEreg), :);
                                    
                                    %% pool the contrasts
                                    con_vec_E_tmp1 = (con_E1_tmp + con_E2_tmp + con_E3_tmp)./3;
                                    
                                    % positive contrast (no need for negative contrast
                                    % for once, because will automatically be created
                                    % by pooling negative R and P together)
                                    jReg = jReg + 1;
                                    switch regType1
                                        case 'ONSET'
                                            con_names{jReg} = regE1Expression(1:(end-1));
                                        case 'REG'
                                            con_names{jReg} = [regE1Expression(1:(end-1)),':',conE1_nm];
                                    end
                                    con_vector(jReg, 1:n_totalRegs) = con_vec_E_tmp1;
                                    
                                end % loop through regressors to pool
                            end % filter: more than 0 regressors then ok
                        end % positive/negative contrast
                    end % ONSET/REG
                    
                    %% split per effort chosen (4 levels)
                elseif (ismember(task_nm_bis,{'Ep','Em'}) && GLMprm.(period_nm4).(task_nm_bis).splitPerE == 2) ||...
                        (GLMprm.(period_nm4).Ep.splitPerE == 2 && GLMprm.(period_nm4).Em.splitPerE == 2)
                    
                    % pool through onsets and regressors
                    for iRegType = 1:nRegTypes
                        regType1 = regTypes{iRegType};
                        
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
                            switch regType1
                                case 'ONSET'
                                    regEch0Expression = [task_nm_bis,' ',posNeg,'ONSET ',period_nm4_bis,' ',RP_nm3,' Ech0'];
                                    regEch1Expression = [task_nm_bis,' ',posNeg,'ONSET ',period_nm4_bis,' ',RP_nm3,' Ech1'];
                                    regEch2Expression = [task_nm_bis,' ',posNeg,'ONSET ',period_nm4_bis,' ',RP_nm3,' Ech2'];
                                    regEch3Expression = [task_nm_bis,' ',posNeg,'ONSET ',period_nm4_bis,' ',RP_nm3,' Ech3'];
                                    areThereEch0regs = strcmp(con_names, regEch0Expression);
                                    areThereEch1regs = strcmp(con_names, regEch1Expression);
                                    areThereEch2regs = strcmp(con_names, regEch2Expression);
                                    areThereEch3regs = strcmp(con_names, regEch3Expression);
                                case 'REG'
                                    regEch0Expression = [task_nm_bis,' ',posNeg,'REG ',period_nm4_bis,' ',RP_nm3,' Ech0'];
                                    regEch1Expression = [task_nm_bis,' ',posNeg,'REG ',period_nm4_bis,' ',RP_nm3,' Ech1'];
                                    regEch2Expression = [task_nm_bis,' ',posNeg,'REG ',period_nm4_bis,' ',RP_nm3,' Ech2'];
                                    regEch3Expression = [task_nm_bis,' ',posNeg,'REG ',period_nm4_bis,' ',RP_nm3,' Ech3'];
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
                                    if strcmp(regType1,'REG')
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
                                    
                                    %% extract individual contrasts
                                    con_Ech0_tmp = con_vector(jEch0regList(iEreg), :);
                                    con_Ech1_tmp = con_vector(jEch1regList(iEreg), :);
                                    con_Ech2_tmp = con_vector(jEch2regList(iEreg), :);
                                    con_Ech3_tmp = con_vector(jEch3regList(iEreg), :);
                                    
                                    %% pool the contrasts
                                    con_vec_E_tmp2 = (con_Ech0_tmp + con_Ech1_tmp + con_Ech2_tmp + con_Ech3_tmp)./4;
                                    
                                    % positive contrast (no need for negative contrast
                                    % for once, because will automatically be created
                                    % by pooling negative R and P together)
                                    jReg = jReg + 1;
                                    switch regType1
                                        case 'ONSET'
                                            con_names{jReg} = regEch0Expression(1:(end-1));
                                        case 'REG'
                                            con_names{jReg} = [regEch0Expression(1:(end-1)),':',conEch0_nm];
                                    end
                                    con_vector(jReg, 1:n_totalRegs) = con_vec_E_tmp2;
                                    
                                end % loop through regressors to pool
                            end % filter: more than 0 regressors then ok
                        end % positive/negative contrast
                    end % ONSET/REG
                    
                    %% split per low vs high effort chosen
                elseif (ismember(task_nm_bis,{'Ep','Em'}) && GLMprm.(period_nm4).(task_nm_bis).splitPerE == 3) ||...
                        (GLMprm.(period_nm4).Ep.splitPerE == 3 && GLMprm.(period_nm4).Em.splitPerE == 3)
                    
                    % pool through onsets and regressors
                    for iRegType = 1:nRegTypes
                        regType1 = regTypes{iRegType};
                        
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
                            switch regType1
                                case 'ONSET'
                                    regLowEchxpression = [task_nm_bis,' ',posNeg,'ONSET ',period_nm4_bis,' ',RP_nm3,' lEch'];
                                    regHighEchxpression = [task_nm_bis,' ',posNeg,'ONSET ',period_nm4_bis,' ',RP_nm3,' hEch'];
                                    areThereLowEchregs = strcmp(con_names, regLowEchxpression);
                                    areThereHighEchregs = strcmp(con_names, regHighEchxpression);
                                case 'REG'
                                    regLowEchxpression = [task_nm_bis,' ',posNeg,'REG ',period_nm4_bis,' ',RP_nm3,' lEch'];
                                    regHighEchxpression = [task_nm_bis,' ',posNeg,'REG ',period_nm4_bis,' ',RP_nm3,' hEch'];
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
                                    if strcmp(regType1,'REG')
                                        conLowEch_nm = lowEchReg_full_nm(n_regNameSize+1:end);
                                        conHighEch_nm = highEchReg_full_nm(n_regNameSize+1:end);
                                        if ~strcmp(conLowEch_nm, conHighEch_nm)
                                            error('problem mismatch low E vs high E chosen regressor for pooling. Please fix it.');
                                        end % comparison R/P regressor
                                    end
                                    
                                    %% extract individual contrasts
                                    con_lowEch_tmp = con_vector(jLowEchRegList(iEreg),:);
                                    con_highEch_tmp = con_vector(jHighEchRegList(iEreg),:);
                                    
                                    %% pool the contrasts
                                    con_vec_E_tmp3 = (con_lowEch_tmp + con_highEch_tmp)./2;
                                    
                                    % positive contrast (no need for negative contrast
                                    % for once, because will automatically be created
                                    % by pooling negative R and P together)
                                    jReg = jReg + 1;
                                    switch regType1
                                        case 'ONSET'
                                            con_names{jReg} = [regLowEchxpression(1:(end-4)),'E'];
                                        case 'REG'
                                            con_names{jReg} = [regLowEchxpression(1:(end-5)),'E:',conLowEch_nm];
                                    end
                                    con_vector(jReg, 1:n_totalRegs) = con_vec_E_tmp3;
                                    
                                    %% compare low vs high effort choice also
                                    con_vec_HighEchMinLowEch_tmp = con_highEch_tmp - con_lowEch_tmp;
                                    jReg = jReg + 1;
                                    switch regType1
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
for iTaskPhase = 1:n_timePeriods
    period_nm5_bis = timePeriod_names{iTaskPhase};
    for iRP = 1:length(RP_cond)
        RP_nm3 = RP_cond{iRP};
        
        for iE = 1:nE_conds
            Econd_nm2 = Esplit_cond{iE};
            
            % money chosen - money unchosen
            money_chosen_ConNm = ['REG ',period_nm5_bis,' ',RP_nm3,' ',Econd_nm2,': money chosen'];
            money_unchosen_ConNm = ['REG ',period_nm5_bis,' ',RP_nm3,' ',Econd_nm2,': money unchosen'];
            % effort chosen - effort unchosen
            effort_chosen_ConNm = ['REG ',period_nm5_bis,' ',RP_nm3,' ',Econd_nm2,': effort chosen'];
            effort_unchosen_ConNm = ['REG ',period_nm5_bis,' ',RP_nm3,' ',Econd_nm2,': effort unchosen'];
            
            for iTask = 1:length(tasks)
                task_nm4 = tasks{iTask}; % Ep/Em
                
                %% money chosen - money unchosen
                if (sum(strcmp(reg_names.(task_nm4), money_chosen_ConNm)) > 0 &&...
                        sum(strcmp(reg_names.(task_nm4), money_unchosen_ConNm)) > 0)
                    
                    %% extract and normalise individual contrasts
                    moneyCh_idx = strcmp(con_names,[task_nm4,' ',money_chosen_ConNm]);
                    moneyUnch_idx = strcmp(con_names,[task_nm4,' ',money_unchosen_ConNm]);
                    con_moneyCh_tmp = con_vector(moneyCh_idx, :);
                    con_moneyUnch_tmp = con_vector(moneyUnch_idx, :);
                    %% check the difference
                    con_moneyCh_min_moneyUnch = con_moneyCh_tmp - con_moneyUnch_tmp;
                    
                    %% positive contrast
                    jReg = jReg + 1;
                    con_names{jReg} = [task_nm4,' REG ',period_nm5_bis,' ',RP_nm3,' ',Econd_nm2,': money chosen - unchosen'];
                    con_vector(jReg, 1:n_totalRegs) = con_moneyCh_min_moneyUnch;
                    
                    %% negative contrast
                    jReg = jReg + 1;
                    con_names{jReg} = [task_nm4,' REG ',period_nm5_bis,' ',RP_nm3,' ',Econd_nm2,': money unchosen - chosen'];
                    con_vector(jReg, 1:n_totalRegs) = -con_moneyCh_min_moneyUnch;
                end % money chosen - money unchosen
                
                %% effort chosen - effort unchosen
                if (sum(strcmp(reg_names.(task_nm4), effort_chosen_ConNm)) > 0 &&...
                        sum(strcmp(reg_names.(task_nm4), effort_unchosen_ConNm)) > 0)
                    %% extract individual contrasts
                    effortCh_idx = strcmp(con_names,[task_nm4,' ',effort_chosen_ConNm]);
                    effortUnch_idx = strcmp(con_names,[task_nm4,' ',effort_unchosen_ConNm]);
                    con_Ech_tmp = con_vector(effortCh_idx, :);
                    con_Eunch_tmp = con_vector(effortUnch_idx, :);
                    %% check the difference
                    con_effortCh_min_effortUnch = con_Ech_tmp - con_Eunch_tmp;
                    
                    %% positive contrast
                    jReg = jReg + 1;
                    con_names{jReg} = [task_nm4,' REG ',period_nm5_bis,' ',RP_nm3,' ',Econd_nm2,': effort chosen - unchosen'];
                    con_vector(jReg, 1:n_totalRegs) = con_effortCh_min_effortUnch;
                    
                    %% negative contrast
                    jReg = jReg + 1;
                    con_names{jReg} = [task_nm4,' REG ',period_nm5_bis,' ',RP_nm3,' ',Econd_nm2,': effort unchosen - chosen'];
                    con_vector(jReg, 1:n_totalRegs) = -con_effortCh_min_effortUnch;
                end % effort chosen - effort unchosen
                
                % Note: possible to compute net value (R-E) and do the same for it
                
            end % task
            
            %% pool both tasks together
            %% money chosen - money unchosen
            if (sum(strcmp(reg_names.Ep, money_chosen_ConNm)) > 0 &&...
                    sum(strcmp(reg_names.Ep, money_unchosen_ConNm)) > 0) &&...
                    (sum(strcmp(reg_names.Em, money_chosen_ConNm)) > 0 &&...
                    sum(strcmp(reg_names.Em, money_unchosen_ConNm)) > 0)
                
                %% extract individual contrasts
                moneyCh_Ep_idx = strcmp(con_names,['Ep ',money_chosen_ConNm]);
                moneyUnch_Ep_idx = strcmp(con_names,['Ep ',money_unchosen_ConNm]);
                moneyCh_Em_idx = strcmp(con_names,['Em ',money_chosen_ConNm]);
                moneyUnch_Em_idx = strcmp(con_names,['Em ',money_unchosen_ConNm]);
                
                con_moneyCh_Ep_tmp = con_vector(moneyCh_Ep_idx, :);
                con_moneyUnch_Ep_tmp = con_vector(moneyUnch_Ep_idx, :);
                con_moneyCh_min_Unch_Ep_tmp = con_moneyCh_Ep_tmp - con_moneyUnch_Ep_tmp;
                con_moneyCh_Em_tmp = con_vector(moneyCh_Em_idx, :);
                con_moneyUnch_Em_tmp = con_vector(moneyUnch_Em_idx, :);
                con_moneyCh_min_Unch_Em_tmp = con_moneyCh_Em_tmp - con_moneyUnch_Em_tmp;
                %% average two tasks
                con_moneyCh_min_moneyUnch_EpEm = (con_moneyCh_min_Unch_Ep_tmp + con_moneyCh_min_Unch_Em_tmp)./2;
                
                %% positive contrast
                jReg = jReg + 1;
                con_names{jReg} = ['Ep+Em REG ',period_nm5_bis,' ',RP_nm3,' ',Econd_nm2,': money chosen - unchosen'];
                con_vector(jReg, 1:n_totalRegs) = con_moneyCh_min_moneyUnch_EpEm;
                
                %% negative contrast
                jReg = jReg + 1;
                con_names{jReg} = ['Ep+Em REG ',period_nm5_bis,' ',RP_nm3,' ',Econd_nm2,': money unchosen - chosen'];
                con_vector(jReg, 1:n_totalRegs) = -con_moneyCh_min_moneyUnch_EpEm;
            end % money chosen - money unchosen
            
            %% effort chosen - effort unchosen
            if (sum(strcmp(reg_names.Ep, effort_chosen_ConNm)) > 0 &&...
                    sum(strcmp(reg_names.Ep, effort_unchosen_ConNm)) > 0) &&...
                    (sum(strcmp(reg_names.Em, effort_chosen_ConNm)) > 0 &&...
                    sum(strcmp(reg_names.Em, effort_unchosen_ConNm)) > 0)
                %% extract and normalise individual contrasts
                effortCh_Ep_idx = strcmp(con_names,['Ep ',effort_chosen_ConNm]);
                effortUnch_Ep_idx = strcmp(con_names,['Ep ',effort_unchosen_ConNm]);
                effortCh_Em_idx = strcmp(con_names,['Em ',effort_chosen_ConNm]);
                effortUnch_Em_idx = strcmp(con_names,['Em ',effort_unchosen_ConNm]);
                
                con_Ech_Ep_tmp = con_vector(effortCh_Ep_idx, :);
                con_Eunch_Ep_tmp = con_vector(effortUnch_Ep_idx, :);
                con_Ech_min_Eunch_Ep_tmp = con_Ech_Ep_tmp - con_Eunch_Ep_tmp;
                con_Ech_Em_tmp = con_vector(effortCh_Em_idx, :);
                con_Eunch_Em_tmp = con_vector(effortUnch_Em_idx, :);
                con_Ech_min_Eunch_Em_tmp = con_Ech_Em_tmp - con_Eunch_Em_tmp;
                con_effortCh_min_effortUnch_EpEm = (con_Ech_min_Eunch_Ep_tmp + con_Ech_min_Eunch_Em_tmp)./2;
                
                %% positive contrast
                jReg = jReg + 1;
                con_names{jReg} = ['Ep+Em REG ',period_nm5_bis,' ',RP_nm3,' ',Econd_nm2,': effort chosen - unchosen'];
                con_vector(jReg, 1:n_totalRegs) = con_effortCh_min_effortUnch_EpEm;
                
                %% negative contrast
                jReg = jReg + 1;
                con_names{jReg} = ['Ep+Em REG ',period_nm5_bis,' ',RP_nm3,' ',Econd_nm2,': effort unchosen - chosen'];
                con_vector(jReg, 1:n_totalRegs) = -con_effortCh_min_effortUnch_EpEm;
            end % effort chosen - effort unchosen
            
            % Note: possible to compute net value (R-E) and do the same for it
        end % Effort loop
    end % RP
end % task phase

end % function