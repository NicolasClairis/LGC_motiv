function[con_names, con_vector] = LGCM_contrasts(study_nm, sub_nm, GLM, computer_root)
% [con_names, con_vector] = LGCM_contrasts(study_nm, sub_nm, GLM, computer_root)
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
% OUTPUTS
% con_names: list of contrast names
%
% con_vector: matrix with all contrasts organized in nb_contrasts (rows)*nb_regressors (columns)


%% working directories
switch study_nm
    case 'fMRI_pilots'
        root = fullfile(computer_root,'fMRI_pilots');
    case 'study1'
        root = fullfile(computer_root,'study1');
    case 'study2'
        root = fullfile(computer_root,'study2');
end
subj_folder             = [root, filesep, sub_nm];
subj_analysis_folder    = [subj_folder, filesep, 'fMRI_analysis' filesep];
resultsFolderName = [subj_analysis_folder 'functional', filesep,...
        'GLM',num2str(GLM),filesep];

%% extract GLM informations
[reg_names, n_regsPerTask] = GLM_details(GLM);

%% define runs based on current subject
runs = runs_definition(study_nm, sub_nm);

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
    disp([sub_nm,' number of regressors and number of betas produced in 1st level - ok'])
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
            
            % define basic contrasts per task (pool runs of the same task)
            if runs.nb_runs.Ep == 2 && runs.nb_runs.Em == 2 % default case = 2 Ep runs & 2 Em runs
                if Ep_runs(1) == 1 && Ep_runs(2) == 0 && Ep_runs(3) == 1 && Ep_runs(4) == 0 &&...
                        Em_runs(1) == 0 && Em_runs(2) == 1 && Em_runs(3) == 0 && Em_runs(4) == 1
                    % Ep/Em/Ep/Em
                    con_vec_Ep_tmp = [strcmp(reg_names.Ep, reg_nm),...
                        zeros(1,n_regsPerTask.Em),...
                        strcmp(reg_names.Ep, reg_nm),...
                        zeros(1,n_regsPerTask.Em),...
                        zeros(1, runs.nb_runs.Ep), zeros(1, runs.nb_runs.Em)] ;
                elseif Ep_runs(1) == 0 && Ep_runs(2) == 1 && Ep_runs(3) == 0 && Ep_runs(4) == 1 &&...
                        Em_runs(1) == 1 && Em_runs(2) == 0 && Em_runs(3) == 1 && Em_runs(4) == 0
                    % Em/Ep/Em/Ep
                    con_vec_Ep_tmp = [zeros(1,n_regsPerTask.Em),...
                        strcmp(reg_names.Ep, reg_nm),...
                        zeros(1,n_regsPerTask.Em),...
                        strcmp(reg_names.Ep, reg_nm),...
                        zeros(1, runs.nb_runs.Ep), zeros(1, runs.nb_runs.Em)] ;
                else
                    error('case not ready yet');
                end
                
            elseif runs.nb_runs.Ep == 1 && runs.nb_runs.Em == 1
                if Ep_runs(1) == 1 && Ep_runs(2) == 0 &&...
                        Em_runs(1) == 0 && Em_runs(2) == 1
                    % Ep/Em
                    con_vec_Ep_tmp = [strcmp(reg_names.Ep, reg_nm),...
                        zeros(1,n_regsPerTask.Em),...
                        zeros(1, runs.nb_runs.Ep), zeros(1, runs.nb_runs.Em)] ;
                elseif Ep_runs(1) == 0 && Ep_runs(2) == 1 &&...
                        Em_runs(1) == 1 && Em_runs(2) == 0
                    % Em/Ep
                    con_vec_Ep_tmp = [zeros(1,n_regsPerTask.Em),...
                        strcmp(reg_names.Ep, reg_nm),...
                        zeros(1, runs.nb_runs.Ep), zeros(1, runs.nb_runs.Em)] ;
                else
                    error('case not ready yet');
                end
                
            elseif runs.nb_runs.Ep == 1 && runs.nb_runs.Em == 0
                if Ep_runs(1) == 1
                    % Ep
                    con_vec_Ep_tmp = [strcmp(reg_names.Ep, reg_nm),...
                        zeros(1, runs.nb_runs.Ep)] ;
                else
                    error('case not ready yet');
                end
            else
                error('case not ready yet');
            end
            
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
            
            % define basic contrasts per task (pool runs of the same task)
            if runs.nb_runs.Ep == 2 && runs.nb_runs.Em == 2 % default case = 2 Ep runs & 2 Em runs
                if Ep_runs(1) == 1 && Ep_runs(2) == 0 && Ep_runs(3) == 1 && Ep_runs(4) == 0 &&...
                        Em_runs(1) == 0 && Em_runs(2) == 1 && Em_runs(3) == 0 && Em_runs(4) == 1
                    % Ep/Em/Ep/Em
                    con_vec_Em_tmp = [zeros(1,n_regsPerTask.Ep),...
                        strcmp(reg_names.Em, reg_nm),...
                        zeros(1,n_regsPerTask.Ep),...
                        strcmp(reg_names.Em, reg_nm),...
                        zeros(1, runs.nb_runs.Ep), zeros(1, runs.nb_runs.Em)] ;
                elseif Ep_runs(1) == 0 && Ep_runs(2) == 1 && Ep_runs(3) == 0 && Ep_runs(4) == 1 &&...
                        Em_runs(1) == 1 && Em_runs(2) == 0 && Em_runs(3) == 1 && Em_runs(4) == 0
                    % Em/Ep/Em/Ep
                    con_vec_Em_tmp = [strcmp(reg_names.Em, reg_nm),...
                        zeros(1,n_regsPerTask.Ep),...
                        strcmp(reg_names.Em, reg_nm),...
                        zeros(1,n_regsPerTask.Ep),...
                        zeros(1, runs.nb_runs.Ep), zeros(1, runs.nb_runs.Em)] ;
                else
                    error('case not ready yet');
                end
                
            elseif runs.nb_runs.Ep == 1 && runs.nb_runs.Em == 1
                if Ep_runs(1) == 1 && Ep_runs(2) == 0 &&...
                        Em_runs(1) == 0 && Em_runs(2) == 1
                    % Ep/Em
                    con_vec_Em_tmp = [zeros(1,n_regsPerTask.Ep),...
                        strcmp(reg_names.Em, reg_nm),...
                        zeros(1, runs.nb_runs.Ep), zeros(1, runs.nb_runs.Em)] ;
                elseif Ep_runs(1) == 0 && Ep_runs(2) == 1 &&...
                        Em_runs(1) == 1 && Em_runs(2) == 0
                    % Em/Ep
                    con_vec_Em_tmp = [strcmp(reg_names.Em, reg_nm),...
                        zeros(1,n_regsPerTask.Ep),...
                        zeros(1, runs.nb_runs.Ep), zeros(1, runs.nb_runs.Em)] ;
                else
                    error('case not ready yet');
                end
                
            elseif runs.nb_runs.Em == 1 && runs.nb_runs.Ep == 0
                if Em_runs(1) == 1
                    % Em
                    con_vec_Em_tmp = [strcmp(reg_names.Em, reg_nm),...
                        zeros(1, runs.nb_runs.Em)] ;
                else
                    error('case not ready yet');
                end
            else
                error('case not ready yet');
            end
            
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
        if sum(strcmp(reg_names.Em, reg_nm)) > 0 &&...
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

end % function