function [reg_names, n_regs] = GLM_details_DCM(GLM, DCM_mode, dispRegs, runs, n_runs)
% [reg_names, n_regs] = GLM_details_DCM(GLM, dispRegs, runs, n_runs)
% GLM_details_DCM will provide the number of regressors, the text with the
% details of the GLM, the order of the regressors in the current GLM for
% each task. This script is similar to GLM_details.m which is made for
% classic GLMs but GLM_details decomposes regressors for each session, while
% GLM_details_DCM will give the total list of regressors directly.
%
% INPUTS
% GLM: number of the GLM
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
% dispRegs: display information on regressors (1) or not (0)
%
% runs: structure with information regarding runs to keep for the current
% subject
%
% n_runs: number of runs to keep for the current subject
%
% OUTPUTS
% reg_names: structure with physical (Ep) and mental (Em) regressors for
% each run
% NB: careful to no include any '|' symbol in reg_names because that may
% interfere with further scripts => try to replace by 'abs()'.
%
% n_regs: structure with number of regressors for each physical (Ep) and
% each mental (Em) run

%% do not display parameters if dispRegs left empty (by default)
if ~exist('dispRegs','var') ||...
        isempty(dispRegs) ||...
        ~ismember(dispRegs,[0,1])
    dispRegs = 0;
end

%% load GLM parameters
GLMprm = which_GLM(GLM);
% check if onsets only GLM do not perform this
if GLMprm.gal.onsets_only == 1
    error('GLM_details.m is not ready for onsets_only ON');
end
% disp GLM number
dispRegFn(['GLM ',num2str(GLM)],dispRegs);

%% initialize the variables of interest
reg_names = deal( {} );
n_regs = 0;

%% main parameters
generalPrm = GLMprm.gal;
dispRegFn('** general parameters of the GLM **',dispRegs);

% grey matter mask
grey_mask = generalPrm.grey_mask;
switch grey_mask
    case 0
        dispRegFn('all voxels included (no grey matter filter)',dispRegs);
    case 1
        dispRegFn('Use of 1st level probability grey mask for each subject.',dispRegs);
    case 2
        dispRegFn('Use of 1st level probability grey mask for each subject.',dispRegs);
end

% temporal/spatial derivative
% check if derivative has been included => if so each regressor will be
% modelled as baseline + temporal and/or spatial derivative
add_drv = generalPrm.add_drv;
switch add_drv
    case 0
        dispRegFn('no derivative included',dispRegs);
    case 1
        dispRegFn('Use of temporal derivative => be aware that each regressor will be doubled in the matrix list.',dispRegs);
    case 2
        dispRegFn('Use of temporal and spatial derivative => be aware that each regressor will be tripled in the matrix list',dispRegs);
end

% zscore variables per run?
z_perRun = generalPrm.zPerRun;
switch z_perRun
    case 1
        dispRegFn('all variables are zscored per run',dispRegs);
end

% orthogonalize variables or not?
orth_vars = generalPrm.orth_vars;
switch orth_vars
    case 0
        dispRegFn('Variables are not orthogonalized',dispRegs);
    case 1
        dispRegFn('Variables are orthogonalized.',dispRegs);
end

% GLM parameters for choice Reward/Punishment pool
task_id_nm = 'Ep';
warning(['use task_id=''',task_id_nm,''' for extraction of GLMprm to avoid modifying which_GLM.']);
warning('If onsets/durations/regressors differ between Ep and Em, consider modifying which_GLM accordingly.');

%% load matlabbatch
% initialize conditions
tasks = {'Ep','Em'};
nTasks = length(tasks);

%% all fixation crosses (pool of pre-choice and pre-effort cross)
if ~strcmp(GLMprm.model_onset.(task_id_nm).allCrosses,'none')
    switch DCM_mode
        case 1 % all sessions independent
            for iRun = 1:n_runs
                run_nb_nm = ['r',num2str(runs.runsToKeep(iRun))];
                n_regs = n_regs + 1;
                reg_names{n_regs} = ['ONSET fixation cross ',run_nb_nm];
                dispRegFn([num2str(n_regs),') ONSET fixation cross ',run_nb_nm,': ',GLMprm.model_onset.(task_id_nm).allCrosses,' '],dispRegs);
                % if derivative added => add derivatives
                n_regs = n_regs + add_drv;
            end % run loop
        case 2 % all tasks independent but sessions pooled
            for iTask = 1:nTasks
                task_nm = tasks{iTask};
                n_regs = n_regs + 1;
                reg_names{n_regs} = ['ONSET fixation cross ',task_nm];
                dispRegFn([num2str(n_regs),') ONSET fixation cross ',task_nm,': ',GLMprm.model_onset.(task_id_nm).allCrosses,' '],dispRegs);
                % if derivative added => add derivatives
                n_regs = n_regs + add_drv;
            end % task loop
        case {3,4,5} % all pooled across sessions
            n_regs = n_regs + 1;
            reg_names{n_regs} = 'ONSET fixation cross';
            dispRegFn([num2str(n_regs),') ONSET fixation cross: ',GLMprm.model_onset.(task_id_nm).allCrosses,' '],dispRegs);
            % if derivative added => add derivatives
            n_regs = n_regs + add_drv;
    end % DCM_mode
end % all crosses

%% pre-choice fixation cross
if ~strcmp(GLMprm.model_onset.(task_id_nm).preChoiceCross,'none')
    switch DCM_mode
        case 1 % all sessions independent
            for iRun = 1:n_runs
                run_nb_nm = ['r',num2str(runs.runsToKeep(iRun))];
                
                % pre-choice fixation cross
                n_regs = n_regs + 1;
                reg_names{n_regs} = ['ONSET preChoice white fixation cross ',run_nb_nm];
                dispRegFn([num2str(n_regs),') ONSET preChoice white cross ',run_nb_nm,': ',GLMprm.model_onset.(task_id_nm).preChoiceCross,' '],dispRegs);
                % if derivative added => add derivatives
                n_regs = n_regs + add_drv;
                
                % RT (first regressor)
                switch GLMprm.preChoiceCross.(task_id_nm).RT
                    case 4
                        n_regs = n_regs + 1;
                        reg_names{n_regs} = ['REG preChoice Cross RP E: RT ',run_nb_nm];
                        dispRegFn([num2str(n_regs),') preChoice Cross ',run_nb_nm,' RT (raw) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs = n_regs + add_drv;
                    case 5
                        n_regs = n_regs + 1;
                        reg_names{n_regs} = ['REG preChoice Cross RP E: RT ',run_nb_nm];
                        dispRegFn([num2str(n_regs),') preChoice Cross ',run_nb_nm,' RT (zscored per run) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs = n_regs + add_drv;
                    case 6
                        n_regs = n_regs + 1;
                        reg_names{n_regs} = ['REG preChoice Cross RP E: RT ',run_nb_nm];
                        dispRegFn([num2str(n_regs),') preChoice Cross ',run_nb_nm,' RT (zscored per subject ie across all runs) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs = n_regs + add_drv;
                end
                
                % binary variable indicating when choice = high effort option
                switch GLMprm.preChoiceCross.(task_id_nm).choiceHighE
                    case 0
                    case {1,2}
                        n_regs = n_regs + 1;
                        reg_names{n_regs} = ['REG preChoice Cross RP E: choice = highE ',run_nb_nm];
                        dispRegFn([num2str(n_regs),') preChoice Cross ',run_nb_nm,' choice hE '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs = n_regs + add_drv;
                    otherwise
                        error('not ready yet');
                end
                
                % effort chosen
                switch GLMprm.preChoiceCross.(task_id_nm).E_chosen
                    case 0
                    case {1,4}
                        n_regs = n_regs + 1;
                        reg_names{n_regs} = ['REG preChoice Cross RP E: effort chosen ',run_nb_nm];
                        dispRegFn([num2str(n_regs),') preChoice Cross ',run_nb_nm,': effort chosen (levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs = n_regs + add_drv;
                    case 2
                        error('case not ready yet');
                    case 3
                        n_regs = n_regs + 1;
                        reg_names{n_regs} = ['REG preChoice Cross RP E: effort-by-choice ',run_nb_nm];
                        dispRegFn([num2str(n_regs),') preChoice Cross ',run_nb_nm,': effort-by-choice '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs = n_regs + add_drv;
                    otherwise
                        error('case not ready yet');
                end
                
                % RT (last regressor)
                switch GLMprm.preChoiceCross.(task_id_nm).RT
                    case 1
                        n_regs = n_regs + 1;
                        reg_names{n_regs} = ['REG preChoice Cross RP E: RT ',run_nb_nm];
                        dispRegFn([num2str(n_regs),') preChoice Cross ',run_nb_nm,': RT (raw) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs = n_regs + add_drv;
                    case 2
                        n_regs = n_regs + 1;
                        reg_names{n_regs} = ['REG preChoice Cross RP E: RT ',run_nb_nm];
                        dispRegFn([num2str(n_regs),') preChoice Cross ',run_nb_nm,': RT (zscored per run) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs = n_regs + add_drv;
                    case 3
                        n_regs = n_regs + 1;
                        reg_names{n_regs} = ['REG preChoice Cross RP E: RT ',run_nb_nm];
                        dispRegFn([num2str(n_regs),') preChoice Cross ',run_nb_nm,': RT (zscored per subject ie across all runs) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs = n_regs + add_drv;
                end
            end % run loop
        case 2 % all tasks independent but sessions pooled
            for iTask = 1:nTasks
                task_nm = tasks{iTask};
                
                % pre-choice fixation cross
                n_regs = n_regs + 1;
                reg_names{n_regs} = ['ONSET preChoice white fixation cross ',task_nm];
                dispRegFn([num2str(n_regs),') ONSET preChoice white cross ',task_nm,': ',GLMprm.model_onset.(task_id_nm).preChoiceCross,' '],dispRegs);
                % if derivative added => add derivatives
                n_regs = n_regs + add_drv;
                
                % RT (first regressor)
                switch GLMprm.preChoiceCross.(task_id_nm).RT
                    case 4
                        n_regs = n_regs + 1;
                        reg_names{n_regs} = ['REG preChoice Cross RP E: RT ',task_nm];
                        dispRegFn([num2str(n_regs),') preChoice Cross ',task_nm,' RT (raw) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs = n_regs + add_drv;
                    case 5
                        n_regs = n_regs + 1;
                        reg_names{n_regs} = ['REG preChoice Cross RP E: RT ',task_nm];
                        dispRegFn([num2str(n_regs),') preChoice Cross ',task_nm,' RT (zscored per run) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs = n_regs + add_drv;
                    case 6
                        n_regs = n_regs + 1;
                        reg_names{n_regs} = ['REG preChoice Cross RP E: RT ',task_nm];
                        dispRegFn([num2str(n_regs),') preChoice Cross ',task_nm,' RT (zscored per subject ie across all runs) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs = n_regs + add_drv;
                end
                
                % binary variable indicating when choice = high effort option
                switch GLMprm.preChoiceCross.(task_id_nm).choiceHighE
                    case 0
                    case {1,2}
                        n_regs = n_regs + 1;
                        reg_names{n_regs} = ['REG preChoice Cross RP E: choice = highE ',task_nm];
                        dispRegFn([num2str(n_regs),') preChoice Cross ',task_nm,' choice hE '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs = n_regs + add_drv;
                    otherwise
                        error('not ready yet');
                end
                
                % effort chosen
                switch GLMprm.preChoiceCross.(task_id_nm).E_chosen
                    case 0
                    case {1,4}
                        n_regs = n_regs + 1;
                        reg_names{n_regs} = ['REG preChoice Cross RP E: effort chosen ',task_nm];
                        dispRegFn([num2str(n_regs),') preChoice Cross ',task_nm,': effort chosen (levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs = n_regs + add_drv;
                    case 2
                        error('case not ready yet');
                    case 3
                        n_regs = n_regs + 1;
                        reg_names{n_regs} = ['REG preChoice Cross RP E: effort-by-choice ',task_nm];
                        dispRegFn([num2str(n_regs),') preChoice Cross ',task_nm,': effort-by-choice '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs = n_regs + add_drv;
                    otherwise
                        error('case not ready yet');
                end
                
                % RT (last regressor)
                switch GLMprm.preChoiceCross.(task_id_nm).RT
                    case 1
                        n_regs = n_regs + 1;
                        reg_names{n_regs} = ['REG preChoice Cross RP E: RT ',task_nm];
                        dispRegFn([num2str(n_regs),') preChoice Cross ',task_nm,': RT (raw) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs = n_regs + add_drv;
                    case 2
                        n_regs = n_regs + 1;
                        reg_names{n_regs} = ['REG preChoice Cross RP E: RT ',task_nm];
                        dispRegFn([num2str(n_regs),') preChoice Cross ',task_nm,': RT (zscored per run) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs = n_regs + add_drv;
                    case 3
                        n_regs = n_regs + 1;
                        reg_names{n_regs} = ['REG preChoice Cross RP E: RT ',task_nm];
                        dispRegFn([num2str(n_regs),') preChoice Cross ',task_nm,': RT (zscored per subject ie across all runs) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs = n_regs + add_drv;
                end
            end % task loop
        case {3,4,5} % all pooled across sessions
            % pre-choice fixation cross
            n_regs = n_regs + 1;
            reg_names{n_regs} = 'ONSET preChoice white fixation cross';
            dispRegFn([num2str(n_regs),') ONSET preChoice white cross: ',GLMprm.model_onset.(task_id_nm).preChoiceCross,' '],dispRegs);
            % if derivative added => add derivatives
            n_regs = n_regs + add_drv;
            
            % RT (first regressor)
            switch GLMprm.preChoiceCross.(task_id_nm).RT
                case 4
                    n_regs = n_regs + 1;
                    reg_names{n_regs} = 'REG preChoice Cross RP E: RT';
                    dispRegFn([num2str(n_regs),') preChoice Cross: RT (raw) '],dispRegs);
                    % if derivative added => add derivatives
                    n_regs = n_regs + add_drv;
                case 5
                    n_regs = n_regs + 1;
                    reg_names{n_regs} = 'REG preChoice Cross RP E: RT';
                    dispRegFn([num2str(n_regs),') preChoice Cross: RT (zscored per run) '],dispRegs);
                    % if derivative added => add derivatives
                    n_regs = n_regs + add_drv;
                case 6
                    n_regs = n_regs + 1;
                    reg_names{n_regs} = 'REG preChoice Cross RP E: RT';
                    dispRegFn([num2str(n_regs),') preChoice Cross: RT (zscored per subject ie across all runs) '],dispRegs);
                    % if derivative added => add derivatives
                    n_regs = n_regs + add_drv;
            end
            
            % binary variable indicating when choice = high effort option
            switch GLMprm.preChoiceCross.(task_id_nm).choiceHighE
                case 0
                case {1,2}
                    n_regs = n_regs + 1;
                    reg_names{n_regs} = 'REG preChoice Cross RP E: choice = highE';
                    dispRegFn([num2str(n_regs),') preChoice Cross: choice hE '],dispRegs);
                    % if derivative added => add derivatives
                    n_regs = n_regs + add_drv;
                otherwise
                    error('not ready yet');
            end
            
            % effort chosen
            switch GLMprm.preChoiceCross.(task_id_nm).E_chosen
                case 0
                case {1,4}
                    n_regs = n_regs + 1;
                    reg_names{n_regs} = 'REG preChoice Cross RP E: effort chosen';
                    dispRegFn([num2str(n_regs),') preChoice Cross: effort chosen (levels) '],dispRegs);
                    % if derivative added => add derivatives
                    n_regs = n_regs + add_drv;
                case 2
                    error('case not ready yet');
                case 3
                    n_regs = n_regs + 1;
                    reg_names{n_regs} = 'REG preChoice Cross RP E: effort-by-choice';
                    dispRegFn([num2str(n_regs),') preChoice Cross: effort-by-choice '],dispRegs);
                    % if derivative added => add derivatives
                    n_regs = n_regs + add_drv;
                otherwise
                    error('case not ready yet');
            end
            
            % RT (last regressor)
            switch GLMprm.preChoiceCross.(task_id_nm).RT
                case 1
                    n_regs = n_regs + 1;
                    reg_names{n_regs} = 'REG preChoice Cross RP E: RT';
                    dispRegFn([num2str(n_regs),') preChoice Cross: RT (raw) '],dispRegs);
                    % if derivative added => add derivatives
                    n_regs = n_regs + add_drv;
                case 2
                    n_regs = n_regs + 1;
                    reg_names{n_regs} = 'REG preChoice Cross RP E: RT';
                    dispRegFn([num2str(n_regs),') preChoice Cross: RT (zscored per run) '],dispRegs);
                    % if derivative added => add derivatives
                    n_regs = n_regs + add_drv;
                case 3
                    n_regs = n_regs + 1;
                    reg_names{n_regs} = 'REG preChoice Cross RP E: RT';
                    dispRegFn([num2str(n_regs),') preChoice Cross: RT (zscored per subject ie across all runs) '],dispRegs);
                    % if derivative added => add derivatives
                    n_regs = n_regs + add_drv;
            end
    end % DCM_mode
end % pre-choice cross

%% choice period
if ~strcmp(GLMprm.model_onset.(task_id_nm).choice,'none')
    % check if trials are split or not according to Reward and
    % Punishment
    if GLMprm.choice.(task_id_nm).RPpool == 1 % pool reward and punishment trials
        RP_dispChoice = {'RP'};
    elseif GLMprm.choice.(task_id_nm).RPpool == 0 % split reward and punishment
        RP_dispChoice = {'R','P'};
    end % RP pool
    n_RP_dispChoice = length(RP_dispChoice);
    
    % check if trials are split or not according to Effort levels
    switch GLMprm.choice.(task_id_nm).splitPerE
        case 0 % pool trials
            splitE_dispChoice = {'E'};
        case 1 % split according to effort proposed
            splitE_dispChoice = {'E1','E2','E3'};
        case 2 % split according to effort chosen
            splitE_dispChoice = {'Ech0','Ech1','Ech2','Ech3'};
        case 3 % split according to option chosen (low/high effort)
            splitE_dispChoice = {'lEch','hEch'};
    end % Effort level pool
    n_splitE_dispChoice = length(splitE_dispChoice);
    
    % loop through conditions for choice period
    for iRP_dispChoice = 1:n_RP_dispChoice
        RP_dispChoice_nm = RP_dispChoice{iRP_dispChoice};
        
        for iE_dispChoice = 1:n_splitE_dispChoice
            splitE_dispChoice_nm = splitE_dispChoice{iE_dispChoice};
            
            switch DCM_mode
                case 1 % all sessions independent
                    for iRun = 1:n_runs
                        run_nb_nm = ['r',num2str(runs.runsToKeep(iRun))];
                        
                        %% choice onset
                        n_regs = n_regs + 1;
                        reg_names{n_regs} = ['ONSET choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm];
                        dispRegFn([num2str(n_regs),') ONSET choice display options ',...
                            RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': ',GLMprm.model_onset.(task_id_nm).choice,' '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs = n_regs + add_drv;
                        
                        %% choice regressors
                        
                        % RT (first regressor)
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).RT
                            case {0,1,2,3}
                            case 4
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',...
                                    RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': RT'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': RT (raw) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 5
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',...
                                    RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': RT'];
                                dispRegFn([num2str(n_regs),') choice: RT (zscored per run) ',run_nb_nm],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 6
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',...
                                    RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': RT'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,' RT (zscored per subject ie across all runs) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % net value chosen option (first regressor)
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).NV_chosen
                            case {0,1,2,3}
                            case {4,6} % NV(chosen)-NV(unchosen)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': NVch-NVunch'];
                                dispRegFn([num2str(n_regs),') choice: NVch-NVunch ',run_nb_nm],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 5 % p(chosen)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': p(chosen)'];
                                dispRegFn([num2str(n_regs),') choice: p(chosen) ',run_nb_nm],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % Reward > Punishment
                        if GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).R_vs_P == 1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': R-P'];
                            dispRegFn([num2str(n_regs),') choice: Reward>Punishment ',run_nb_nm],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        end
                        
                        % binary variable indicating when choice = high effort option
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).choiceHighE
                            case 0
                            case {1,2}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': choice = highE'];
                                dispRegFn([num2str(n_regs),') choice: choice hE ',run_nb_nm],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % reward variable option
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).R_varOption
                            case 0
                            case {1,3}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': R amount highE'];
                                dispRegFn([num2str(n_regs),') choice: R amount hE ',run_nb_nm],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case {2,4}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': R level highE'];
                                dispRegFn([num2str(n_regs),') choice: R level hE ',run_nb_nm],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % reward chosen option
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).R_chosen
                            case 0
                            case {1,3}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': R amount chosen'];
                                dispRegFn([num2str(n_regs),') choice: R amount chosen ',run_nb_nm],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case {2,4}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': R level chosen'];
                                dispRegFn([num2str(n_regs),') choice: R level chosen ',run_nb_nm],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % punishment variable option
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).P_varOption
                            case 0
                            case {1,3}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': P amount highE'];
                                dispRegFn([num2str(n_regs),') choice: P amount hE ',run_nb_nm],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case {2,4}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': P level highE'];
                                dispRegFn([num2str(n_regs),') choice: P level hE ',run_nb_nm],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % punishment chosen option
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).P_chosen
                            case 0
                            case {1,3}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': P amount chosen'];
                                dispRegFn([num2str(n_regs),') choice: P amount chosen ',run_nb_nm],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case {2,4}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': P level chosen'];
                                dispRegFn([num2str(n_regs),') choice: P level chosen ',run_nb_nm],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % money left
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_left
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': money left'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': money left (amounts)'],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': money left'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': |money left| (amounts)'],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': money left'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': money left (levels)'],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': money left'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': |money left| (levels)'],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % money right
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_right
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': money right'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': money right (amounts)'],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': money right'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': |money right| (amounts)'],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': money right'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': money right (levels)'],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': money right'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': |money right| (levels)'],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % money chosen option
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_chosen
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': money chosen'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': money chosen (amounts)'],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': money chosen'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': |money chosen| (amounts)'],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': money chosen'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': money chosen (levels)'],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': money chosen'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': |money chosen| (levels)'],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % money unchosen option
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_unchosen
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': money unchosen'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': money unchosen (amounts)'],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': money unchosen'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': |money unchosen| (amounts)'],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': money unchosen'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': money unchosen (levels)'],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': money unchosen'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': |money unchosen| (levels)'],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % money associated to the option which varies (the non-default
                        % option)
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_varOption
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': money hE'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': money hE (amount) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': money hE'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': |money| hE (amount) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': money hE'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': money hE (level) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': money hE'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': |money| hE (level) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % money chosen - money unchosen
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_ch_min_unch
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': money ch-unch'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': money ch-unch '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % money chosen - money default
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_ch_min_fixOption
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': money ch-def'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': money chosen-default (amount) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': money ch-def'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': money chosen-default (level) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % money sum of both options
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_sum
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': money sum'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': money sum '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % effort left
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).E_left
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': effort left'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': effort left (level) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': effort left'];
                                switch task_id_nm
                                    case 'Ep'
                                        dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': effort left (durations) '],dispRegs);
                                    case 'Em'
                                        dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': effort left (nb answers to give) '],dispRegs);
                                end
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % effort right
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).E_right
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': effort right'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': effort right (level) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': effort right'];
                                switch task_id_nm
                                    case 'Ep'
                                        dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': effort right (durations) '],dispRegs);
                                    case 'Em'
                                        dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': effort right (nb answers to give) '],dispRegs);
                                end
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % effort chosen option
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).E_chosen
                            case 0
                            case {1,4}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': effort chosen'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': effort chosen (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': effort chosen'];
                                switch task_id_nm
                                    case 'Ep'
                                        dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': effort chosen (durations) '],dispRegs);
                                    case 'Em'
                                        dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': effort chosen (nb answers to give) '],dispRegs);
                                end
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': effort-by-choice'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': effort-by-choice '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % effort unchosen option
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).E_unchosen
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': effort unchosen'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': effort unchosen (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': effort unchosen'];
                                switch task_id_nm
                                    case 'Ep'
                                        dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': effort unchosen (durations) '],dispRegs);
                                    case 'Em'
                                        dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': effort unchosen (nb answers to give) '],dispRegs);
                                end
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % effort level hE
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).E_varOption
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': effort level hE'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': effort level hE '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': effort amount hE'];
                                switch task_id_nm
                                    case 'Ep'
                                        dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': effort amount hE (durations) '],dispRegs);
                                    case 'Em'
                                        dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': effort amount hE (nb answers to give) '],dispRegs);
                                end
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % effort chosen - unchosen option
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).E_ch_min_unch
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': effort ch-unch'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': effort ch-unch (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': effort ch-unch'];
                                switch task_id_nm
                                    case 'Ep'
                                        dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': effort ch-unch (durations) '],dispRegs);
                                    case 'Em'
                                        dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': effort ch-unch (nb answers to give) '],dispRegs);
                                end
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % effort chosen - fixed low effort option
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).E_ch_min_fixOption
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': effort ch-def'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': effort ch-def (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % effort sum of both options
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).E_sum
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': effort sum'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': effort sum (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': effort sum'];
                                switch task_id_nm
                                    case 'Ep'
                                        dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': effort sum (durations) '],dispRegs);
                                    case 'Em'
                                        dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': effort sum (nb answers to give) '],dispRegs);
                                end
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % (money)*(effort) for high effort option
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_level_x_E_varOption
                            case 0
                            case 1 % (money level)*(effort level) high E option
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': money x effort hE'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': money x effort hE '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % (money)*(effort) for chosen option
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_level_x_E_chosen
                            case 0
                            case 1 % (money level)*(effort level) chosen option
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': money x effort chosen'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': money x effort chosen '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % (reward)*(effort) for high effort option
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).R_level_x_E_varOption
                            case 0
                            case 1 % (R level)*(effort level) high E option
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': R x E hE'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': R x E hE '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % (reward)*(effort) for chosen option
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).R_level_x_E_chosen
                            case 0
                            case 1 % (R level)*(effort level) chosen option
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': R x E chosen'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': R x E chosen '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % (punishment)*(effort) for high effort option
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).P_level_x_E_varOption
                            case 0
                            case 1 % (P level)*(effort level) high E option
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': P x E hE'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': P x E hE '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % (punishment)*(effort) for chosen option
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).P_level_x_E_chosen
                            case 0
                            case 1 % (P level)*(effort level) chosen option
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': P x E chosen'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': P x E chosen '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % net value chosen option
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).NV_chosen
                            case {0,4,5,6}
                            case {1,3} % NV(chosen)-NV(unchosen)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': NVch-NVunch'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': NVch-NVunch '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2 % p(chosen)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': p(chosen)'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': p(chosen) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % net value variable option
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).NV_varOption
                            case 0
                            case 1 % NV(high E) - NV(low E)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': net value high E - low E'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': net value high E - low E '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2 % |NV(high E) - NV(low E)|
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': abs(NVhE-NVlE)'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': |net value high E - low E| '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3 % p(choice = high E)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': p(choice = high E)'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': p(choice = high E) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4 % NV(high E) - NV(low E) + bias
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': net value high E - low E + bias'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': net value high E - low E + bias '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 5 % |NV(high E) - NV(low E) + bias|
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': abs(NVhE-NVlE+bias)'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': |net value high E - low E + bias| '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % net value variable option bis
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).NV_varOption_bis
                            case 0
                            case 1 % NV(high E) - NV(low E)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': net value high E - low E'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': net value high E - low E '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2 % |NV(high E) - NV(low E)|
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': abs(NVhE-NVlE)'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': |net value high E - low E| '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3 % p(choice = high E)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': p(choice = high E)'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': p(choice = high E) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4 % NV(high E) - NV(low E) + bias
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': net value high E - low E + bias'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': net value high E - low E + bias '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 5 % |NV(high E) - NV(low E) + bias|
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': abs(NVhE-NVlE+bias)'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': |net value high E - low E + bias| '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        if strcmp(task_id_nm, 'Ep') % physical effort only
                            % area under the curve of the force that is gonna be produced
                            switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).F_integral
                                case 0
                                case {1,3,5,7}
                                    n_regs = n_regs + 1;
                                    reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': effort integral'];
                                    dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': effort integral '],dispRegs);
                                    % if derivative added => add derivatives
                                    n_regs = n_regs + add_drv;
                                case {2,4,6,8}
                                    n_regs = n_regs + 1;
                                    reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': effort integral overshoot'];
                                    dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': effort integral overshoot '],dispRegs);
                                    % if derivative added => add derivatives
                                    n_regs = n_regs + add_drv;
                                otherwise
                                    error('not ready yet');
                            end
                            
                            % fatigue
                            switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).fatigue
                                case 0
                                case {1,2}
                                    n_regs = n_regs + 1;
                                    reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': fatigue'];
                                    dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': fatigue '],dispRegs);
                                    % if derivative added => add derivatives
                                    n_regs = n_regs + add_drv;
                                otherwise
                                    error('not ready yet');
                            end
                            
                            % Ech*fatigue
                            switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).Ech_x_fatigue
                                case 0
                                case 1
                                    n_regs = n_regs + 1;
                                    reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': Ech_x_fatigue'];
                                    dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': Ech_x_fatigue '],dispRegs);
                                    % if derivative added => add derivatives
                                    n_regs = n_regs + add_drv;
                                otherwise
                                    error('not ready yet');
                            end
                        end % physical effort filter
                        
                        if strcmp(task_id_nm,'Em') % mental effort only
                            % efficacy of the next trial
                            switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).efficacy
                                case 0
                                case {1,2,3,4,5,6,7,8}
                                    n_regs = n_regs + 1;
                                    reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': efficacy'];
                                    dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': efficacy '],dispRegs);
                                    % if derivative added => add derivatives
                                    n_regs = n_regs + add_drv;
                                otherwise
                                    error('not ready yet');
                            end
                            
                            % efficacy during the previous trial
                            switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).prevEfficacy
                                case 0
                                case {1,2,3,4,5,6,7,8}
                                    n_regs = n_regs + 1;
                                    reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': previous efficacy'];
                                    dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': previous efficacy '],dispRegs);
                                    % if derivative added => add derivatives
                                    n_regs = n_regs + add_drv;
                                otherwise
                                    error('not ready yet');
                            end
                            
                            % (effort chosen)*(efficacy during the previous trial)
                            switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).Ech_x_prevEfficacy
                                case 0
                                case {1,2,3,4}
                                    n_regs = n_regs + 1;
                                    reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': (Ech)x(previous efficacy)'];
                                    dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': (Ech)x(previous efficacy) '],dispRegs);
                                    % if derivative added => add derivatives
                                    n_regs = n_regs + add_drv;
                                otherwise
                                    error('not ready yet');
                            end
                        end % mental effort filter
                        
                        % trial number
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).trialN
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': trial number'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': trial number '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': (trial number)x(Echosen-Eunchosen) (E levels)'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': (trial number)x(effort chosen-effort unchosen) (effort levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': trial number x (EnonDef-Edef) (E levels)'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': (trial number)x(effort non-default - effort default) (effort levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': (trial number)x(EnonDef) (E levels)'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': (trial number)x(effort non-default) (effort levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % confidence
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).confidence
                            case 0
                            case 1 % confidence ratings 0/1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': confidence'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': confidence (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case {2,3,4} % confidence inferred by the model
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': confidence'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': confidence (inferred by the model) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % RT (last regressor)
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).RT
                            case {0,4,5,6}
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': RT'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': RT (raw) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': RT'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': RT (zscored per run) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',run_nb_nm,': RT'];
                                dispRegFn([num2str(n_regs),') choice ',run_nb_nm,': RT (zscored per subject ie across all runs) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                    end % run loop
                case {2,4} % all tasks independent but sessions pooled
                    for iTask = 1:nTasks
                        task_nm = tasks{iTask};
                        
                        %% choice onset
                        n_regs = n_regs + 1;
                        reg_names{n_regs} = ['ONSET choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm];
                        dispRegFn([num2str(n_regs),') ONSET choice display options ',...
                            RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': ',GLMprm.model_onset.(task_id_nm).choice,' '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs = n_regs + add_drv;
                        
                        %% choice regressors
                        
                        % RT (first regressor)
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).RT
                            case {0,1,2,3}
                            case 4
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',...
                                    RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': RT'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': RT (raw) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 5
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',...
                                    RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': RT'];
                                dispRegFn([num2str(n_regs),') choice: RT (zscored per run) ',task_nm],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 6
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',...
                                    RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': RT'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,' RT (zscored per subject ie across all runs) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % net value chosen option (first regressor)
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).NV_chosen
                            case {0,1,2,3}
                            case {4,6} % NV(chosen)-NV(unchosen)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': NVch-NVunch'];
                                dispRegFn([num2str(n_regs),') choice: NVch-NVunch ',task_nm],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 5 % p(chosen)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': p(chosen)'];
                                dispRegFn([num2str(n_regs),') choice: p(chosen) ',task_nm],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % Reward > Punishment
                        if GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).R_vs_P == 1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': R-P'];
                            dispRegFn([num2str(n_regs),') choice: Reward>Punishment ',task_nm],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        end
                        
                        % binary variable indicating when choice = high effort option
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).choiceHighE
                            case 0
                            case {1,2}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': choice = highE'];
                                dispRegFn([num2str(n_regs),') choice: choice hE ',task_nm],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % reward variable option
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).R_varOption
                            case 0
                            case {1,3}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': R amount highE'];
                                dispRegFn([num2str(n_regs),') choice: R amount hE ',task_nm],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case {2,4}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': R level highE'];
                                dispRegFn([num2str(n_regs),') choice: R level hE ',task_nm],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % reward chosen option
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).R_chosen
                            case 0
                            case {1,3}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': R amount chosen'];
                                dispRegFn([num2str(n_regs),') choice: R amount chosen ',task_nm],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case {2,4}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': R level chosen'];
                                dispRegFn([num2str(n_regs),') choice: R level chosen ',task_nm],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % punishment variable option
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).P_varOption
                            case 0
                            case {1,3}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': P amount highE'];
                                dispRegFn([num2str(n_regs),') choice: P amount hE ',task_nm],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case {2,4}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': P level highE'];
                                dispRegFn([num2str(n_regs),') choice: P level hE ',task_nm],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % punishment chosen option
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).P_chosen
                            case 0
                            case {1,3}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': P amount chosen'];
                                dispRegFn([num2str(n_regs),') choice: P amount chosen ',task_nm],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case {2,4}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': P level chosen'];
                                dispRegFn([num2str(n_regs),') choice: P level chosen ',task_nm],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % money left
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_left
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': money left'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': money left (amounts)'],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': money left'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': |money left| (amounts)'],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': money left'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': money left (levels)'],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': money left'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': |money left| (levels)'],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % money right
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_right
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': money right'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': money right (amounts)'],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': money right'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': |money right| (amounts)'],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': money right'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': money right (levels)'],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': money right'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': |money right| (levels)'],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % money chosen option
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_chosen
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': money chosen'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': money chosen (amounts)'],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': money chosen'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': |money chosen| (amounts)'],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': money chosen'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': money chosen (levels)'],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': money chosen'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': |money chosen| (levels)'],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % money unchosen option
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_unchosen
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': money unchosen'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': money unchosen (amounts)'],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': money unchosen'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': |money unchosen| (amounts)'],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': money unchosen'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': money unchosen (levels)'],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': money unchosen'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': |money unchosen| (levels)'],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % money associated to the option which varies (the non-default
                        % option)
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_varOption
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': money hE'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': money hE (amount) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': money hE'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': |money| hE (amount) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': money hE'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': money hE (level) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': money hE'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': |money| hE (level) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % money chosen - money unchosen
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_ch_min_unch
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': money ch-unch'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': money ch-unch '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % money chosen - money default
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_ch_min_fixOption
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': money ch-def'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': money chosen-default (amount) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': money ch-def'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': money chosen-default (level) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % money sum of both options
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_sum
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': money sum'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': money sum '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % effort left
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).E_left
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': effort left'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': effort left (level) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': effort left'];
                                switch task_id_nm
                                    case 'Ep'
                                        dispRegFn([num2str(n_regs),') choice ',task_nm,': effort left (durations) '],dispRegs);
                                    case 'Em'
                                        dispRegFn([num2str(n_regs),') choice ',task_nm,': effort left (nb answers to give) '],dispRegs);
                                end
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % effort right
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).E_right
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': effort right'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': effort right (level) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': effort right'];
                                switch task_id_nm
                                    case 'Ep'
                                        dispRegFn([num2str(n_regs),') choice ',task_nm,': effort right (durations) '],dispRegs);
                                    case 'Em'
                                        dispRegFn([num2str(n_regs),') choice ',task_nm,': effort right (nb answers to give) '],dispRegs);
                                end
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % effort chosen option
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).E_chosen
                            case 0
                            case {1,4}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': effort chosen'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': effort chosen (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': effort chosen'];
                                switch task_id_nm
                                    case 'Ep'
                                        dispRegFn([num2str(n_regs),') choice ',task_nm,': effort chosen (durations) '],dispRegs);
                                    case 'Em'
                                        dispRegFn([num2str(n_regs),') choice ',task_nm,': effort chosen (nb answers to give) '],dispRegs);
                                end
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': effort-by-choice'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': effort-by-choice '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % effort unchosen option
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).E_unchosen
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': effort unchosen'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': effort unchosen (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': effort unchosen'];
                                switch task_id_nm
                                    case 'Ep'
                                        dispRegFn([num2str(n_regs),') choice ',task_nm,': effort unchosen (durations) '],dispRegs);
                                    case 'Em'
                                        dispRegFn([num2str(n_regs),') choice ',task_nm,': effort unchosen (nb answers to give) '],dispRegs);
                                end
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % effort level hE
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).E_varOption
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': effort level hE'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': effort level hE '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': effort amount hE'];
                                switch task_id_nm
                                    case 'Ep'
                                        dispRegFn([num2str(n_regs),') choice ',task_nm,': effort amount hE (durations) '],dispRegs);
                                    case 'Em'
                                        dispRegFn([num2str(n_regs),') choice ',task_nm,': effort amount hE (nb answers to give) '],dispRegs);
                                end
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % effort chosen - unchosen option
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).E_ch_min_unch
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': effort ch-unch'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': effort ch-unch (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': effort ch-unch'];
                                switch task_id_nm
                                    case 'Ep'
                                        dispRegFn([num2str(n_regs),') choice ',task_nm,': effort ch-unch (durations) '],dispRegs);
                                    case 'Em'
                                        dispRegFn([num2str(n_regs),') choice ',task_nm,': effort ch-unch (nb answers to give) '],dispRegs);
                                end
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % effort chosen - fixed low effort option
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).E_ch_min_fixOption
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': effort ch-def'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': effort ch-def (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % effort sum of both options
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).E_sum
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': effort sum'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': effort sum (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': effort sum'];
                                switch task_id_nm
                                    case 'Ep'
                                        dispRegFn([num2str(n_regs),') choice ',task_nm,': effort sum (durations) '],dispRegs);
                                    case 'Em'
                                        dispRegFn([num2str(n_regs),') choice ',task_nm,': effort sum (nb answers to give) '],dispRegs);
                                end
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % (money)*(effort) for high effort option
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_level_x_E_varOption
                            case 0
                            case 1 % (money level)*(effort level) high E option
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': money x effort hE'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': money x effort hE '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % (money)*(effort) for chosen option
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_level_x_E_chosen
                            case 0
                            case 1 % (money level)*(effort level) chosen option
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': money x effort chosen'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': money x effort chosen '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % (reward)*(effort) for high effort option
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).R_level_x_E_varOption
                            case 0
                            case 1 % (R level)*(effort level) high E option
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': R x E hE'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': R x E hE '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % (reward)*(effort) for chosen option
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).R_level_x_E_chosen
                            case 0
                            case 1 % (R level)*(effort level) chosen option
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': R x E chosen'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': R x E chosen '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % (punishment)*(effort) for high effort option
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).P_level_x_E_varOption
                            case 0
                            case 1 % (P level)*(effort level) high E option
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': P x E hE'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': P x E hE '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % (punishment)*(effort) for chosen option
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).P_level_x_E_chosen
                            case 0
                            case 1 % (P level)*(effort level) chosen option
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': P x E chosen'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': P x E chosen '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % net value chosen option
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).NV_chosen
                            case {0,4,5,6}
                            case {1,3} % NV(chosen)-NV(unchosen)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': NVch-NVunch'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': NVch-NVunch '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2 % p(chosen)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': p(chosen)'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': p(chosen) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % net value variable option
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).NV_varOption
                            case 0
                            case 1 % NV(high E) - NV(low E)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': net value high E - low E'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': net value high E - low E '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2 % |NV(high E) - NV(low E)|
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': abs(NVhE-NVlE)'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': |net value high E - low E| '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3 % p(choice = high E)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': p(choice = high E)'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': p(choice = high E) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4 % NV(high E) - NV(low E) + bias
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': net value high E - low E + bias'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': net value high E - low E + bias '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 5 % |NV(high E) - NV(low E) + bias|
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': abs(NVhE-NVlE+bias)'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': |net value high E - low E + bias| '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % net value variable option bis
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).NV_varOption_bis
                            case 0
                            case 1 % NV(high E) - NV(low E)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': net value high E - low E'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': net value high E - low E '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2 % |NV(high E) - NV(low E)|
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': abs(NVhE-NVlE)'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': |net value high E - low E| '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3 % p(choice = high E)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': p(choice = high E)'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': p(choice = high E) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4 % NV(high E) - NV(low E) + bias
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': net value high E - low E + bias'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': net value high E - low E + bias '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 5 % |NV(high E) - NV(low E) + bias|
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': abs(NVhE-NVlE+bias)'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': |net value high E - low E + bias| '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        if strcmp(task_id_nm, 'Ep') % physical effort only
                            % area under the curve of the force that is gonna be produced
                            switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).F_integral
                                case 0
                                case {1,3,5,7}
                                    n_regs = n_regs + 1;
                                    reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': effort integral'];
                                    dispRegFn([num2str(n_regs),') choice ',task_nm,': effort integral '],dispRegs);
                                    % if derivative added => add derivatives
                                    n_regs = n_regs + add_drv;
                                case {2,4,6,8}
                                    n_regs = n_regs + 1;
                                    reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': effort integral overshoot'];
                                    dispRegFn([num2str(n_regs),') choice ',task_nm,': effort integral overshoot '],dispRegs);
                                    % if derivative added => add derivatives
                                    n_regs = n_regs + add_drv;
                                otherwise
                                    error('not ready yet');
                            end
                            
                            % fatigue
                            switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).fatigue
                                case 0
                                case {1,2}
                                    n_regs = n_regs + 1;
                                    reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': fatigue'];
                                    dispRegFn([num2str(n_regs),') choice ',task_nm,': fatigue '],dispRegs);
                                    % if derivative added => add derivatives
                                    n_regs = n_regs + add_drv;
                                otherwise
                                    error('not ready yet');
                            end
                            
                            % Ech*fatigue
                            switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).Ech_x_fatigue
                                case 0
                                case 1
                                    n_regs = n_regs + 1;
                                    reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': Ech_x_fatigue'];
                                    dispRegFn([num2str(n_regs),') choice ',task_nm,': Ech_x_fatigue '],dispRegs);
                                    % if derivative added => add derivatives
                                    n_regs = n_regs + add_drv;
                                otherwise
                                    error('not ready yet');
                            end
                        end % physical effort filter
                        
                        if strcmp(task_id_nm,'Em') % mental effort only
                            % efficacy of the next trial
                            switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).efficacy
                                case 0
                                case {1,2,3,4,5,6,7,8}
                                    n_regs = n_regs + 1;
                                    reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': efficacy'];
                                    dispRegFn([num2str(n_regs),') choice ',task_nm,': efficacy '],dispRegs);
                                    % if derivative added => add derivatives
                                    n_regs = n_regs + add_drv;
                                otherwise
                                    error('not ready yet');
                            end
                            
                            % efficacy during the previous trial
                            switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).prevEfficacy
                                case 0
                                case {1,2,3,4,5,6,7,8}
                                    n_regs = n_regs + 1;
                                    reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': previous efficacy'];
                                    dispRegFn([num2str(n_regs),') choice ',task_nm,': previous efficacy '],dispRegs);
                                    % if derivative added => add derivatives
                                    n_regs = n_regs + add_drv;
                                otherwise
                                    error('not ready yet');
                            end
                            
                            % (effort chosen)*(efficacy during the previous trial)
                            switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).Ech_x_prevEfficacy
                                case 0
                                case {1,2,3,4}
                                    n_regs = n_regs + 1;
                                    reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': (Ech)x(previous efficacy)'];
                                    dispRegFn([num2str(n_regs),') choice ',task_nm,': (Ech)x(previous efficacy) '],dispRegs);
                                    % if derivative added => add derivatives
                                    n_regs = n_regs + add_drv;
                                otherwise
                                    error('not ready yet');
                            end
                        end % mental effort filter
                        
                        % trial number
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).trialN
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': trial number'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': trial number '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': (trial number)x(Echosen-Eunchosen) (E levels)'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': (trial number)x(effort chosen-effort unchosen) (effort levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': trial number x (EnonDef-Edef) (E levels)'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': (trial number)x(effort non-default - effort default) (effort levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': (trial number)x(EnonDef) (E levels)'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': (trial number)x(effort non-default) (effort levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % confidence
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).confidence
                            case 0
                            case 1 % confidence ratings 0/1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': confidence'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': confidence (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case {2,3,4} % confidence inferred by the model
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': confidence'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': confidence (inferred by the model) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % RT (last regressor)
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).RT
                            case {0,4,5,6}
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': RT'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': RT (raw) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': RT'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': RT (zscored per run) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,' ',task_nm,': RT'];
                                dispRegFn([num2str(n_regs),') choice ',task_nm,': RT (zscored per subject ie across all runs) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                    end % task loop
                case {3,5} % all pooled across sessions
                    
                    %% choice onset
                    n_regs = n_regs + 1;
                    reg_names{n_regs} = ['ONSET choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm];
                    dispRegFn([num2str(n_regs),') ONSET choice display options ',...
                        RP_dispChoice_nm,' ',splitE_dispChoice_nm,': ',GLMprm.model_onset.(task_id_nm).choice,' '],dispRegs);
                    % if derivative added => add derivatives
                    n_regs = n_regs + add_drv;
                    
                    %% choice regressors
                    
                    % RT (first regressor)
                    switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).RT
                        case {0,1,2,3}
                        case 4
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',...
                                RP_dispChoice_nm,' ',splitE_dispChoice_nm,': RT'];
                            dispRegFn([num2str(n_regs),') choice: RT (raw) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 5
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',...
                                RP_dispChoice_nm,' ',splitE_dispChoice_nm,': RT'];
                            dispRegFn([num2str(n_regs),') choice: RT (zscored per run) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 6
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',...
                                RP_dispChoice_nm,' ',splitE_dispChoice_nm,': RT'];
                            dispRegFn([num2str(n_regs),') choice: RT (zscored per subject ie across all runs) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % net value chosen option (first regressor)
                    switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).NV_chosen
                        case {0,1,2,3}
                        case {4,6} % NV(chosen)-NV(unchosen)
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': NVch-NVunch'];
                            dispRegFn([num2str(n_regs),') choice: NVch-NVunch '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 5 % p(chosen)
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': p(chosen)'];
                            dispRegFn([num2str(n_regs),') choice: p(chosen) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % Reward > Punishment
                    if GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).R_vs_P == 1
                        n_regs = n_regs + 1;
                        reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': R-P'];
                        dispRegFn([num2str(n_regs),') choice: Reward>Punishment '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs = n_regs + add_drv;
                    end
                    
                    % binary variable indicating when choice = high effort option
                    switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).choiceHighE
                        case 0
                        case {1,2}
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': choice = highE'];
                            dispRegFn([num2str(n_regs),') choice: choice hE '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % reward variable option
                    switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).R_varOption
                        case 0
                        case {1,3}
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': R amount highE'];
                            dispRegFn([num2str(n_regs),') choice: R amount hE '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case {2,4}
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': R level highE'];
                            dispRegFn([num2str(n_regs),') choice: R level hE '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % reward chosen option
                    switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).R_chosen
                        case 0
                        case {1,3}
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': R amount chosen'];
                            dispRegFn([num2str(n_regs),') choice: R amount chosen '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case {2,4}
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': R level chosen'];
                            dispRegFn([num2str(n_regs),') choice: R level chosen '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % punishment variable option
                    switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).P_varOption
                        case 0
                        case {1,3}
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': P amount highE'];
                            dispRegFn([num2str(n_regs),') choice: P amount hE '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case {2,4}
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': P level highE'];
                            dispRegFn([num2str(n_regs),') choice: P level hE '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % punishment chosen option
                    switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).P_chosen
                        case 0
                        case {1,3}
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': P amount chosen'];
                            dispRegFn([num2str(n_regs),') choice: P amount chosen '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case {2,4}
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': P level chosen'];
                            dispRegFn([num2str(n_regs),') choice: P level chosen '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % money left
                    switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_left
                        case 0
                        case 1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money left'];
                            dispRegFn([num2str(n_regs),') choice: money left (amounts)'],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money left'];
                            dispRegFn([num2str(n_regs),') choice: |money left| (amounts)'],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 3
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money left'];
                            dispRegFn([num2str(n_regs),') choice: money left (levels)'],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 4
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money left'];
                            dispRegFn([num2str(n_regs),') choice: |money left| (levels)'],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % money right
                    switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_right
                        case 0
                        case 1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money right'];
                            dispRegFn([num2str(n_regs),') choice: money right (amounts)'],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money right'];
                            dispRegFn([num2str(n_regs),') choice: |money right| (amounts)'],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 3
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money right'];
                            dispRegFn([num2str(n_regs),') choice: money right (levels)'],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 4
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money right'];
                            dispRegFn([num2str(n_regs),') choice: |money right| (levels)'],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % money chosen option
                    switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_chosen
                        case 0
                        case 1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money chosen'];
                            dispRegFn([num2str(n_regs),') choice: money chosen (amounts)'],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money chosen'];
                            dispRegFn([num2str(n_regs),') choice: |money chosen| (amounts)'],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 3
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money chosen'];
                            dispRegFn([num2str(n_regs),') choice: money chosen (levels)'],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 4
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money chosen'];
                            dispRegFn([num2str(n_regs),') choice: |money chosen| (levels)'],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % money unchosen option
                    switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_unchosen
                        case 0
                        case 1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money unchosen'];
                            dispRegFn([num2str(n_regs),') choice: money unchosen (amounts)'],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money unchosen'];
                            dispRegFn([num2str(n_regs),') choice: |money unchosen| (amounts)'],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 3
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money unchosen'];
                            dispRegFn([num2str(n_regs),') choice: money unchosen (levels)'],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 4
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money unchosen'];
                            dispRegFn([num2str(n_regs),') choice: |money unchosen| (levels)'],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % money associated to the option which varies (the non-default
                    % option)
                    switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_varOption
                        case 0
                        case 1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money hE'];
                            dispRegFn([num2str(n_regs),') choice: money hE (amount) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money hE'];
                            dispRegFn([num2str(n_regs),') choice: |money| hE (amount) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 3
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money hE'];
                            dispRegFn([num2str(n_regs),') choice: money hE (level) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 4
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money hE'];
                            dispRegFn([num2str(n_regs),') choice: |money| hE (level) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % money chosen - money unchosen
                    switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_ch_min_unch
                        case 0
                        case 1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money ch-unch'];
                            dispRegFn([num2str(n_regs),') choice: money ch-unch '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % money chosen - money default
                    switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_ch_min_fixOption
                        case 0
                        case 1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money ch-def'];
                            dispRegFn([num2str(n_regs),') choice: money chosen-default (amount) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money ch-def'];
                            dispRegFn([num2str(n_regs),') choice: money chosen-default (level) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % money sum of both options
                    switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_sum
                        case 0
                        case 1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money sum'];
                            dispRegFn([num2str(n_regs),') choice: money sum '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % effort left
                    switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).E_left
                        case 0
                        case 1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': effort left'];
                            dispRegFn([num2str(n_regs),') choice: effort left (level) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': effort left'];
                            switch task_id_nm
                                case 'Ep'
                                    dispRegFn([num2str(n_regs),') choice: effort left (durations) '],dispRegs);
                                case 'Em'
                                    dispRegFn([num2str(n_regs),') choice: effort left (nb answers to give) '],dispRegs);
                            end
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % effort right
                    switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).E_right
                        case 0
                        case 1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': effort right'];
                            dispRegFn([num2str(n_regs),') choice: effort right (level) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': effort right'];
                            switch task_id_nm
                                case 'Ep'
                                    dispRegFn([num2str(n_regs),') choice: effort right (durations) '],dispRegs);
                                case 'Em'
                                    dispRegFn([num2str(n_regs),') choice: effort right (nb answers to give) '],dispRegs);
                            end
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % effort chosen option
                    switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).E_chosen
                        case 0
                        case {1,4}
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': effort chosen'];
                            dispRegFn([num2str(n_regs),') choice: effort chosen (levels) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': effort chosen'];
                            switch task_id_nm
                                case 'Ep'
                                    dispRegFn([num2str(n_regs),') choice: effort chosen (durations) '],dispRegs);
                                case 'Em'
                                    dispRegFn([num2str(n_regs),') choice: effort chosen (nb answers to give) '],dispRegs);
                            end
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 3
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': effort-by-choice'];
                            dispRegFn([num2str(n_regs),') choice: effort-by-choice '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % effort unchosen option
                    switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).E_unchosen
                        case 0
                        case 1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': effort unchosen'];
                            dispRegFn([num2str(n_regs),') choice: effort unchosen (levels) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': effort unchosen'];
                            switch task_id_nm
                                case 'Ep'
                                    dispRegFn([num2str(n_regs),') choice: effort unchosen (durations) '],dispRegs);
                                case 'Em'
                                    dispRegFn([num2str(n_regs),') choice: effort unchosen (nb answers to give) '],dispRegs);
                            end
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % effort level hE
                    switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).E_varOption
                        case 0
                        case 1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': effort level hE'];
                            dispRegFn([num2str(n_regs),') choice: effort level hE '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': effort amount hE'];
                            switch task_id_nm
                                case 'Ep'
                                    dispRegFn([num2str(n_regs),') choice: effort amount hE (durations) '],dispRegs);
                                case 'Em'
                                    dispRegFn([num2str(n_regs),') choice: effort amount hE (nb answers to give) '],dispRegs);
                            end
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % effort chosen - unchosen option
                    switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).E_ch_min_unch
                        case 0
                        case 1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': effort ch-unch'];
                            dispRegFn([num2str(n_regs),') choice: effort ch-unch (levels) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': effort ch-unch'];
                            switch task_id_nm
                                case 'Ep'
                                    dispRegFn([num2str(n_regs),') choice: effort ch-unch (durations) '],dispRegs);
                                case 'Em'
                                    dispRegFn([num2str(n_regs),') choice: effort ch-unch (nb answers to give) '],dispRegs);
                            end
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % effort chosen - fixed low effort option
                    switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).E_ch_min_fixOption
                        case 0
                        case 1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': effort ch-def'];
                            dispRegFn([num2str(n_regs),') choice: effort ch-def (levels) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % effort sum of both options
                    switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).E_sum
                        case 0
                        case 1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': effort sum'];
                            dispRegFn([num2str(n_regs),') choice: effort sum (levels) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': effort sum'];
                            switch task_id_nm
                                case 'Ep'
                                    dispRegFn([num2str(n_regs),') choice: effort sum (durations) '],dispRegs);
                                case 'Em'
                                    dispRegFn([num2str(n_regs),') choice: effort sum (nb answers to give) '],dispRegs);
                            end
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % (money)*(effort) for high effort option
                    switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_level_x_E_varOption
                        case 0
                        case 1 % (money level)*(effort level) high E option
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money x effort hE'];
                            dispRegFn([num2str(n_regs),') choice: money x effort hE '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % (money)*(effort) for chosen option
                    switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_level_x_E_chosen
                        case 0
                        case 1 % (money level)*(effort level) chosen option
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money x effort chosen'];
                            dispRegFn([num2str(n_regs),') choice: money x effort chosen '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % (reward)*(effort) for high effort option
                    switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).R_level_x_E_varOption
                        case 0
                        case 1 % (R level)*(effort level) high E option
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': R x E hE'];
                            dispRegFn([num2str(n_regs),') choice: R x E hE '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % (reward)*(effort) for chosen option
                    switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).R_level_x_E_chosen
                        case 0
                        case 1 % (R level)*(effort level) chosen option
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': R x E chosen'];
                            dispRegFn([num2str(n_regs),') choice: R x E chosen '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % (punishment)*(effort) for high effort option
                    switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).P_level_x_E_varOption
                        case 0
                        case 1 % (P level)*(effort level) high E option
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': P x E hE'];
                            dispRegFn([num2str(n_regs),') choice: P x E hE '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % (punishment)*(effort) for chosen option
                    switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).P_level_x_E_chosen
                        case 0
                        case 1 % (P level)*(effort level) chosen option
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': P x E chosen'];
                            dispRegFn([num2str(n_regs),') choice: P x E chosen '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % net value chosen option
                    switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).NV_chosen
                        case {0,4,5,6}
                        case {1,3} % NV(chosen)-NV(unchosen)
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': NVch-NVunch'];
                            dispRegFn([num2str(n_regs),') choice: NVch-NVunch '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2 % p(chosen)
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': p(chosen)'];
                            dispRegFn([num2str(n_regs),') choice: p(chosen) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % net value variable option
                    switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).NV_varOption
                        case 0
                        case 1 % NV(high E) - NV(low E)
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': net value high E - low E'];
                            dispRegFn([num2str(n_regs),') choice: net value high E - low E '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2 % |NV(high E) - NV(low E)|
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': abs(NVhE-NVlE)'];
                            dispRegFn([num2str(n_regs),') choice: |net value high E - low E| '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 3 % p(choice = high E)
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': p(choice = high E)'];
                            dispRegFn([num2str(n_regs),') choice: p(choice = high E) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 4 % NV(high E) - NV(low E) + bias
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': net value high E - low E + bias'];
                            dispRegFn([num2str(n_regs),') choice: net value high E - low E + bias '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 5 % |NV(high E) - NV(low E) + bias|
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': abs(NVhE-NVlE+bias)'];
                            dispRegFn([num2str(n_regs),') choice: |net value high E - low E + bias| '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % net value variable option bis
                    switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).NV_varOption_bis
                        case 0
                        case 1 % NV(high E) - NV(low E)
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': net value high E - low E'];
                            dispRegFn([num2str(n_regs),') choice: net value high E - low E '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2 % |NV(high E) - NV(low E)|
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': abs(NVhE-NVlE)'];
                            dispRegFn([num2str(n_regs),') choice: |net value high E - low E| '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 3 % p(choice = high E)
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': p(choice = high E)'];
                            dispRegFn([num2str(n_regs),') choice: p(choice = high E) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 4 % NV(high E) - NV(low E) + bias
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': net value high E - low E + bias'];
                            dispRegFn([num2str(n_regs),') choice: net value high E - low E + bias '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 5 % |NV(high E) - NV(low E) + bias|
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': abs(NVhE-NVlE+bias)'];
                            dispRegFn([num2str(n_regs),') choice: |net value high E - low E + bias| '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    if strcmp(task_id_nm, 'Ep') % physical effort only
                        % area under the curve of the force that is gonna be produced
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).F_integral
                            case 0
                            case {1,3,5,7}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': effort integral'];
                                dispRegFn([num2str(n_regs),') choice: effort integral '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case {2,4,6,8}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': effort integral overshoot'];
                                dispRegFn([num2str(n_regs),') choice: effort integral overshoot '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % fatigue
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).fatigue
                            case 0
                            case {1,2}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': fatigue'];
                                dispRegFn([num2str(n_regs),') choice: fatigue '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % Ech*fatigue
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).Ech_x_fatigue
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': Ech_x_fatigue'];
                                dispRegFn([num2str(n_regs),') choice: Ech_x_fatigue '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                    end % physical effort filter
                    
                    if strcmp(task_id_nm,'Em') % mental effort only
                        % efficacy of the next trial
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).efficacy
                            case 0
                            case {1,2,3,4,5,6,7,8}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': efficacy'];
                                dispRegFn([num2str(n_regs),') choice: efficacy '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % efficacy during the previous trial
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).prevEfficacy
                            case 0
                            case {1,2,3,4,5,6,7,8}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': previous efficacy'];
                                dispRegFn([num2str(n_regs),') choice: previous efficacy '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % (effort chosen)*(efficacy during the previous trial)
                        switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).Ech_x_prevEfficacy
                            case 0
                            case {1,2,3,4}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': (Ech)x(previous efficacy)'];
                                dispRegFn([num2str(n_regs),') choice: (Ech)x(previous efficacy) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                    end % mental effort filter
                    
                    % trial number
                    switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).trialN
                        case 0
                        case 1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': trial number'];
                            dispRegFn([num2str(n_regs),') choice: trial number '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': (trial number)x(Echosen-Eunchosen) (E levels)'];
                            dispRegFn([num2str(n_regs),') choice: (trial number)x(effort chosen-effort unchosen) (effort levels) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 3
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': trial number x (EnonDef-Edef) (E levels)'];
                            dispRegFn([num2str(n_regs),') choice: (trial number)x(effort non-default - effort default) (effort levels) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 4
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': (trial number)x(EnonDef) (E levels)'];
                            dispRegFn([num2str(n_regs),') choice: (trial number)x(effort non-default) (effort levels) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % confidence
                    switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).confidence
                        case 0
                        case 1 % confidence ratings 0/1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': confidence'];
                            dispRegFn([num2str(n_regs),') choice: confidence (levels) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case {2,3,4} % confidence inferred by the model
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': confidence'];
                            dispRegFn([num2str(n_regs),') choice: confidence (inferred by the model) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % RT (last regressor)
                    switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).RT
                        case {0,4,5,6}
                        case 1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': RT'];
                            dispRegFn([num2str(n_regs),') choice: RT (raw) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': RT'];
                            dispRegFn([num2str(n_regs),') choice: RT (zscored per run) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 3
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': RT'];
                            dispRegFn([num2str(n_regs),') choice: RT (zscored per subject ie across all runs) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
            end % DCM mode
        end % loop effort levels
    end % loop reward/punishment
end % choice onset

%% chosen period
if ~strcmp(GLMprm.model_onset.(task_id_nm).chosen,'none')
    % check if trials are split or not
    if GLMprm.chosen.(task_id_nm).RPpool == 1 % pool reward and punishment trials
        RP_chosen = {'RP'};
    elseif GLMprm.chosen.(task_id_nm).RPpool == 0 % split reward and punishment
        RP_chosen = {'R','P'};
    end % RP pool
    n_RP_chosen = length(RP_chosen);
    
    % check if trials are split or not according to Effort levels
    switch GLMprm.chosen.(task_id_nm).splitPerE
        case 0 % pool trials
            splitE_chosen = {'E'};
        case 1 % split according to effort proposed
            splitE_chosen = {'E1','E2','E3'};
        case 2 % split according to effort chosen
            splitE_chosen = {'Ech0','Ech1','Ech2','Ech3'};
        case 3 % split according to option chosen (low/high effort)
            splitE_chosen = {'lEch','hEch'};
    end % Effort level pool
    n_splitE_chosen = length(splitE_chosen);
    
    % loop through conditions for chosen period
    for iRP_chosen = 1:n_RP_chosen
        RP_dispChosen_nm = RP_chosen{iRP_chosen};
        
        for iE_chosen = 1:n_splitE_chosen
            splitE_dispChosen_nm = splitE_chosen{iE_chosen};
            
            switch DCM_mode
                case 1 % all sessions independent
                    for iRun = 1:n_runs
                        run_nb_nm = ['r',num2str(runs.runsToKeep(iRun))];
                        
                        %% chosen onset
                        n_regs = n_regs + 1;
                        reg_names{n_regs} = ['ONSET chosen ',...
                            RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm];
                        dispRegFn([num2str(n_regs),') ONSET chosen option display ',...
                            RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': ',GLMprm.model_onset.(task_id_nm).chosen,' '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs = n_regs + add_drv;
                        
                        %% chosen regressors
                        
                        % RT (first regressor)
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).RT
                            case {0,1,2,3}
                            case 4
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': RT'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': RT (raw) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 5
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': RT'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': RT (zscored per run) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 6
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': RT'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': RT (zscored per subject ie across all runs) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % net value chosen option (first regressor)
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).NV_chosen
                            case {0,1,2,3}
                            case {4,6}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': NVch-NVunch'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': NVch-NVunch '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 5 % p(chosen)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': p(chosen)'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': p(chosen) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % reward/punishment trial
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).R_vs_P
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': R-P'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': R vs P '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % binary variable indicating when choice = high effort option
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).choiceHighE
                            case 0
                            case {1,2}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': choice = highE'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': choice hE '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % reward variable option
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).R_varOption
                            case 0
                            case {1,3}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': R amount highE'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': R amount hE '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case {2,4}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': R level highE'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': R level hE '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % reward chosen option
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).R_chosen
                            case 0
                            case {1,3}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': R amount chosen'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': R amount chosen '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case {2,4}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': R level chosen'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': R level chosen '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % punishment variable option
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).P_varOption
                            case 0
                            case {1,3}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': P amount highE'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': P amount hE '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case {2,4}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': P level highE'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': P level hE '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % punishment chosen option
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).P_chosen
                            case 0
                            case {1,3}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': P amount chosen'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': P amount chosen '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case {2,4}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': P level chosen'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': P level chosen '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % money left
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_left
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': money left'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': money left (amounts) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': money left'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': |money left| (amounts) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': money left'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': money left (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': money left'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': |money left| (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % money right
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_right
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': money right'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': money right (amounts) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': money right'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': |money right| (amounts) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': money right'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': money right (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': money right'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': |money right| (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % money chosen
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_chosen
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': money chosen'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': money chosen (amounts) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': money chosen'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': |money chosen| (amounts) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': money chosen'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': money chosen (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': money chosen'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': |money chosen| (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % money unchosen
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_unchosen
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': money unchosen'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': money unchosen (amounts) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': money unchosen'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': |money unchosen| (amounts) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': money unchosen'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': money unchosen (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': money unchosen'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': |money unchosen| (amounts) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % money associated to the option which varies (the non-default
                        % option)
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_varOption
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': money hE'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': money hE option (amount) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': money hE'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': |money| hE (amount) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': money hE'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': money hE (level) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': money hE'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': |money| hE (level) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % money chosen - money unchosen
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_ch_min_unch
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': money ch-unch'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': money ch-unch '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % money chosen - money default
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_ch_min_fixOption
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': money ch-def'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': money chosen-default (amount) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': money ch-def'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': money chosen-default (level) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % money sum of both options
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_sum
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': money sum'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': money sum '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % effort chosen
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).E_chosen
                            case 0
                            case {1,4}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': effort chosen'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': effort chosen (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                error('case not ready yet');
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': effort-by-choice'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': effort-by-choice '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % effort unchosen
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).E_unchosen
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': effort unchosen'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': effort unchosen (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('case not ready yet');
                        end
                        
                        % effort hE
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).E_varOption
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': effort level hE'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': effort level hE '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('case not ready yet');
                        end
                        
                        % effort chosen - effort unchosen
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).E_ch_min_unch
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': effort ch-unch'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': effort ch-unch (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': effort ch-unch'];
                                switch task_id_nm
                                    case 'Ep'
                                        dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': effort ch-unch (durations) '],dispRegs);
                                    case 'Em'
                                        dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': effort ch-unch (nb answers to give) '],dispRegs);
                                end
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('case not ready yet');
                        end
                        
                        % effort chosen - fixed low effort option
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).E_ch_min_fixOption
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': effort ch-def'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': effort ch-def (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('case not ready yet');
                        end
                        
                        % effort sum of both options
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).E_sum
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': effort sum'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': effort sum (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': effort sum'];
                                switch task_id_nm
                                    case 'Ep'
                                        dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': effort sum (durations) '],dispRegs);
                                    case 'Em'
                                        dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': effort sum (nb answers to give) '],dispRegs);
                                end
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('case not ready yet');
                        end
                        
                        % (money)*(effort) for high effort option
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_level_x_E_varOption
                            case 0
                            case 1 % (money level)*(effort level) high E option
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': money x effort non-default'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': money x effort hE '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % (money)*(effort) for chosen option
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_level_x_E_chosen
                            case 0
                            case 1 % (money level)*(effort level) chosen option
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': money x effort chosen'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': money x effort chosen '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % (reward)*(effort) for high effort option
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).R_level_x_E_varOption
                            case 0
                            case 1 % (R level)*(effort level) high E option
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': R x E non-default'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': R x E hE '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % (reward)*(effort) for chosen option
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).R_level_x_E_chosen
                            case 0
                            case 1 % (R level)*(effort level) chosen option
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': R x E chosen'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': R x E chosen '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % (punishment)*(effort) for high effort option
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).P_level_x_E_varOption
                            case 0
                            case 1 % (P level)*(effort level) high E option
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': P x E non-default'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': P x E hE '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % (punishment)*(effort) for chosen option
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).P_level_x_E_chosen
                            case 0
                            case 1 % (P level)*(effort level) chosen option
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': P x E chosen'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': P x E chosen '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % net value chosen option
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).NV_chosen
                            case {0,4,5,6}
                            case {1,3}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': NVch-NVunch'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': NVch-NVunch '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2 % p(chosen)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': p(chosen)'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': p(chosen) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % net value variable option
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).NV_varOption
                            case 0
                            case 1 % NV(high E) - NV(low E)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': net value high E - low E'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': net value high E - low E '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2 % |NV(high E) - NV(low E)|
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': abs(NVhE-NVlE)'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': |net value high E - low E| '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3 % p(choice = high E)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': p(choice = high E)'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': p(choice = high E) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4 % NV(high E) - NV(low E) + bias
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': net value high E - low E + bias'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': net value high E - low E + bias '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 5 % |NV(high E) - NV(low E) + bias|
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': abs(NVhE-NVlE+bias)'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': |net value high E - low E + bias| '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % net value variable option bis
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).NV_varOption_bis
                            case 0
                            case 1 % NV(high E) - NV(low E)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': net value high E - low E'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': net value high E - low E '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2 % |NV(high E) - NV(low E)|
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': abs(NVhE-NVlE)'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': |net value high E - low E| '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3 % p(choice = high E)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': p(choice = high E)'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': p(choice = high E) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4 % NV(high E) - NV(low E) + bias
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': net value high E - low E + bias'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': net value high E - low E + bias '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 5 % |NV(high E) - NV(low E) + bias|
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': abs(NVhE-NVlE+bias)'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': |net value high E - low E + bias| '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        if strcmp(task_id_nm, 'Ep') % physical effort only
                            % area under the curve of the force that is gonna be produced
                            switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).F_integral
                                case 0
                                case {1,3,5,7}
                                    n_regs = n_regs + 1;
                                    reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': force integral'];
                                    dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': effort integral '],dispRegs);
                                    % if derivative added => add derivatives
                                    n_regs = n_regs + add_drv;
                                case {2,4,6,8}
                                    n_regs = n_regs + 1;
                                    reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': force integral overshoot'];
                                    dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': effort integral overshoot '],dispRegs);
                                    % if derivative added => add derivatives
                                    n_regs = n_regs + add_drv;
                                otherwise
                                    error('not ready yet');
                            end
                            
                            % fatigue
                            switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).fatigue
                                case 0
                                case {1,2}
                                    n_regs = n_regs + 1;
                                    reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': fatigue'];
                                    dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': fatigue '],dispRegs);
                                    % if derivative added => add derivatives
                                    n_regs = n_regs + add_drv;
                                otherwise
                                    error('not ready yet');
                            end
                            
                            % Ech*fatigue
                            switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).Ech_x_fatigue
                                case 0
                                case 1
                                    n_regs = n_regs + 1;
                                    reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': Ech_x_fatigue'];
                                    dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': Ech_x_fatigue '],dispRegs);
                                    % if derivative added => add derivatives
                                    n_regs = n_regs + add_drv;
                                otherwise
                                    error('not ready yet');
                            end
                        end % physical effort filter
                        
                        if strcmp(task_id_nm,'Em') % mental effort only
                            % efficacy of the next trial
                            switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).efficacy
                                case 0
                                case {1,2,3,4,5,6,7,8}
                                    n_regs = n_regs + 1;
                                    reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': efficacy'];
                                    dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': efficacy '],dispRegs);
                                    % if derivative added => add derivatives
                                    n_regs = n_regs + add_drv;
                                otherwise
                                    error('not ready yet');
                            end
                            
                            % efficacy during the previous trial
                            switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).prevEfficacy
                                case 0
                                case {1,2,3,4,5,6,7,8}
                                    n_regs = n_regs + 1;
                                    reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': previous efficacy'];
                                    dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': previous efficacy '],dispRegs);
                                    % if derivative added => add derivatives
                                    n_regs = n_regs + add_drv;
                                otherwise
                                    error('not ready yet');
                            end
                            
                            % (effort chosen)*(efficacy during the previous trial)
                            switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).Ech_x_prevEfficacy
                                case 0
                                case {1,2,3,4}
                                    n_regs = n_regs + 1;
                                    reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': (Ech)x(previous efficacy)'];
                                    dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': (Ech)x(previous efficacy) '],dispRegs);
                                    % if derivative added => add derivatives
                                    n_regs = n_regs + add_drv;
                                otherwise
                                    error('not ready yet');
                            end
                        end % mental effort filter
                        
                        % trial number
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).trialN
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': trial number'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': trial number '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': (trial number)x(Echosen-Eunchosen) (E levels)'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': (trial number)x(effort chosen-effort unchosen) (effort levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': (trial number)x(EnonDef-Edef) (E levels)'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': (trial number)x(effort non-default - effort default) (effort levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': (trial number)x(EnonDef) (E levels)'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': (trial number)x(effort non-default) (effort levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % confidence
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).confidence
                            case 0
                            case 1 % confidence rating by the subjects
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': confidence'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': confidence (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case {2,3,4} % confidence inferred by the model
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': confidence'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': confidence (inferred by the model) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % RT (last regressor)
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).RT
                            case {0,4,5,6}
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': RT'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': RT (raw) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': RT'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': RT (zscored per run) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',run_nb_nm,': RT'];
                                dispRegFn([num2str(n_regs),') chosen ',run_nb_nm,': RT (zscored per subject ie across all runs) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                    end % run loop
                case 2 % all tasks independent but sessions pooled
                    for iTask = 1:nTasks
                        task_nm = tasks{iTask};
                        
                        %% chosen onset
                        n_regs = n_regs + 1;
                        reg_names{n_regs} = ['ONSET chosen ',...
                            RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm];
                        dispRegFn([num2str(n_regs),') ONSET chosen option display ',...
                            RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': ',GLMprm.model_onset.(task_id_nm).chosen,' '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs = n_regs + add_drv;
                        
                        %% chosen regressors
                        
                        % RT (first regressor)
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).RT
                            case {0,1,2,3}
                            case 4
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': RT'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': RT (raw) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 5
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': RT'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': RT (zscored per run) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 6
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': RT'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': RT (zscored per subject ie across all runs) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % net value chosen option (first regressor)
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).NV_chosen
                            case {0,1,2,3}
                            case {4,6}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': NVch-NVunch'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': NVch-NVunch '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 5 % p(chosen)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': p(chosen)'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': p(chosen) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % reward/punishment trial
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).R_vs_P
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': R-P'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': R vs P '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % binary variable indicating when choice = high effort option
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).choiceHighE
                            case 0
                            case {1,2}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': choice = highE'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': choice hE '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % reward variable option
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).R_varOption
                            case 0
                            case {1,3}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': R amount highE'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': R amount hE '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case {2,4}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': R level highE'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': R level hE '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % reward chosen option
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).R_chosen
                            case 0
                            case {1,3}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': R amount chosen'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': R amount chosen '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case {2,4}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': R level chosen'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': R level chosen '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % punishment variable option
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).P_varOption
                            case 0
                            case {1,3}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': P amount highE'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': P amount hE '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case {2,4}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': P level highE'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': P level hE '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % punishment chosen option
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).P_chosen
                            case 0
                            case {1,3}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': P amount chosen'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': P amount chosen '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case {2,4}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': P level chosen'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': P level chosen '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % money left
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_left
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': money left'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': money left (amounts) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': money left'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': |money left| (amounts) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': money left'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': money left (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': money left'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': |money left| (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % money right
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_right
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': money right'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': money right (amounts) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': money right'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': |money right| (amounts) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': money right'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': money right (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': money right'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': |money right| (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % money chosen
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_chosen
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': money chosen'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': money chosen (amounts) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': money chosen'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': |money chosen| (amounts) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': money chosen'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': money chosen (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': money chosen'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': |money chosen| (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % money unchosen
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_unchosen
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': money unchosen'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': money unchosen (amounts) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': money unchosen'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': |money unchosen| (amounts) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': money unchosen'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': money unchosen (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': money unchosen'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': |money unchosen| (amounts) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % money associated to the option which varies (the non-default
                        % option)
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_varOption
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': money hE'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': money hE option (amount) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': money hE'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': |money| hE (amount) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': money hE'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': money hE (level) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': money hE'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': |money| hE (level) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % money chosen - money unchosen
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_ch_min_unch
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': money ch-unch'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': money ch-unch '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % money chosen - money default
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_ch_min_fixOption
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': money ch-def'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': money chosen-default (amount) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': money ch-def'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': money chosen-default (level) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % money sum of both options
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_sum
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': money sum'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': money sum '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % effort chosen
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).E_chosen
                            case 0
                            case {1,4}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': effort chosen'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': effort chosen (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                error('case not ready yet');
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': effort-by-choice'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': effort-by-choice '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % effort unchosen
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).E_unchosen
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': effort unchosen'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': effort unchosen (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('case not ready yet');
                        end
                        
                        % effort hE
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).E_varOption
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': effort level hE'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': effort level hE '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('case not ready yet');
                        end
                        
                        % effort chosen - effort unchosen
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).E_ch_min_unch
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': effort ch-unch'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': effort ch-unch (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': effort ch-unch'];
                                switch task_id_nm
                                    case 'Ep'
                                        dispRegFn([num2str(n_regs),') chosen ',task_nm,': effort ch-unch (durations) '],dispRegs);
                                    case 'Em'
                                        dispRegFn([num2str(n_regs),') chosen ',task_nm,': effort ch-unch (nb answers to give) '],dispRegs);
                                end
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('case not ready yet');
                        end
                        
                        % effort chosen - fixed low effort option
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).E_ch_min_fixOption
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': effort ch-def'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': effort ch-def (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('case not ready yet');
                        end
                        
                        % effort sum of both options
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).E_sum
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': effort sum'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': effort sum (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': effort sum'];
                                switch task_id_nm
                                    case 'Ep'
                                        dispRegFn([num2str(n_regs),') chosen ',task_nm,': effort sum (durations) '],dispRegs);
                                    case 'Em'
                                        dispRegFn([num2str(n_regs),') chosen ',task_nm,': effort sum (nb answers to give) '],dispRegs);
                                end
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('case not ready yet');
                        end
                        
                        % (money)*(effort) for high effort option
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_level_x_E_varOption
                            case 0
                            case 1 % (money level)*(effort level) high E option
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': money x effort non-default'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': money x effort hE '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % (money)*(effort) for chosen option
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_level_x_E_chosen
                            case 0
                            case 1 % (money level)*(effort level) chosen option
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': money x effort chosen'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': money x effort chosen '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % (reward)*(effort) for high effort option
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).R_level_x_E_varOption
                            case 0
                            case 1 % (R level)*(effort level) high E option
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': R x E non-default'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': R x E hE '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % (reward)*(effort) for chosen option
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).R_level_x_E_chosen
                            case 0
                            case 1 % (R level)*(effort level) chosen option
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': R x E chosen'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': R x E chosen '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % (punishment)*(effort) for high effort option
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).P_level_x_E_varOption
                            case 0
                            case 1 % (P level)*(effort level) high E option
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': P x E non-default'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': P x E hE '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % (punishment)*(effort) for chosen option
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).P_level_x_E_chosen
                            case 0
                            case 1 % (P level)*(effort level) chosen option
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': P x E chosen'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': P x E chosen '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % net value chosen option
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).NV_chosen
                            case {0,4,5,6}
                            case {1,3}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': NVch-NVunch'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': NVch-NVunch '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2 % p(chosen)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': p(chosen)'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': p(chosen) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % net value variable option
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).NV_varOption
                            case 0
                            case 1 % NV(high E) - NV(low E)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': net value high E - low E'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': net value high E - low E '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2 % |NV(high E) - NV(low E)|
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': abs(NVhE-NVlE)'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': |net value high E - low E| '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3 % p(choice = high E)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': p(choice = high E)'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': p(choice = high E) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4 % NV(high E) - NV(low E) + bias
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': net value high E - low E + bias'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': net value high E - low E + bias '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 5 % |NV(high E) - NV(low E) + bias|
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': abs(NVhE-NVlE+bias)'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': |net value high E - low E + bias| '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % net value variable option bis
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).NV_varOption_bis
                            case 0
                            case 1 % NV(high E) - NV(low E)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': net value high E - low E'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': net value high E - low E '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2 % |NV(high E) - NV(low E)|
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': abs(NVhE-NVlE)'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': |net value high E - low E| '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3 % p(choice = high E)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': p(choice = high E)'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': p(choice = high E) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4 % NV(high E) - NV(low E) + bias
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': net value high E - low E + bias'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': net value high E - low E + bias '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 5 % |NV(high E) - NV(low E) + bias|
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': abs(NVhE-NVlE+bias)'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': |net value high E - low E + bias| '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        if strcmp(task_id_nm, 'Ep') % physical effort only
                            % area under the curve of the force that is gonna be produced
                            switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).F_integral
                                case 0
                                case {1,3,5,7}
                                    n_regs = n_regs + 1;
                                    reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': force integral'];
                                    dispRegFn([num2str(n_regs),') chosen ',task_nm,': effort integral '],dispRegs);
                                    % if derivative added => add derivatives
                                    n_regs = n_regs + add_drv;
                                case {2,4,6,8}
                                    n_regs = n_regs + 1;
                                    reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': force integral overshoot'];
                                    dispRegFn([num2str(n_regs),') chosen ',task_nm,': effort integral overshoot '],dispRegs);
                                    % if derivative added => add derivatives
                                    n_regs = n_regs + add_drv;
                                otherwise
                                    error('not ready yet');
                            end
                            
                            % fatigue
                            switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).fatigue
                                case 0
                                case {1,2}
                                    n_regs = n_regs + 1;
                                    reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': fatigue'];
                                    dispRegFn([num2str(n_regs),') chosen ',task_nm,': fatigue '],dispRegs);
                                    % if derivative added => add derivatives
                                    n_regs = n_regs + add_drv;
                                otherwise
                                    error('not ready yet');
                            end
                            
                            % Ech*fatigue
                            switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).Ech_x_fatigue
                                case 0
                                case 1
                                    n_regs = n_regs + 1;
                                    reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': Ech_x_fatigue'];
                                    dispRegFn([num2str(n_regs),') chosen ',task_nm,': Ech_x_fatigue '],dispRegs);
                                    % if derivative added => add derivatives
                                    n_regs = n_regs + add_drv;
                                otherwise
                                    error('not ready yet');
                            end
                        end % physical effort filter
                        
                        if strcmp(task_id_nm,'Em') % mental effort only
                            % efficacy of the next trial
                            switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).efficacy
                                case 0
                                case {1,2,3,4,5,6,7,8}
                                    n_regs = n_regs + 1;
                                    reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': efficacy'];
                                    dispRegFn([num2str(n_regs),') chosen ',task_nm,': efficacy '],dispRegs);
                                    % if derivative added => add derivatives
                                    n_regs = n_regs + add_drv;
                                otherwise
                                    error('not ready yet');
                            end
                            
                            % efficacy during the previous trial
                            switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).prevEfficacy
                                case 0
                                case {1,2,3,4,5,6,7,8}
                                    n_regs = n_regs + 1;
                                    reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': previous efficacy'];
                                    dispRegFn([num2str(n_regs),') chosen ',task_nm,': previous efficacy '],dispRegs);
                                    % if derivative added => add derivatives
                                    n_regs = n_regs + add_drv;
                                otherwise
                                    error('not ready yet');
                            end
                            
                            % (effort chosen)*(efficacy during the previous trial)
                            switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).Ech_x_prevEfficacy
                                case 0
                                case {1,2,3,4}
                                    n_regs = n_regs + 1;
                                    reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': (Ech)x(previous efficacy)'];
                                    dispRegFn([num2str(n_regs),') chosen ',task_nm,': (Ech)x(previous efficacy) '],dispRegs);
                                    % if derivative added => add derivatives
                                    n_regs = n_regs + add_drv;
                                otherwise
                                    error('not ready yet');
                            end
                        end % mental effort filter
                        
                        % trial number
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).trialN
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': trial number'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': trial number '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': (trial number)x(Echosen-Eunchosen) (E levels)'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': (trial number)x(effort chosen-effort unchosen) (effort levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': (trial number)x(EnonDef-Edef) (E levels)'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': (trial number)x(effort non-default - effort default) (effort levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': (trial number)x(EnonDef) (E levels)'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': (trial number)x(effort non-default) (effort levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % confidence
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).confidence
                            case 0
                            case 1 % confidence rating by the subjects
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': confidence'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': confidence (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case {2,3,4} % confidence inferred by the model
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': confidence'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': confidence (inferred by the model) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % RT (last regressor)
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).RT
                            case {0,4,5,6}
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': RT'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': RT (raw) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': RT'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': RT (zscored per run) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,' ',task_nm,': RT'];
                                dispRegFn([num2str(n_regs),') chosen ',task_nm,': RT (zscored per subject ie across all runs) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                    end % task loop
                case {3,4,5} % all pooled across sessions
                    
                    %% chosen onset
                    n_regs = n_regs + 1;
                    reg_names{n_regs} = ['ONSET chosen ',...
                        RP_dispChosen_nm,' ',splitE_dispChosen_nm];
                    dispRegFn([num2str(n_regs),') ONSET chosen option display ',...
                        RP_dispChosen_nm,' ',splitE_dispChosen_nm,': ',GLMprm.model_onset.(task_id_nm).chosen,' '],dispRegs);
                    % if derivative added => add derivatives
                    n_regs = n_regs + add_drv;
                    
                    %% chosen regressors
                    
                    % RT (first regressor)
                    switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).RT
                        case {0,1,2,3}
                        case 4
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': RT'];
                            dispRegFn([num2str(n_regs),') chosen: RT (raw) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 5
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': RT'];
                            dispRegFn([num2str(n_regs),') chosen: RT (zscored per run) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 6
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': RT'];
                            dispRegFn([num2str(n_regs),') chosen: RT (zscored per subject ie across all runs) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % net value chosen option (first regressor)
                    switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).NV_chosen
                        case {0,1,2,3}
                        case {4,6}
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': NVch-NVunch'];
                            dispRegFn([num2str(n_regs),') chosen: NVch-NVunch '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 5 % p(chosen)
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': p(chosen)'];
                            dispRegFn([num2str(n_regs),') chosen: p(chosen) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % reward/punishment trial
                    switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).R_vs_P
                        case 0
                        case 1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': R-P'];
                            dispRegFn([num2str(n_regs),') chosen: R vs P '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % binary variable indicating when choice = high effort option
                    switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).choiceHighE
                        case 0
                        case {1,2}
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': choice = highE'];
                            dispRegFn([num2str(n_regs),') chosen: choice hE '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % reward variable option
                    switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).R_varOption
                        case 0
                        case {1,3}
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': R amount highE'];
                            dispRegFn([num2str(n_regs),') chosen: R amount hE '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case {2,4}
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': R level highE'];
                            dispRegFn([num2str(n_regs),') chosen: R level hE '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % reward chosen option
                    switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).R_chosen
                        case 0
                        case {1,3}
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': R amount chosen'];
                            dispRegFn([num2str(n_regs),') chosen: R amount chosen '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case {2,4}
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': R level chosen'];
                            dispRegFn([num2str(n_regs),') chosen: R level chosen '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % punishment variable option
                    switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).P_varOption
                        case 0
                        case {1,3}
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': P amount highE'];
                            dispRegFn([num2str(n_regs),') chosen: P amount hE '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case {2,4}
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': P level highE'];
                            dispRegFn([num2str(n_regs),') chosen: P level hE '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % punishment chosen option
                    switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).P_chosen
                        case 0
                        case {1,3}
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': P amount chosen'];
                            dispRegFn([num2str(n_regs),') chosen: P amount chosen '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case {2,4}
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': P level chosen'];
                            dispRegFn([num2str(n_regs),') chosen: P level chosen '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % money left
                    switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_left
                        case 0
                        case 1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money left'];
                            dispRegFn([num2str(n_regs),') chosen: money left (amounts) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money left'];
                            dispRegFn([num2str(n_regs),') chosen: |money left| (amounts) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 3
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money left'];
                            dispRegFn([num2str(n_regs),') chosen: money left (levels) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 4
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money left'];
                            dispRegFn([num2str(n_regs),') chosen: |money left| (levels) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % money right
                    switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_right
                        case 0
                        case 1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money right'];
                            dispRegFn([num2str(n_regs),') chosen: money right (amounts) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money right'];
                            dispRegFn([num2str(n_regs),') chosen: |money right| (amounts) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 3
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money right'];
                            dispRegFn([num2str(n_regs),') chosen: money right (levels) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 4
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money right'];
                            dispRegFn([num2str(n_regs),') chosen: |money right| (levels) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                    end
                    
                    % money chosen
                    switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_chosen
                        case 0
                        case 1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money chosen'];
                            dispRegFn([num2str(n_regs),') chosen: money chosen (amounts) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money chosen'];
                            dispRegFn([num2str(n_regs),') chosen: |money chosen| (amounts) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 3
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money chosen'];
                            dispRegFn([num2str(n_regs),') chosen: money chosen (levels) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 4
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money chosen'];
                            dispRegFn([num2str(n_regs),') chosen: |money chosen| (levels) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                    end
                    
                    % money unchosen
                    switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_unchosen
                        case 0
                        case 1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money unchosen'];
                            dispRegFn([num2str(n_regs),') chosen: money unchosen (amounts) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money unchosen'];
                            dispRegFn([num2str(n_regs),') chosen: |money unchosen| (amounts) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 3
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money unchosen'];
                            dispRegFn([num2str(n_regs),') chosen: money unchosen (levels) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 4
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money unchosen'];
                            dispRegFn([num2str(n_regs),') chosen: |money unchosen| (amounts) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                    end
                    
                    % money associated to the option which varies (the non-default
                    % option)
                    switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_varOption
                        case 0
                        case 1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money hE'];
                            dispRegFn([num2str(n_regs),') chosen: money hE option (amount) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money hE'];
                            dispRegFn([num2str(n_regs),') chosen: |money| hE (amount) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 3
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money hE'];
                            dispRegFn([num2str(n_regs),') chosen: money hE (level) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 4
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money hE'];
                            dispRegFn([num2str(n_regs),') chosen: |money| hE (level) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                    end
                    
                    % money chosen - money unchosen
                    switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_ch_min_unch
                        case 0
                        case 1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money ch-unch'];
                            dispRegFn([num2str(n_regs),') chosen: money ch-unch '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                    end
                    
                    % money chosen - money default
                    switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_ch_min_fixOption
                        case 0
                        case 1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money ch-def'];
                            dispRegFn([num2str(n_regs),') chosen: money chosen-default (amount) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money ch-def'];
                            dispRegFn([num2str(n_regs),') chosen: money chosen-default (level) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                    end
                    
                    % money sum of both options
                    switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_sum
                        case 0
                        case 1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money sum'];
                            dispRegFn([num2str(n_regs),') chosen: money sum '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                    end
                    
                    % effort chosen
                    switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).E_chosen
                        case 0
                        case {1,4}
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': effort chosen'];
                            dispRegFn([num2str(n_regs),') chosen: effort chosen (levels) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2
                            error('case not ready yet');
                        case 3
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': effort-by-choice'];
                            dispRegFn([num2str(n_regs),') chosen: effort-by-choice '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                    end
                    
                    % effort unchosen
                    switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).E_unchosen
                        case 0
                        case 1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': effort unchosen'];
                            dispRegFn([num2str(n_regs),') chosen: effort unchosen (levels) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('case not ready yet');
                    end
                    
                    % effort hE
                    switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).E_varOption
                        case 0
                        case 1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': effort level hE'];
                            dispRegFn([num2str(n_regs),') chosen: effort level hE '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('case not ready yet');
                    end
                    
                    % effort chosen - effort unchosen
                    switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).E_ch_min_unch
                        case 0
                        case 1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': effort ch-unch'];
                            dispRegFn([num2str(n_regs),') chosen: effort ch-unch (levels) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': effort ch-unch'];
                            switch task_id_nm
                                case 'Ep'
                                    dispRegFn([num2str(n_regs),') chosen: effort ch-unch (durations) '],dispRegs);
                                case 'Em'
                                    dispRegFn([num2str(n_regs),') chosen: effort ch-unch (nb answers to give) '],dispRegs);
                            end
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('case not ready yet');
                    end
                    
                    % effort chosen - fixed low effort option
                    switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).E_ch_min_fixOption
                        case 0
                        case 1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': effort ch-def'];
                            dispRegFn([num2str(n_regs),') chosen: effort ch-def (levels) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('case not ready yet');
                    end
                    
                    % effort sum of both options
                    switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).E_sum
                        case 0
                        case 1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': effort sum'];
                            dispRegFn([num2str(n_regs),') chosen: effort sum (levels) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': effort sum'];
                            switch task_id_nm
                                case 'Ep'
                                    dispRegFn([num2str(n_regs),') chosen: effort sum (durations) '],dispRegs);
                                case 'Em'
                                    dispRegFn([num2str(n_regs),') chosen: effort sum (nb answers to give) '],dispRegs);
                            end
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('case not ready yet');
                    end
                    
                    % (money)*(effort) for high effort option
                    switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_level_x_E_varOption
                        case 0
                        case 1 % (money level)*(effort level) high E option
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money x effort non-default'];
                            dispRegFn([num2str(n_regs),') chosen: money x effort hE '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % (money)*(effort) for chosen option
                    switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_level_x_E_chosen
                        case 0
                        case 1 % (money level)*(effort level) chosen option
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money x effort chosen'];
                            dispRegFn([num2str(n_regs),') chosen: money x effort chosen '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % (reward)*(effort) for high effort option
                    switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).R_level_x_E_varOption
                        case 0
                        case 1 % (R level)*(effort level) high E option
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': R x E non-default'];
                            dispRegFn([num2str(n_regs),') chosen: R x E hE '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % (reward)*(effort) for chosen option
                    switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).R_level_x_E_chosen
                        case 0
                        case 1 % (R level)*(effort level) chosen option
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': R x E chosen'];
                            dispRegFn([num2str(n_regs),') chosen: R x E chosen '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % (punishment)*(effort) for high effort option
                    switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).P_level_x_E_varOption
                        case 0
                        case 1 % (P level)*(effort level) high E option
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': P x E non-default'];
                            dispRegFn([num2str(n_regs),') chosen: P x E hE '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % (punishment)*(effort) for chosen option
                    switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).P_level_x_E_chosen
                        case 0
                        case 1 % (P level)*(effort level) chosen option
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': P x E chosen'];
                            dispRegFn([num2str(n_regs),') chosen: P x E chosen '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % net value chosen option
                    switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).NV_chosen
                        case {0,4,5,6}
                        case {1,3}
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': NVch-NVunch'];
                            dispRegFn([num2str(n_regs),') chosen: NVch-NVunch '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2 % p(chosen)
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': p(chosen)'];
                            dispRegFn([num2str(n_regs),') chosen: p(chosen) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % net value variable option
                    switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).NV_varOption
                        case 0
                        case 1 % NV(high E) - NV(low E)
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': net value high E - low E'];
                            dispRegFn([num2str(n_regs),') chosen: net value high E - low E '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2 % |NV(high E) - NV(low E)|
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': abs(NVhE-NVlE)'];
                            dispRegFn([num2str(n_regs),') chosen: |net value high E - low E| '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 3 % p(choice = high E)
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': p(choice = high E)'];
                            dispRegFn([num2str(n_regs),') chosen: p(choice = high E) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 4 % NV(high E) - NV(low E) + bias
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': net value high E - low E + bias'];
                            dispRegFn([num2str(n_regs),') chosen: net value high E - low E + bias '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 5 % |NV(high E) - NV(low E) + bias|
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': abs(NVhE-NVlE+bias)'];
                            dispRegFn([num2str(n_regs),') chosen: |net value high E - low E + bias| '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % net value variable option bis
                    switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).NV_varOption_bis
                        case 0
                        case 1 % NV(high E) - NV(low E)
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': net value high E - low E'];
                            dispRegFn([num2str(n_regs),') chosen: net value high E - low E '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2 % |NV(high E) - NV(low E)|
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': abs(NVhE-NVlE)'];
                            dispRegFn([num2str(n_regs),') chosen: |net value high E - low E| '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 3 % p(choice = high E)
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': p(choice = high E)'];
                            dispRegFn([num2str(n_regs),') chosen: p(choice = high E) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 4 % NV(high E) - NV(low E) + bias
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': net value high E - low E + bias'];
                            dispRegFn([num2str(n_regs),') chosen: net value high E - low E + bias '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 5 % |NV(high E) - NV(low E) + bias|
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': abs(NVhE-NVlE+bias)'];
                            dispRegFn([num2str(n_regs),') chosen: |net value high E - low E + bias| '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    if strcmp(task_id_nm, 'Ep') % physical effort only
                        % area under the curve of the force that is gonna be produced
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).F_integral
                            case 0
                            case {1,3,5,7}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': force integral'];
                                dispRegFn([num2str(n_regs),') chosen: effort integral '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case {2,4,6,8}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': force integral overshoot'];
                                dispRegFn([num2str(n_regs),') chosen: effort integral overshoot '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % fatigue
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).fatigue
                            case 0
                            case {1,2}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': fatigue'];
                                dispRegFn([num2str(n_regs),') chosen: fatigue '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % Ech*fatigue
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).Ech_x_fatigue
                            case 0
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': Ech_x_fatigue'];
                                dispRegFn([num2str(n_regs),') chosen: Ech_x_fatigue '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                    end % physical effort filter
                    
                    if strcmp(task_id_nm,'Em') % mental effort only
                        % efficacy of the next trial
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).efficacy
                            case 0
                            case {1,2,3,4,5,6,7,8}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': efficacy'];
                                dispRegFn([num2str(n_regs),') chosen: efficacy '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % efficacy during the previous trial
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).prevEfficacy
                            case 0
                            case {1,2,3,4,5,6,7,8}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': previous efficacy'];
                                dispRegFn([num2str(n_regs),') chosen: previous efficacy '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % (effort chosen)*(efficacy during the previous trial)
                        switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).Ech_x_prevEfficacy
                            case 0
                            case {1,2,3,4}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': (Ech)x(previous efficacy)'];
                                dispRegFn([num2str(n_regs),') chosen: (Ech)x(previous efficacy) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                    end % mental effort filter
                    
                    % trial number
                    switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).trialN
                        case 0
                        case 1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': trial number'];
                            dispRegFn([num2str(n_regs),') chosen: trial number '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': (trial number)x(Echosen-Eunchosen) (E levels)'];
                            dispRegFn([num2str(n_regs),') chosen: (trial number)x(effort chosen-effort unchosen) (effort levels) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 3
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': (trial number)x(EnonDef-Edef) (E levels)'];
                            dispRegFn([num2str(n_regs),') chosen: (trial number)x(effort non-default - effort default) (effort levels) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 4
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': (trial number)x(EnonDef) (E levels)'];
                            dispRegFn([num2str(n_regs),') chosen: (trial number)x(effort non-default) (effort levels) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % confidence
                    switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).confidence
                        case 0
                        case 1 % confidence rating by the subjects
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': confidence'];
                            dispRegFn([num2str(n_regs),') chosen: confidence (levels) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case {2,3,4} % confidence inferred by the model
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': confidence'];
                            dispRegFn([num2str(n_regs),') chosen: confidence (inferred by the model) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % RT (last regressor)
                    switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).RT
                        case {0,4,5,6}
                        case 1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': RT'];
                            dispRegFn([num2str(n_regs),') chosen: RT (raw) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': RT'];
                            dispRegFn([num2str(n_regs),') chosen: RT (zscored per run) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 3
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': RT'];
                            dispRegFn([num2str(n_regs),') chosen: RT (zscored per subject ie across all runs) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
            end % DCM mode
        end % loop effort level
    end % loop reward/punishment
end % chosen onset

%% pre-effort fixation cross
if ~strcmp(GLMprm.model_onset.(task_id_nm).preEffortCross,'none')
    
    % check if trials are split or not
    if GLMprm.preEffortCross.(task_id_nm).RPpool == 1 % pool reward and punishment trials
        RP_preEcross = {'RP'};
    elseif GLMprm.preEffortCross.(task_id_nm).RPpool == 0 % split reward and punishment
        RP_preEcross = {'R','P'};
    end % RP pool
    n_RP_preEcross = length(RP_preEcross);
    
    % check if trials are split or not according to Effort levels
    switch GLMprm.preEffortCross.(task_id_nm).splitPerE
        case 0 % pool trials
            splitE_preEcross = {'E'};
        case 1 % split according to effort proposed
            splitE_preEcross = {'E1','E2','E3'};
        case 2 % split according to effort chosen
            splitE_preEcross = {'Ech0','Ech1','Ech2','Ech3'};
        case 3 % split according to option chosen (low/high effort)
            splitE_preEcross = {'lEch','hEch'};
    end % Effort level pool
    n_splitE_preEcross = length(splitE_preEcross);
    
    % loop through conditions for choice period
    for iRP_preEcross = 1:n_RP_preEcross
        RP_preEcross_nm = RP_preEcross{iRP_preEcross};
        
        for iE_preEcross = 1:n_splitE_preEcross
            splitE_preEcross_nm = splitE_preEcross{iE_preEcross};
            
            switch DCM_mode
                case 1 % all sessions independent
                    for iRun = 1:n_runs
                        run_nb_nm = ['r',num2str(runs.runsToKeep(iRun))];
                        
                        %% effort period onset
                        n_regs = n_regs + 1;
                        reg_names{n_regs} = ['ONSET preEffort black cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',run_nb_nm];
                        dispRegFn([num2str(n_regs),') ONSET preEffort black cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',run_nb_nm,': ',...
                            GLMprm.model_onset.(task_id_nm).preEffortCross,' '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs = n_regs + add_drv;
                        
                        %% effort period regressors
                        % binary variable indicating when choice = high effort option
                        switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).choiceHighE
                            case 0
                            case {1,2}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',run_nb_nm,': choice = highE'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',run_nb_nm,': choice hE '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % money chosen
                        switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).money_chosen
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',run_nb_nm,': money chosen'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',run_nb_nm,': money chosen (amounts) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',run_nb_nm,': money chosen'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',run_nb_nm,': |money chosen| (amounts) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',run_nb_nm,': money chosen'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',run_nb_nm,': money chosen (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',run_nb_nm,': money chosen'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',run_nb_nm,': |money chosen| (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % effort chosen
                        switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).E_chosen
                            case 0
                            case {1,4}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',run_nb_nm,': effort chosen'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',run_nb_nm,': effort chosen (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                error('case not ready yet');
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',run_nb_nm,': effort-by-choice'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',run_nb_nm,': effort-by-choice '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('case not ready yet');
                        end
                        
                        % force peak
                        switch task_id_nm
                            case 'Ep'
                                switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).F_peak
                                    case {1,2,3,4}
                                        n_regs = n_regs + 1;
                                        reg_names{n_regs} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',run_nb_nm,': force peak'];
                                        dispRegFn([num2str(n_regs),') pre-effort cross ',run_nb_nm,': force peak '],dispRegs);
                                        % if derivative added => add derivatives
                                        n_regs = n_regs + add_drv;
                                end
                                
                                % force integral
                                switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).F_integral
                                    case {1,3,5,7}
                                        n_regs = n_regs + 1;
                                        reg_names{n_regs} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',run_nb_nm,': force integral'];
                                        dispRegFn([num2str(n_regs),') pre-effort cross ',run_nb_nm,': force integral '],dispRegs);
                                        % if derivative added => add derivatives
                                        n_regs = n_regs + add_drv;
                                    case {2,4,6,8}
                                        error('please update');
                                end
                                
                            case 'Em'
                                % average RT for all numbers of each trial
                                switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).RT_avg
                                    case 1
                                        n_regs = n_regs + 1;
                                        reg_names{n_regs} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',run_nb_nm,': average RT effort'];
                                        dispRegFn([num2str(n_regs),') pre-effort cross ',run_nb_nm,': average RT effort '],dispRegs);
                                        % if derivative added => add derivatives
                                        n_regs = n_regs + add_drv;
                                end
                                
                                % number of correct answers
                                switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).n_correct
                                    case 1
                                        n_regs = n_regs + 1;
                                        reg_names{n_regs} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',run_nb_nm,': number of correct answers'];
                                        dispRegFn([num2str(n_regs),') pre-effort cross ',run_nb_nm,': number of correct answers '],dispRegs);
                                        % if derivative added => add derivatives
                                        n_regs = n_regs + add_drv;
                                end
                                
                                % number of errors
                                switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).n_errors
                                    case 1
                                        n_regs = n_regs + 1;
                                        reg_names{n_regs} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',run_nb_nm,': number of errors'];
                                        dispRegFn([num2str(n_regs),') pre-effort cross ',run_nb_nm,': number of errors '],dispRegs);
                                        % if derivative added => add derivatives
                                        n_regs = n_regs + add_drv;
                                end
                        end
                        
                        % net value chosen option
                        switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).NV_chosen
                            case 0
                            case {1,3}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',run_nb_nm,': NVch-NVunch'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',run_nb_nm,': NVch-NVunch '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',run_nb_nm,': p(chosen)'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',run_nb_nm,': p(chosen) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('problem with pre-E cross NV_chosen value');
                        end
                        
                        % net value hE
                        switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).NV_varOption
                            case 1 % NV(high E) - NV(low E)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',run_nb_nm,': net value high E - low E'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',run_nb_nm,': net value high E - low E '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2 % |NV(high E) - NV(low E)|
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',run_nb_nm,': abs(NVhE-NVlE)'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',run_nb_nm,': |net value high E - low E| '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3 % p(choice = high E)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',run_nb_nm,': p(choice = high E)'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',run_nb_nm,': p(choice = high E) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4 % NV(high E) - NV(low E) + bias
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',run_nb_nm,': net value high E - low E + bias'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',run_nb_nm,': net value high E - low E + bias '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 5 % |NV(high E) - NV(low E) + bias|
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',run_nb_nm,': abs(NVhE-NVlE+bias)'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',run_nb_nm,': |net value high E - low E + bias| '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % net value hE bis
                        switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).NV_varOption_bis
                            case 1 % NV(high E) - NV(low E)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',run_nb_nm,': net value high E - low E'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',run_nb_nm,': net value high E - low E '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2 % |NV(high E) - NV(low E)|
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',run_nb_nm,': abs(NVhE-NVlE)'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',run_nb_nm,': |net value high E - low E| '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3 % p(choice = high E)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',run_nb_nm,': p(choice = high E)'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',run_nb_nm,': p(choice = high E) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4 % NV(high E) - NV(low E) + bias
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',run_nb_nm,': net value high E - low E + bias'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',run_nb_nm,': net value high E - low E + bias '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 5 % |NV(high E) - NV(low E) + bias|
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',run_nb_nm,': abs(NVhE-NVlE+bias)'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',run_nb_nm,': |net value high E - low E + bias| '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % RT first answer
                        switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).RT_1stAnswer
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',run_nb_nm,': RT 1st answer'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',run_nb_nm,': RT 1st answer (raw) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % trial number
                        switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).trialN
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',run_nb_nm,': trial number'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',run_nb_nm,': trial number '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',run_nb_nm,': (trial number)x(Echosen-Eunchosen) (E levels)'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',run_nb_nm,': (trial number)x(effort chosen-effort unchosen) (effort levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',run_nb_nm,': (trial number)x(EnonDef-Edef) (E levels)'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',run_nb_nm,': (trial number)x(effort non-default - effort default) (effort levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',run_nb_nm,': (trial number)x(EnonDef) (E levels)'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',run_nb_nm,': (trial number)x(effort non-default) (effort levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % confidence
                        switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).confidence
                            case 1 % confidence ratings 0/1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',run_nb_nm,': confidence'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',run_nb_nm,': confidence (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case {2,3,4} % confidence inferred by the model
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',run_nb_nm,': confidence'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',run_nb_nm,': confidence (inferred by the model) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                    end % run loop
                case 2 % all tasks independent but sessions pooled
                    for iTask = 1:nTasks
                        task_nm = tasks{iTask};
                        
                        %% effort period onset
                        n_regs = n_regs + 1;
                        reg_names{n_regs} = ['ONSET preEffort black cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',task_nm];
                        dispRegFn([num2str(n_regs),') ONSET preEffort black cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',task_nm,': ',...
                            GLMprm.model_onset.(task_id_nm).preEffortCross,' '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs = n_regs + add_drv;
                        
                        %% effort period regressors
                        % binary variable indicating when choice = high effort option
                        switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).choiceHighE
                            case 0
                            case {1,2}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',task_nm,': choice = highE'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',task_nm,': choice hE '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % money chosen
                        switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).money_chosen
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',task_nm,': money chosen'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',task_nm,': money chosen (amounts) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',task_nm,': money chosen'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',task_nm,': |money chosen| (amounts) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',task_nm,': money chosen'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',task_nm,': money chosen (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',task_nm,': money chosen'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',task_nm,': |money chosen| (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % effort chosen
                        switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).E_chosen
                            case 0
                            case {1,4}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',task_nm,': effort chosen'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',task_nm,': effort chosen (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                error('case not ready yet');
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',task_nm,': effort-by-choice'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',task_nm,': effort-by-choice '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('case not ready yet');
                        end
                        
                        % force peak
                        switch task_id_nm
                            case 'Ep'
                                switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).F_peak
                                    case {1,2,3,4}
                                        n_regs = n_regs + 1;
                                        reg_names{n_regs} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',task_nm,': force peak'];
                                        dispRegFn([num2str(n_regs),') pre-effort cross ',task_nm,': force peak '],dispRegs);
                                        % if derivative added => add derivatives
                                        n_regs = n_regs + add_drv;
                                end
                                
                                % force integral
                                switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).F_integral
                                    case {1,3,5,7}
                                        n_regs = n_regs + 1;
                                        reg_names{n_regs} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',task_nm,': force integral'];
                                        dispRegFn([num2str(n_regs),') pre-effort cross ',task_nm,': force integral '],dispRegs);
                                        % if derivative added => add derivatives
                                        n_regs = n_regs + add_drv;
                                    case {2,4,6,8}
                                        error('please update');
                                end
                                
                            case 'Em'
                                % average RT for all numbers of each trial
                                switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).RT_avg
                                    case 1
                                        n_regs = n_regs + 1;
                                        reg_names{n_regs} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',task_nm,': average RT effort'];
                                        dispRegFn([num2str(n_regs),') pre-effort cross ',task_nm,': average RT effort '],dispRegs);
                                        % if derivative added => add derivatives
                                        n_regs = n_regs + add_drv;
                                end
                                
                                % number of correct answers
                                switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).n_correct
                                    case 1
                                        n_regs = n_regs + 1;
                                        reg_names{n_regs} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',task_nm,': number of correct answers'];
                                        dispRegFn([num2str(n_regs),') pre-effort cross ',task_nm,': number of correct answers '],dispRegs);
                                        % if derivative added => add derivatives
                                        n_regs = n_regs + add_drv;
                                end
                                
                                % number of errors
                                switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).n_errors
                                    case 1
                                        n_regs = n_regs + 1;
                                        reg_names{n_regs} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',task_nm,': number of errors'];
                                        dispRegFn([num2str(n_regs),') pre-effort cross ',task_nm,': number of errors '],dispRegs);
                                        % if derivative added => add derivatives
                                        n_regs = n_regs + add_drv;
                                end
                        end
                        
                        % net value chosen option
                        switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).NV_chosen
                            case 0
                            case {1,3}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',task_nm,': NVch-NVunch'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',task_nm,': NVch-NVunch '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',task_nm,': p(chosen)'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',task_nm,': p(chosen) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('problem with pre-E cross NV_chosen value');
                        end
                        
                        % net value hE
                        switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).NV_varOption
                            case 1 % NV(high E) - NV(low E)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',task_nm,': net value high E - low E'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',task_nm,': net value high E - low E '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2 % |NV(high E) - NV(low E)|
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',task_nm,': abs(NVhE-NVlE)'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',task_nm,': |net value high E - low E| '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3 % p(choice = high E)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',task_nm,': p(choice = high E)'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',task_nm,': p(choice = high E) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4 % NV(high E) - NV(low E) + bias
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',task_nm,': net value high E - low E + bias'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',task_nm,': net value high E - low E + bias '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 5 % |NV(high E) - NV(low E) + bias|
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',task_nm,': abs(NVhE-NVlE+bias)'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',task_nm,': |net value high E - low E + bias| '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % net value hE bis
                        switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).NV_varOption_bis
                            case 1 % NV(high E) - NV(low E)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',task_nm,': net value high E - low E'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',task_nm,': net value high E - low E '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2 % |NV(high E) - NV(low E)|
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',task_nm,': abs(NVhE-NVlE)'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',task_nm,': |net value high E - low E| '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3 % p(choice = high E)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',task_nm,': p(choice = high E)'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',task_nm,': p(choice = high E) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4 % NV(high E) - NV(low E) + bias
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',task_nm,': net value high E - low E + bias'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',task_nm,': net value high E - low E + bias '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 5 % |NV(high E) - NV(low E) + bias|
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',task_nm,': abs(NVhE-NVlE+bias)'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',task_nm,': |net value high E - low E + bias| '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % RT first answer
                        switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).RT_1stAnswer
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',task_nm,': RT 1st answer'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',task_nm,': RT 1st answer (raw) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % trial number
                        switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).trialN
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',task_nm,': trial number'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',task_nm,': trial number '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',task_nm,': (trial number)x(Echosen-Eunchosen) (E levels)'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',task_nm,': (trial number)x(effort chosen-effort unchosen) (effort levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',task_nm,': (trial number)x(EnonDef-Edef) (E levels)'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',task_nm,': (trial number)x(effort non-default - effort default) (effort levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',task_nm,': (trial number)x(EnonDef) (E levels)'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',task_nm,': (trial number)x(effort non-default) (effort levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % confidence
                        switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).confidence
                            case 1 % confidence ratings 0/1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',task_nm,': confidence'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',task_nm,': confidence (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case {2,3,4} % confidence inferred by the model
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,' ',task_nm,': confidence'];
                                dispRegFn([num2str(n_regs),') pre-effort cross ',task_nm,': confidence (inferred by the model) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                    end % task loop
                case {3,4,5} % all pooled across sessions
                    
                    %% effort period onset
                    n_regs = n_regs + 1;
                    reg_names{n_regs} = ['ONSET preEffort black cross ',RP_preEcross_nm,' ',splitE_preEcross_nm];
                    dispRegFn([num2str(n_regs),') ONSET preEffort black cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,': ',...
                        GLMprm.model_onset.(task_id_nm).preEffortCross,' '],dispRegs);
                    % if derivative added => add derivatives
                    n_regs = n_regs + add_drv;
                    
                    %% effort period regressors
                    % binary variable indicating when choice = high effort option
                    switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).choiceHighE
                        case 0
                        case {1,2}
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,': choice = highE'];
                            dispRegFn([num2str(n_regs),') pre-effort cross: choice hE '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % money chosen
                    switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).money_chosen
                        case 1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,': money chosen'];
                            dispRegFn([num2str(n_regs),') pre-effort cross: money chosen (amounts) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,': money chosen'];
                            dispRegFn([num2str(n_regs),') pre-effort cross: |money chosen| (amounts) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 3
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,': money chosen'];
                            dispRegFn([num2str(n_regs),') pre-effort cross: money chosen (levels) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 4
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,': money chosen'];
                            dispRegFn([num2str(n_regs),') pre-effort cross: |money chosen| (levels) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                    end
                    
                    % effort chosen
                    switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).E_chosen
                        case 0
                        case {1,4}
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,': effort chosen'];
                            dispRegFn([num2str(n_regs),') pre-effort cross: effort chosen (levels) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2
                            error('case not ready yet');
                        case 3
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,': effort-by-choice'];
                            dispRegFn([num2str(n_regs),') pre-effort cross: effort-by-choice '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('case not ready yet');
                    end
                    
                    % force peak
                    switch task_id_nm
                        case 'Ep'
                            switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).F_peak
                                case {1,2,3,4}
                                    n_regs = n_regs + 1;
                                    reg_names{n_regs} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,': force peak'];
                                    dispRegFn([num2str(n_regs),') pre-effort cross: force peak '],dispRegs);
                                    % if derivative added => add derivatives
                                    n_regs = n_regs + add_drv;
                            end
                            
                            % force integral
                            switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).F_integral
                                case {1,3,5,7}
                                    n_regs = n_regs + 1;
                                    reg_names{n_regs} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,': force integral'];
                                    dispRegFn([num2str(n_regs),') pre-effort cross: force integral '],dispRegs);
                                    % if derivative added => add derivatives
                                    n_regs = n_regs + add_drv;
                                case {2,4,6,8}
                                    error('please update');
                            end
                            
                        case 'Em'
                            % average RT for all numbers of each trial
                            switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).RT_avg
                                case 1
                                    n_regs = n_regs + 1;
                                    reg_names{n_regs} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,': average RT effort'];
                                    dispRegFn([num2str(n_regs),') pre-effort cross: average RT effort '],dispRegs);
                                    % if derivative added => add derivatives
                                    n_regs = n_regs + add_drv;
                            end
                            
                            % number of correct answers
                            switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).n_correct
                                case 1
                                    n_regs = n_regs + 1;
                                    reg_names{n_regs} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,': number of correct answers'];
                                    dispRegFn([num2str(n_regs),') pre-effort cross: number of correct answers '],dispRegs);
                                    % if derivative added => add derivatives
                                    n_regs = n_regs + add_drv;
                            end
                            
                            % number of errors
                            switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).n_errors
                                case 1
                                    n_regs = n_regs + 1;
                                    reg_names{n_regs} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,': number of errors'];
                                    dispRegFn([num2str(n_regs),') pre-effort cross: number of errors '],dispRegs);
                                    % if derivative added => add derivatives
                                    n_regs = n_regs + add_drv;
                            end
                    end
                    
                    % net value chosen option
                    switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).NV_chosen
                        case 0
                        case {1,3}
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,': NVch-NVunch'];
                            dispRegFn([num2str(n_regs),') pre-effort cross: NVch-NVunch '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,': p(chosen)'];
                            dispRegFn([num2str(n_regs),') pre-effort cross: p(chosen) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('problem with pre-E cross NV_chosen value');
                    end
                    
                    % net value hE
                    switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).NV_varOption
                        case 1 % NV(high E) - NV(low E)
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,': net value high E - low E'];
                            dispRegFn([num2str(n_regs),') pre-effort cross: net value high E - low E '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2 % |NV(high E) - NV(low E)|
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,': abs(NVhE-NVlE)'];
                            dispRegFn([num2str(n_regs),') pre-effort cross: |net value high E - low E| '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 3 % p(choice = high E)
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,': p(choice = high E)'];
                            dispRegFn([num2str(n_regs),') pre-effort cross: p(choice = high E) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 4 % NV(high E) - NV(low E) + bias
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,': net value high E - low E + bias'];
                            dispRegFn([num2str(n_regs),') pre-effort cross: net value high E - low E + bias '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 5 % |NV(high E) - NV(low E) + bias|
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,': abs(NVhE-NVlE+bias)'];
                            dispRegFn([num2str(n_regs),') pre-effort cross: |net value high E - low E + bias| '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                    end
                    
                    % net value hE bis
                    switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).NV_varOption_bis
                        case 1 % NV(high E) - NV(low E)
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,': net value high E - low E'];
                            dispRegFn([num2str(n_regs),') pre-effort cross: net value high E - low E '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2 % |NV(high E) - NV(low E)|
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,': abs(NVhE-NVlE)'];
                            dispRegFn([num2str(n_regs),') pre-effort cross: |net value high E - low E| '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 3 % p(choice = high E)
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,': p(choice = high E)'];
                            dispRegFn([num2str(n_regs),') pre-effort cross: p(choice = high E) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 4 % NV(high E) - NV(low E) + bias
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,': net value high E - low E + bias'];
                            dispRegFn([num2str(n_regs),') pre-effort cross: net value high E - low E + bias '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 5 % |NV(high E) - NV(low E) + bias|
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,': abs(NVhE-NVlE+bias)'];
                            dispRegFn([num2str(n_regs),') pre-effort cross: |net value high E - low E + bias| '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                    end
                    
                    % RT first answer
                    switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).RT_1stAnswer
                        case 1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,': RT 1st answer'];
                            dispRegFn([num2str(n_regs),') pre-effort cross: RT 1st answer (raw) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                    end
                    
                    % trial number
                    switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).trialN
                        case 1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,': trial number'];
                            dispRegFn([num2str(n_regs),') pre-effort cross: trial number '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,': (trial number)x(Echosen-Eunchosen) (E levels)'];
                            dispRegFn([num2str(n_regs),') pre-effort cross: (trial number)x(effort chosen-effort unchosen) (effort levels) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 3
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,': (trial number)x(EnonDef-Edef) (E levels)'];
                            dispRegFn([num2str(n_regs),') pre-effort cross: (trial number)x(effort non-default - effort default) (effort levels) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 4
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,': (trial number)x(EnonDef) (E levels)'];
                            dispRegFn([num2str(n_regs),') pre-effort cross: (trial number)x(effort non-default) (effort levels) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                    end
                    
                    % confidence
                    switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).confidence
                        case 1 % confidence ratings 0/1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,': confidence'];
                            dispRegFn([num2str(n_regs),') pre-effort cross: confidence (levels) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case {2,3,4} % confidence inferred by the model
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,': confidence'];
                            dispRegFn([num2str(n_regs),') pre-effort cross: confidence (inferred by the model) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                    end
                    
            end % DCM mode
        end % effort level loop
    end % loop reward/punishment
end % pre-effort cross

%% effort performance period
if ~strcmp(GLMprm.model_onset.(task_id_nm).Eperf,'none')
    % check if trials are split or not
    if GLMprm.Eperf.(task_id_nm).RPpool == 1 % pool reward and punishment trials
        RP_Eperf = {'RP'};
    elseif GLMprm.Eperf.(task_id_nm).RPpool == 0 % split reward and punishment
        RP_Eperf = {'R','P'};
    end % RP pool
    n_RP_Eperf = length(RP_Eperf);
    
    % check if trials are split or not according to Effort levels
    switch GLMprm.Eperf.(task_id_nm).splitPerE
        case 0 % pool trials
            splitE_Eperf = {'E'};
        case 1 % split according to effort proposed
            splitE_Eperf = {'E1','E2','E3'};
        case 2 % split according to effort chosen
            splitE_Eperf = {'Ech0','Ech1','Ech2','Ech3'};
        case 3 % split according to option chosen (low/high effort)
            splitE_Eperf = {'lEch','hEch'};
    end % Effort level pool
    n_splitE_Eperf = length(splitE_Eperf);
    
    % loop through conditions for choice period
    for iRP_Eperf = 1:n_RP_Eperf
        RP_Eperf_nm = RP_Eperf{iRP_Eperf};
        
        for iE_Eperf = 1:n_splitE_Eperf
            splitE_Eperf_nm = splitE_Eperf{iE_Eperf};
            
            switch DCM_mode
                case 1 % all sessions independent
                    for iRun = 1:n_runs
                        run_nb_nm = ['r',num2str(runs.runsToKeep(iRun))];
                        
                        %% effort period onset
                        n_regs = n_regs + 1;
                        reg_names{n_regs} = ['ONSET effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',run_nb_nm];
                        dispRegFn([num2str(n_regs),') ONSET effort period ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',run_nb_nm,': ',...
                            GLMprm.model_onset.(task_id_nm).Eperf,' '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs = n_regs + add_drv;
                        
                        %% effort period regressors
                        % binary variable indicating when choice = high effort option
                        switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).choiceHighE
                            case 0
                            case {1,2}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',run_nb_nm,': choice = highE'];
                                dispRegFn([num2str(n_regs),') effort period ',run_nb_nm,': choice hE '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % money chosen
                        switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).money_chosen
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',run_nb_nm,': money chosen'];
                                dispRegFn([num2str(n_regs),') effort period ',run_nb_nm,': money chosen (amounts) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',run_nb_nm,': money chosen'];
                                dispRegFn([num2str(n_regs),') effort period ',run_nb_nm,': |money chosen| (amounts) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',run_nb_nm,': money chosen'];
                                dispRegFn([num2str(n_regs),') effort period ',run_nb_nm,': money chosen (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',run_nb_nm,': money chosen'];
                                dispRegFn([num2str(n_regs),') effort period ',run_nb_nm,': |money chosen| (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % effort chosen
                        switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).E_chosen
                            case 0
                            case {1,4}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',run_nb_nm,': effort chosen'];
                                dispRegFn([num2str(n_regs),') effort period ',run_nb_nm,': effort chosen (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                error('case not ready yet');
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',run_nb_nm,': effort-by-choice'];
                                dispRegFn([num2str(n_regs),') effort period ',run_nb_nm,': effort-by-choice '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('case not ready yet');
                        end
                        
                        % force peak
                        switch task_id_nm
                            case 'Ep'
                                % force peak
                                switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).F_peak
                                    case {1,2,3,4}
                                        n_regs = n_regs + 1;
                                        reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',run_nb_nm,': force peak'];
                                        dispRegFn([num2str(n_regs),') effort period ',run_nb_nm,': force peak '],dispRegs);
                                        % if derivative added => add derivatives
                                        n_regs = n_regs + add_drv;
                                end
                                
                                % force integral
                                switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).F_integral
                                    case {1,3,5,7}
                                        n_regs = n_regs + 1;
                                        reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',run_nb_nm,': force integral'];
                                        dispRegFn([num2str(n_regs),') effort period ',run_nb_nm,': force integral '],dispRegs);
                                        % if derivative added => add derivatives
                                        n_regs = n_regs + add_drv;
                                    case {2,4,6,8}
                                        n_regs = n_regs + 1;
                                        reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',run_nb_nm,': force integral overshoot'];
                                        dispRegFn([num2str(n_regs),') effort period ',run_nb_nm,': force integral overshoot '],dispRegs);
                                        % if derivative added => add derivatives
                                        n_regs = n_regs + add_drv;
                                end
                                
                            case 'Em'
                                % efficacy
                                switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).efficacy
                                    case {1,2,3,4,5,6,7,8}
                                        n_regs = n_regs + 1;
                                        reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',run_nb_nm,': efficacy'];
                                        dispRegFn([num2str(n_regs),') effort period ',run_nb_nm,': efficacy '],dispRegs);
                                        % if derivative added => add derivatives
                                        n_regs = n_regs + add_drv;
                                end
                                
                                % average RT for all numbers of each trial
                                switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).RT_avg
                                    case 1
                                        n_regs = n_regs + 1;
                                        reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',run_nb_nm,': average RT effort'];
                                        dispRegFn([num2str(n_regs),') effort period ',run_nb_nm,': average RT effort '],dispRegs);
                                        % if derivative added => add derivatives
                                        n_regs = n_regs + add_drv;
                                end
                                
                                % number of correct answers
                                switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).n_correct
                                    case 1
                                        n_regs = n_regs + 1;
                                        reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',run_nb_nm,': number of correct answers'];
                                        dispRegFn([num2str(n_regs),') effort period ',run_nb_nm,': number of correct answers '],dispRegs);
                                        % if derivative added => add derivatives
                                        n_regs = n_regs + add_drv;
                                end
                                
                                % number of errors
                                switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).n_errors
                                    case 1
                                        n_regs = n_regs + 1;
                                        reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',run_nb_nm,': number of errors'];
                                        dispRegFn([num2str(n_regs),') effort period ',run_nb_nm,': number of errors '],dispRegs);
                                        % if derivative added => add derivatives
                                        n_regs = n_regs + add_drv;
                                end
                        end
                        
                        % net value chosen option
                        switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).NV_chosen
                            case {1,3}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',run_nb_nm,': NVch-NVunch'];
                                dispRegFn([num2str(n_regs),') effort period ',run_nb_nm,': NVch-NVunch '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2 % p(chosen)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',run_nb_nm,': p(chosen)'];
                                dispRegFn([num2str(n_regs),') effort period ',run_nb_nm,': p(chosen) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % net value hE
                        switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).NV_varOption
                            case 1 % NV(high E) - NV(low E)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',run_nb_nm,': net value high E - low E'];
                                dispRegFn([num2str(n_regs),') effort period ',run_nb_nm,': net value high E - low E '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2 % |NV(high E) - NV(low E)|
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',run_nb_nm,': abs(NVhE-NVlE)'];
                                dispRegFn([num2str(n_regs),') effort period ',run_nb_nm,': |net value high E - low E| '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3 % p(choice = high E)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',run_nb_nm,': p(choice = high E)'];
                                dispRegFn([num2str(n_regs),') effort period ',run_nb_nm,': p(choice = high E) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4 % NV(high E) - NV(low E) + bias
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',run_nb_nm,': net value high E - low E + bias'];
                                dispRegFn([num2str(n_regs),') effort period ',run_nb_nm,': net value high E - low E + bias '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 5 % |NV(high E) - NV(low E) + bias|
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',run_nb_nm,': abs(NVhE-NVlE+bias)'];
                                dispRegFn([num2str(n_regs),') effort period ',run_nb_nm,': |net value high E - low E + bias| '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % net value hE bis
                        switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).NV_varOption_bis
                            case 1 % NV(high E) - NV(low E)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',run_nb_nm,': net value high E - low E'];
                                dispRegFn([num2str(n_regs),') effort period ',run_nb_nm,': net value high E - low E '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2 % |NV(high E) - NV(low E)|
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',run_nb_nm,': abs(NVhE-NVlE)'];
                                dispRegFn([num2str(n_regs),') effort period ',run_nb_nm,': |net value high E - low E| '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3 % p(choice = high E)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',run_nb_nm,': p(choice = high E)'];
                                dispRegFn([num2str(n_regs),') effort period ',run_nb_nm,': p(choice = high E) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4 % NV(high E) - NV(low E) + bias
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',run_nb_nm,': net value high E - low E + bias'];
                                dispRegFn([num2str(n_regs),') effort period ',run_nb_nm,': net value high E - low E + bias '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 5 % |NV(high E) - NV(low E) + bias|
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',run_nb_nm,': abs(NVhE-NVlE+bias)'];
                                dispRegFn([num2str(n_regs),') effort period ',run_nb_nm,': |net value high E - low E + bias| '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % RT first answer
                        switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).RT_1stAnswer
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',run_nb_nm,': effort latency'];
                                dispRegFn([num2str(n_regs),') effort period ',run_nb_nm,': effort latency (raw) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % fatigue
                        switch task_id_nm
                            case 'Ep'
                                % fatigue
                                switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).fatigue
                                    case {1,2}
                                        n_regs = n_regs + 1;
                                        reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',run_nb_nm,': fatigue'];
                                        dispRegFn([num2str(n_regs),') effort period ',run_nb_nm,': fatigue '],dispRegs);
                                        % if derivative added => add derivatives
                                        n_regs = n_regs + add_drv;
                                end
                                
                                % Ech*fatigue
                                switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).Ech_x_fatigue
                                    case 1
                                        n_regs = n_regs + 1;
                                        reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',run_nb_nm,': Ech_x_fatigue'];
                                        dispRegFn([num2str(n_regs),') effort period ',run_nb_nm,': Ech_x_fatigue '],dispRegs);
                                        % if derivative added => add derivatives
                                        n_regs = n_regs + add_drv;
                                end
                                
                            case 'Em'
                                % efficacy during the previous trial (learning boost)
                                switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).prevEfficacy
                                    case {1,2,3,4,5,6,7,8}
                                        n_regs = n_regs + 1;
                                        reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',run_nb_nm,': previous efficacy'];
                                        dispRegFn([num2str(n_regs),') effort period ',run_nb_nm,': previous efficacy '],dispRegs);
                                        % if derivative added => add derivatives
                                        n_regs = n_regs + add_drv;
                                end
                                
                                % (effort chosen)*(efficacy during the previous trial)
                                switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).Ech_x_prevEfficacy
                                    case {1,2,3,4}
                                        n_regs = n_regs + 1;
                                        reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',run_nb_nm,': (Ech)x(previous efficacy)'];
                                        dispRegFn([num2str(n_regs),') effort period ',run_nb_nm,': (Ech)x(previous efficacy) '],dispRegs);
                                        % if derivative added => add derivatives
                                        n_regs = n_regs + add_drv;
                                end
                        end
                        
                        % trial number
                        switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).trialN
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',run_nb_nm,': trial number'];
                                dispRegFn([num2str(n_regs),') effort period ',run_nb_nm,': trial number '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',run_nb_nm,': (trial number)x(Echosen-Eunchosen) (E levels)'];
                                dispRegFn([num2str(n_regs),') effort period ',run_nb_nm,': (trial number)x(effort chosen-effort unchosen) (effort levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',run_nb_nm,': (trial number)x(EnonDef-Edef) (E levels)'];
                                dispRegFn([num2str(n_regs),') effort period ',run_nb_nm,': (trial number)x(effort non-default - effort default) (effort levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',run_nb_nm,': (trial number)x(EnonDef) (E levels)'];
                                dispRegFn([num2str(n_regs),') effort period ',run_nb_nm,': (trial number)x(effort non-default) (effort levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % confidence
                        switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).confidence
                            case 1 % confidence ratings 0/1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',run_nb_nm,': confidence'];
                                dispRegFn([num2str(n_regs),') effort period ',run_nb_nm,': confidence (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case {2,3,4} % confidence inferred by the model
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',run_nb_nm,': confidence'];
                                dispRegFn([num2str(n_regs),') effort period ',run_nb_nm,': confidence (inferred by the model) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                    end % run loop
                case {2,4,5} % all tasks independent but sessions pooled
                    for iTask = 1:nTasks
                        task_nm = tasks{iTask};
                        
                        %% effort period onset
                        n_regs = n_regs + 1;
                        reg_names{n_regs} = ['ONSET effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',task_nm];
                        dispRegFn([num2str(n_regs),') ONSET effort period ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',task_nm,': ',...
                            GLMprm.model_onset.(task_id_nm).Eperf,' '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs = n_regs + add_drv;
                        
                        %% effort period regressors
                        % binary variable indicating when choice = high effort option
                        switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).choiceHighE
                            case 0
                            case {1,2}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',task_nm,': choice = highE'];
                                dispRegFn([num2str(n_regs),') effort period ',task_nm,': choice hE '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % money chosen
                        switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).money_chosen
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',task_nm,': money chosen'];
                                dispRegFn([num2str(n_regs),') effort period ',task_nm,': money chosen (amounts) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',task_nm,': money chosen'];
                                dispRegFn([num2str(n_regs),') effort period ',task_nm,': |money chosen| (amounts) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',task_nm,': money chosen'];
                                dispRegFn([num2str(n_regs),') effort period ',task_nm,': money chosen (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',task_nm,': money chosen'];
                                dispRegFn([num2str(n_regs),') effort period ',task_nm,': |money chosen| (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % effort chosen
                        switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).E_chosen
                            case 0
                            case {1,4}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',task_nm,': effort chosen'];
                                dispRegFn([num2str(n_regs),') effort period ',task_nm,': effort chosen (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                error('case not ready yet');
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',task_nm,': effort-by-choice'];
                                dispRegFn([num2str(n_regs),') effort period ',task_nm,': effort-by-choice '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('case not ready yet');
                        end
                        
                        % force peak
                        switch task_id_nm
                            case 'Ep'
                                % force peak
                                switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).F_peak
                                    case {1,2,3,4}
                                        n_regs = n_regs + 1;
                                        reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',task_nm,': force peak'];
                                        dispRegFn([num2str(n_regs),') effort period ',task_nm,': force peak '],dispRegs);
                                        % if derivative added => add derivatives
                                        n_regs = n_regs + add_drv;
                                end
                                
                                % force integral
                                switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).F_integral
                                    case {1,3,5,7}
                                        n_regs = n_regs + 1;
                                        reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',task_nm,': force integral'];
                                        dispRegFn([num2str(n_regs),') effort period ',task_nm,': force integral '],dispRegs);
                                        % if derivative added => add derivatives
                                        n_regs = n_regs + add_drv;
                                    case {2,4,6,8}
                                        n_regs = n_regs + 1;
                                        reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',task_nm,': force integral overshoot'];
                                        dispRegFn([num2str(n_regs),') effort period ',task_nm,': force integral overshoot '],dispRegs);
                                        % if derivative added => add derivatives
                                        n_regs = n_regs + add_drv;
                                end
                                
                            case 'Em'
                                % efficacy
                                switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).efficacy
                                    case {1,2,3,4,5,6,7,8}
                                        n_regs = n_regs + 1;
                                        reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',task_nm,': efficacy'];
                                        dispRegFn([num2str(n_regs),') effort period ',task_nm,': efficacy '],dispRegs);
                                        % if derivative added => add derivatives
                                        n_regs = n_regs + add_drv;
                                end
                                
                                % average RT for all numbers of each trial
                                switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).RT_avg
                                    case 1
                                        n_regs = n_regs + 1;
                                        reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',task_nm,': average RT effort'];
                                        dispRegFn([num2str(n_regs),') effort period ',task_nm,': average RT effort '],dispRegs);
                                        % if derivative added => add derivatives
                                        n_regs = n_regs + add_drv;
                                end
                                
                                % number of correct answers
                                switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).n_correct
                                    case 1
                                        n_regs = n_regs + 1;
                                        reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',task_nm,': number of correct answers'];
                                        dispRegFn([num2str(n_regs),') effort period ',task_nm,': number of correct answers '],dispRegs);
                                        % if derivative added => add derivatives
                                        n_regs = n_regs + add_drv;
                                end
                                
                                % number of errors
                                switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).n_errors
                                    case 1
                                        n_regs = n_regs + 1;
                                        reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',task_nm,': number of errors'];
                                        dispRegFn([num2str(n_regs),') effort period ',task_nm,': number of errors '],dispRegs);
                                        % if derivative added => add derivatives
                                        n_regs = n_regs + add_drv;
                                end
                        end
                        
                        % net value chosen option
                        switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).NV_chosen
                            case {1,3}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',task_nm,': NVch-NVunch'];
                                dispRegFn([num2str(n_regs),') effort period ',task_nm,': NVch-NVunch '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2 % p(chosen)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',task_nm,': p(chosen)'];
                                dispRegFn([num2str(n_regs),') effort period ',task_nm,': p(chosen) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % net value hE
                        switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).NV_varOption
                            case 1 % NV(high E) - NV(low E)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',task_nm,': net value high E - low E'];
                                dispRegFn([num2str(n_regs),') effort period ',task_nm,': net value high E - low E '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2 % |NV(high E) - NV(low E)|
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',task_nm,': abs(NVhE-NVlE)'];
                                dispRegFn([num2str(n_regs),') effort period ',task_nm,': |net value high E - low E| '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3 % p(choice = high E)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',task_nm,': p(choice = high E)'];
                                dispRegFn([num2str(n_regs),') effort period ',task_nm,': p(choice = high E) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4 % NV(high E) - NV(low E) + bias
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',task_nm,': net value high E - low E + bias'];
                                dispRegFn([num2str(n_regs),') effort period ',task_nm,': net value high E - low E + bias '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 5 % |NV(high E) - NV(low E) + bias|
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',task_nm,': abs(NVhE-NVlE+bias)'];
                                dispRegFn([num2str(n_regs),') effort period ',task_nm,': |net value high E - low E + bias| '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % net value hE bis
                        switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).NV_varOption_bis
                            case 1 % NV(high E) - NV(low E)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',task_nm,': net value high E - low E'];
                                dispRegFn([num2str(n_regs),') effort period ',task_nm,': net value high E - low E '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2 % |NV(high E) - NV(low E)|
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',task_nm,': abs(NVhE-NVlE)'];
                                dispRegFn([num2str(n_regs),') effort period ',task_nm,': |net value high E - low E| '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3 % p(choice = high E)
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',task_nm,': p(choice = high E)'];
                                dispRegFn([num2str(n_regs),') effort period ',task_nm,': p(choice = high E) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4 % NV(high E) - NV(low E) + bias
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',task_nm,': net value high E - low E + bias'];
                                dispRegFn([num2str(n_regs),') effort period ',task_nm,': net value high E - low E + bias '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 5 % |NV(high E) - NV(low E) + bias|
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',task_nm,': abs(NVhE-NVlE+bias)'];
                                dispRegFn([num2str(n_regs),') effort period ',task_nm,': |net value high E - low E + bias| '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % RT first answer
                        switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).RT_1stAnswer
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',task_nm,': effort latency'];
                                dispRegFn([num2str(n_regs),') effort period ',task_nm,': effort latency (raw) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % fatigue
                        switch task_id_nm
                            case 'Ep'
                                % fatigue
                                switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).fatigue
                                    case {1,2}
                                        n_regs = n_regs + 1;
                                        reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',task_nm,': fatigue'];
                                        dispRegFn([num2str(n_regs),') effort period ',task_nm,': fatigue '],dispRegs);
                                        % if derivative added => add derivatives
                                        n_regs = n_regs + add_drv;
                                end
                                
                                % Ech*fatigue
                                switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).Ech_x_fatigue
                                    case 1
                                        n_regs = n_regs + 1;
                                        reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',task_nm,': Ech_x_fatigue'];
                                        dispRegFn([num2str(n_regs),') effort period ',task_nm,': Ech_x_fatigue '],dispRegs);
                                        % if derivative added => add derivatives
                                        n_regs = n_regs + add_drv;
                                end
                                
                            case 'Em'
                                % efficacy during the previous trial (learning boost)
                                switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).prevEfficacy
                                    case {1,2,3,4,5,6,7,8}
                                        n_regs = n_regs + 1;
                                        reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',task_nm,': previous efficacy'];
                                        dispRegFn([num2str(n_regs),') effort period ',task_nm,': previous efficacy '],dispRegs);
                                        % if derivative added => add derivatives
                                        n_regs = n_regs + add_drv;
                                end
                                
                                % (effort chosen)*(efficacy during the previous trial)
                                switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).Ech_x_prevEfficacy
                                    case {1,2,3,4}
                                        n_regs = n_regs + 1;
                                        reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',task_nm,': (Ech)x(previous efficacy)'];
                                        dispRegFn([num2str(n_regs),') effort period ',task_nm,': (Ech)x(previous efficacy) '],dispRegs);
                                        % if derivative added => add derivatives
                                        n_regs = n_regs + add_drv;
                                end
                        end
                        
                        % trial number
                        switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).trialN
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',task_nm,': trial number'];
                                dispRegFn([num2str(n_regs),') effort period ',task_nm,': trial number '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',task_nm,': (trial number)x(Echosen-Eunchosen) (E levels)'];
                                dispRegFn([num2str(n_regs),') effort period ',task_nm,': (trial number)x(effort chosen-effort unchosen) (effort levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',task_nm,': (trial number)x(EnonDef-Edef) (E levels)'];
                                dispRegFn([num2str(n_regs),') effort period ',task_nm,': (trial number)x(effort non-default - effort default) (effort levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 4
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',task_nm,': (trial number)x(EnonDef) (E levels)'];
                                dispRegFn([num2str(n_regs),') effort period ',task_nm,': (trial number)x(effort non-default) (effort levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % confidence
                        switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).confidence
                            case 1 % confidence ratings 0/1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',task_nm,': confidence'];
                                dispRegFn([num2str(n_regs),') effort period ',task_nm,': confidence (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case {2,3,4} % confidence inferred by the model
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,' ',task_nm,': confidence'];
                                dispRegFn([num2str(n_regs),') effort period ',task_nm,': confidence (inferred by the model) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                    end % task loop
                case 3 % all pooled across sessions
                    
                    %% effort period onset
                    n_regs = n_regs + 1;
                    reg_names{n_regs} = ['ONSET effort ',RP_Eperf_nm,' ',splitE_Eperf_nm];
                    dispRegFn([num2str(n_regs),') ONSET effort period ',RP_Eperf_nm,' ',splitE_Eperf_nm,': ',...
                        GLMprm.model_onset.(task_id_nm).Eperf,' '],dispRegs);
                    % if derivative added => add derivatives
                    n_regs = n_regs + add_drv;
                    
                    %% effort period regressors
                    % binary variable indicating when choice = high effort option
                    switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).choiceHighE
                        case 0
                        case {1,2}
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': choice = highE'];
                            dispRegFn([num2str(n_regs),') effort period: choice hE '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % money chosen
                    switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).money_chosen
                        case 1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': money chosen'];
                            dispRegFn([num2str(n_regs),') effort period: money chosen (amounts) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': money chosen'];
                            dispRegFn([num2str(n_regs),') effort period: |money chosen| (amounts) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 3
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': money chosen'];
                            dispRegFn([num2str(n_regs),') effort period: money chosen (levels) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 4
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': money chosen'];
                            dispRegFn([num2str(n_regs),') effort period: |money chosen| (levels) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                    end
                    
                    % effort chosen
                    switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).E_chosen
                        case 0
                        case {1,4}
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': effort chosen'];
                            dispRegFn([num2str(n_regs),') effort period: effort chosen (levels) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2
                            error('case not ready yet');
                        case 3
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': effort-by-choice'];
                            dispRegFn([num2str(n_regs),') effort period: effort-by-choice '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('case not ready yet');
                    end
                    
                    % force peak
                    switch task_id_nm
                        case 'Ep'
                            % force peak
                            switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).F_peak
                                case {1,2,3,4}
                                    n_regs = n_regs + 1;
                                    reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': force peak'];
                                    dispRegFn([num2str(n_regs),') effort period: force peak '],dispRegs);
                                    % if derivative added => add derivatives
                                    n_regs = n_regs + add_drv;
                            end
                            
                            % force integral
                            switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).F_integral
                                case {1,3,5,7}
                                    n_regs = n_regs + 1;
                                    reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': force integral'];
                                    dispRegFn([num2str(n_regs),') effort period: force integral '],dispRegs);
                                    % if derivative added => add derivatives
                                    n_regs = n_regs + add_drv;
                                case {2,4,6,8}
                                    n_regs = n_regs + 1;
                                    reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': force integral overshoot'];
                                    dispRegFn([num2str(n_regs),') effort period: force integral overshoot '],dispRegs);
                                    % if derivative added => add derivatives
                                    n_regs = n_regs + add_drv;
                            end
                            
                        case 'Em'
                            % efficacy
                            switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).efficacy
                                case {1,2,3,4,5,6,7,8}
                                    n_regs = n_regs + 1;
                                    reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': efficacy'];
                                    dispRegFn([num2str(n_regs),') effort period: efficacy '],dispRegs);
                                    % if derivative added => add derivatives
                                    n_regs = n_regs + add_drv;
                            end
                            
                            % average RT for all numbers of each trial
                            switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).RT_avg
                                case 1
                                    n_regs = n_regs + 1;
                                    reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': average RT effort'];
                                    dispRegFn([num2str(n_regs),') effort period: average RT effort '],dispRegs);
                                    % if derivative added => add derivatives
                                    n_regs = n_regs + add_drv;
                            end
                            
                            % number of correct answers
                            switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).n_correct
                                case 1
                                    n_regs = n_regs + 1;
                                    reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': number of correct answers'];
                                    dispRegFn([num2str(n_regs),') effort period: number of correct answers '],dispRegs);
                                    % if derivative added => add derivatives
                                    n_regs = n_regs + add_drv;
                            end
                            
                            % number of errors
                            switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).n_errors
                                case 1
                                    n_regs = n_regs + 1;
                                    reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': number of errors'];
                                    dispRegFn([num2str(n_regs),') effort period: number of errors '],dispRegs);
                                    % if derivative added => add derivatives
                                    n_regs = n_regs + add_drv;
                            end
                    end
                    
                    % net value chosen option
                    switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).NV_chosen
                        case {1,3}
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': NVch-NVunch'];
                            dispRegFn([num2str(n_regs),') effort period: NVch-NVunch '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2 % p(chosen)
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': p(chosen)'];
                            dispRegFn([num2str(n_regs),') effort period: p(chosen) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                    end
                    
                    % net value hE
                    switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).NV_varOption
                        case 1 % NV(high E) - NV(low E)
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': net value high E - low E'];
                            dispRegFn([num2str(n_regs),') effort period: net value high E - low E '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2 % |NV(high E) - NV(low E)|
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': abs(NVhE-NVlE)'];
                            dispRegFn([num2str(n_regs),') effort period: |net value high E - low E| '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 3 % p(choice = high E)
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': p(choice = high E)'];
                            dispRegFn([num2str(n_regs),') effort period: p(choice = high E) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 4 % NV(high E) - NV(low E) + bias
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': net value high E - low E + bias'];
                            dispRegFn([num2str(n_regs),') effort period: net value high E - low E + bias '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 5 % |NV(high E) - NV(low E) + bias|
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': abs(NVhE-NVlE+bias)'];
                            dispRegFn([num2str(n_regs),') effort period: |net value high E - low E + bias| '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                    end
                    
                    % net value hE bis
                    switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).NV_varOption_bis
                        case 1 % NV(high E) - NV(low E)
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': net value high E - low E'];
                            dispRegFn([num2str(n_regs),') effort period: net value high E - low E '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2 % |NV(high E) - NV(low E)|
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': abs(NVhE-NVlE)'];
                            dispRegFn([num2str(n_regs),') effort period: |net value high E - low E| '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 3 % p(choice = high E)
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': p(choice = high E)'];
                            dispRegFn([num2str(n_regs),') effort period: p(choice = high E) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 4 % NV(high E) - NV(low E) + bias
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': net value high E - low E + bias'];
                            dispRegFn([num2str(n_regs),') effort period: net value high E - low E + bias '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 5 % |NV(high E) - NV(low E) + bias|
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': abs(NVhE-NVlE+bias)'];
                            dispRegFn([num2str(n_regs),') effort period: |net value high E - low E + bias| '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                    end
                    
                    % RT first answer
                    switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).RT_1stAnswer
                        case 1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': effort latency'];
                            dispRegFn([num2str(n_regs),') effort period: effort latency (raw) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                    end
                    
                    % fatigue
                    switch task_id_nm
                        case 'Ep'
                            % fatigue
                            switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).fatigue
                                case {1,2}
                                    n_regs = n_regs + 1;
                                    reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': fatigue'];
                                    dispRegFn([num2str(n_regs),') effort period: fatigue '],dispRegs);
                                    % if derivative added => add derivatives
                                    n_regs = n_regs + add_drv;
                            end
                            
                            % Ech*fatigue
                            switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).Ech_x_fatigue
                                case 1
                                    n_regs = n_regs + 1;
                                    reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': Ech_x_fatigue'];
                                    dispRegFn([num2str(n_regs),') effort period: Ech_x_fatigue '],dispRegs);
                                    % if derivative added => add derivatives
                                    n_regs = n_regs + add_drv;
                            end
                            
                        case 'Em'
                            % efficacy during the previous trial (learning boost)
                            switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).prevEfficacy
                                case {1,2,3,4,5,6,7,8}
                                    n_regs = n_regs + 1;
                                    reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': previous efficacy'];
                                    dispRegFn([num2str(n_regs),') effort period: previous efficacy '],dispRegs);
                                    % if derivative added => add derivatives
                                    n_regs = n_regs + add_drv;
                            end
                            
                            % (effort chosen)*(efficacy during the previous trial)
                            switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).Ech_x_prevEfficacy
                                case {1,2,3,4}
                                    n_regs = n_regs + 1;
                                    reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': (Ech)x(previous efficacy)'];
                                    dispRegFn([num2str(n_regs),') effort period: (Ech)x(previous efficacy) '],dispRegs);
                                    % if derivative added => add derivatives
                                    n_regs = n_regs + add_drv;
                            end
                    end
                    
                    % trial number
                    switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).trialN
                        case 1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': trial number'];
                            dispRegFn([num2str(n_regs),') effort period: trial number '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': (trial number)x(Echosen-Eunchosen) (E levels)'];
                            dispRegFn([num2str(n_regs),') effort period: (trial number)x(effort chosen-effort unchosen) (effort levels) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 3
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': (trial number)x(EnonDef-Edef) (E levels)'];
                            dispRegFn([num2str(n_regs),') effort period: (trial number)x(effort non-default - effort default) (effort levels) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 4
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': (trial number)x(EnonDef) (E levels)'];
                            dispRegFn([num2str(n_regs),') effort period: (trial number)x(effort non-default) (effort levels) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                    end
                    
                    % confidence
                    switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).confidence
                        case 1 % confidence ratings 0/1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': confidence'];
                            dispRegFn([num2str(n_regs),') effort period: confidence (levels) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case {2,3,4} % confidence inferred by the model
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': confidence'];
                            dispRegFn([num2str(n_regs),') effort period: confidence (inferred by the model) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                    end
                    
            end % DCM mode
        end % loop effort levels
    end % loop reward/punishment
end % choice onset

%% feedback period
if ~strcmp(GLMprm.model_onset.(task_id_nm).fbk,'none')
    % check if trials are split or not
    if GLMprm.fbk.(task_id_nm).RPpool == 1 % pool reward and punishment trials
        RP_fbk = {'RP'};
    elseif GLMprm.fbk.(task_id_nm).RPpool == 0 % split reward and punishment
        RP_fbk = {'R','P'};
    end % RP pool
    n_RP_fbk = length(RP_fbk);
    
    % check if trials are split or not according to Effort levels
    switch GLMprm.fbk.(task_id_nm).splitPerE
        case 0 % pool trials
            splitE_fbk = {'E'};
        case 1 % split according to effort proposed
            splitE_fbk = {'E1','E2','E3'};
        case 2 % split according to effort chosen
            splitE_fbk = {'Ech0','Ech1','Ech2','Ech3'};
        case 3 % split according to option chosen (low/high effort)
            splitE_fbk = {'lEch','hEch'};
    end % Effort level pool
    n_splitE_fbk = length(splitE_fbk);
    
    % loop through conditions for feedback period
    for iRP_fbk = 1:n_RP_fbk
        RP_fbk_nm = RP_fbk{iRP_fbk};
        
        for iE_fbk = 1:n_splitE_fbk
            splitE_fbk_nm = splitE_fbk{iE_fbk};
            
            switch DCM_mode
                case 1 % all sessions independent
                    for iRun = 1:n_runs
                        run_nb_nm = ['r',num2str(runs.runsToKeep(iRun))];
                        
                        %% feedback onset
                        n_regs = n_regs + 1;
                        reg_names{n_regs} = ['ONSET feedback ',RP_fbk_nm,' ',splitE_fbk_nm,' ',run_nb_nm];
                        dispRegFn([num2str(n_regs),') ONSET feedback ',run_nb_nm,': ',GLMprm.model_onset.(task_id_nm).fbk,' '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs = n_regs + add_drv;
                        
                        %% feedback regressors
                        % win vs loss
                        switch GLMprm.fbk.(task_id_nm).(RP_fbk_nm).(splitE_fbk_nm).win_vs_loss
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG feedback ',RP_fbk_nm,' ',splitE_fbk_nm,' ',run_nb_nm,': win-loss'];
                                dispRegFn([num2str(n_regs),') feedback ',run_nb_nm,': win-loss '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % binary variable indicating when choice = high effort option
                        switch GLMprm.fbk.(task_id_nm).(RP_fbk_nm).(splitE_fbk_nm).choiceHighE
                            case 0
                            case {1,2}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG feedback ',RP_fbk_nm,' ',splitE_fbk_nm,' ',run_nb_nm,': choice = highE'];
                                dispRegFn([num2str(n_regs),') feedback ',run_nb_nm,': choice hE '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % money obtained
                        switch GLMprm.fbk.(task_id_nm).(RP_fbk_nm).(splitE_fbk_nm).money_obtained
                            case 1 % money amount
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG feedback ',RP_fbk_nm,' ',splitE_fbk_nm,' ',run_nb_nm,': money obtained'];
                                dispRegFn([num2str(n_regs),') feedback ',run_nb_nm,': money obtained (amount) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2 % |money amount|
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG feedback ',RP_fbk_nm,' ',splitE_fbk_nm,' ',run_nb_nm,': abs(money obtained)'];
                                dispRegFn([num2str(n_regs),') feedback ',run_nb_nm,': |money obtained| (amount) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % effort performed
                        switch GLMprm.fbk.(task_id_nm).(RP_fbk_nm).(splitE_fbk_nm).E_made
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG feedback ',RP_fbk_nm,' ',splitE_fbk_nm,' ',run_nb_nm,': effort performed'];
                                dispRegFn([num2str(n_regs),') feedback ',run_nb_nm,': effort performed (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % trial number
                        switch GLMprm.fbk.(task_id_nm).(RP_fbk_nm).(splitE_fbk_nm).trialN
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG feedback ',RP_fbk_nm,' ',splitE_fbk_nm,' ',run_nb_nm,': trial number'];
                                dispRegFn([num2str(n_regs),') feedback ',run_nb_nm,': trial number '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG feedback ',RP_fbk_nm,' ',splitE_fbk_nm,' ',run_nb_nm,': trial number x (Echosen-Eunchosen) (E levels)'];
                                dispRegFn([num2str(n_regs),') feedback ',run_nb_nm,': trial number x (effort chosen-effort unchosen) (effort levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG feedback ',RP_fbk_nm,' ',splitE_fbk_nm,' ',run_nb_nm,': trial number x (EnonDef-Edef) (E levels)'];
                                dispRegFn([num2str(n_regs),') feedback ',run_nb_nm,': trial number x (effort non-default - effort default) (effort levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % confidence
                        switch GLMprm.fbk.(task_id_nm).(RP_fbk_nm).(splitE_fbk_nm).confidence
                            case 1 % confidence rated by the subjects
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG feedback ',RP_fbk_nm,' ',splitE_fbk_nm,' ',run_nb_nm,': confidence'];
                                dispRegFn([num2str(n_regs),') feedback ',run_nb_nm,': confidence (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case {2,3,4} % confidence inferred by the model
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG feedback ',RP_fbk_nm,' ',splitE_fbk_nm,' ',run_nb_nm,': confidence'];
                                dispRegFn([num2str(n_regs),') feedback ',run_nb_nm,': confidence (inferred by the model) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                    end % run loop
                case 2 % all tasks independent but sessions pooled
                    for iTask = 1:nTasks
                        task_nm = tasks{iTask};
                        
                        %% feedback onset
                        n_regs = n_regs + 1;
                        reg_names{n_regs} = ['ONSET feedback ',RP_fbk_nm,' ',splitE_fbk_nm,' ',task_nm];
                        dispRegFn([num2str(n_regs),') ONSET feedback ',task_nm,': ',GLMprm.model_onset.(task_id_nm).fbk,' '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs = n_regs + add_drv;
                        
                        %% feedback regressors
                        % win vs loss
                        switch GLMprm.fbk.(task_id_nm).(RP_fbk_nm).(splitE_fbk_nm).win_vs_loss
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG feedback ',RP_fbk_nm,' ',splitE_fbk_nm,' ',task_nm,': win-loss'];
                                dispRegFn([num2str(n_regs),') feedback ',task_nm,': win-loss '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % binary variable indicating when choice = high effort option
                        switch GLMprm.fbk.(task_id_nm).(RP_fbk_nm).(splitE_fbk_nm).choiceHighE
                            case 0
                            case {1,2}
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG feedback ',RP_fbk_nm,' ',splitE_fbk_nm,' ',task_nm,': choice = highE'];
                                dispRegFn([num2str(n_regs),') feedback ',task_nm,': choice hE '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            otherwise
                                error('not ready yet');
                        end
                        
                        % money obtained
                        switch GLMprm.fbk.(task_id_nm).(RP_fbk_nm).(splitE_fbk_nm).money_obtained
                            case 1 % money amount
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG feedback ',RP_fbk_nm,' ',splitE_fbk_nm,' ',task_nm,': money obtained'];
                                dispRegFn([num2str(n_regs),') feedback ',task_nm,': money obtained (amount) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2 % |money amount|
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG feedback ',RP_fbk_nm,' ',splitE_fbk_nm,' ',task_nm,': abs(money obtained)'];
                                dispRegFn([num2str(n_regs),') feedback ',task_nm,': |money obtained| (amount) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % effort performed
                        switch GLMprm.fbk.(task_id_nm).(RP_fbk_nm).(splitE_fbk_nm).E_made
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG feedback ',RP_fbk_nm,' ',splitE_fbk_nm,' ',task_nm,': effort performed'];
                                dispRegFn([num2str(n_regs),') feedback ',task_nm,': effort performed (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % trial number
                        switch GLMprm.fbk.(task_id_nm).(RP_fbk_nm).(splitE_fbk_nm).trialN
                            case 1
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG feedback ',RP_fbk_nm,' ',splitE_fbk_nm,' ',task_nm,': trial number'];
                                dispRegFn([num2str(n_regs),') feedback ',task_nm,': trial number '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 2
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG feedback ',RP_fbk_nm,' ',splitE_fbk_nm,' ',task_nm,': trial number x (Echosen-Eunchosen) (E levels)'];
                                dispRegFn([num2str(n_regs),') feedback ',task_nm,': trial number x (effort chosen-effort unchosen) (effort levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case 3
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG feedback ',RP_fbk_nm,' ',splitE_fbk_nm,' ',task_nm,': trial number x (EnonDef-Edef) (E levels)'];
                                dispRegFn([num2str(n_regs),') feedback ',task_nm,': trial number x (effort non-default - effort default) (effort levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                        
                        % confidence
                        switch GLMprm.fbk.(task_id_nm).(RP_fbk_nm).(splitE_fbk_nm).confidence
                            case 1 % confidence rated by the subjects
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG feedback ',RP_fbk_nm,' ',splitE_fbk_nm,' ',task_nm,': confidence'];
                                dispRegFn([num2str(n_regs),') feedback ',task_nm,': confidence (levels) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                            case {2,3,4} % confidence inferred by the model
                                n_regs = n_regs + 1;
                                reg_names{n_regs} = ['REG feedback ',RP_fbk_nm,' ',splitE_fbk_nm,' ',task_nm,': confidence'];
                                dispRegFn([num2str(n_regs),') feedback ',task_nm,': confidence (inferred by the model) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs = n_regs + add_drv;
                        end
                    end % task loop
                case {3,4,5} % all pooled across sessions
                    
                    %% feedback onset
                    n_regs = n_regs + 1;
                    reg_names{n_regs} = ['ONSET feedback ',RP_fbk_nm,' ',splitE_fbk_nm];
                    dispRegFn([num2str(n_regs),') ONSET feedback: ',GLMprm.model_onset.(task_id_nm).fbk,' '],dispRegs);
                    % if derivative added => add derivatives
                    n_regs = n_regs + add_drv;
                    
                    %% feedback regressors
                    % win vs loss
                    switch GLMprm.fbk.(task_id_nm).(RP_fbk_nm).(splitE_fbk_nm).win_vs_loss
                        case 1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG feedback ',RP_fbk_nm,' ',splitE_fbk_nm,': win-loss'];
                            dispRegFn([num2str(n_regs),') feedback: win-loss '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                    end
                    
                    % binary variable indicating when choice = high effort option
                    switch GLMprm.fbk.(task_id_nm).(RP_fbk_nm).(splitE_fbk_nm).choiceHighE
                        case 0
                        case {1,2}
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG feedback ',RP_fbk_nm,' ',splitE_fbk_nm,': choice = highE'];
                            dispRegFn([num2str(n_regs),') feedback: choice hE '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % money obtained
                    switch GLMprm.fbk.(task_id_nm).(RP_fbk_nm).(splitE_fbk_nm).money_obtained
                        case 1 % money amount
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG feedback ',RP_fbk_nm,' ',splitE_fbk_nm,': money obtained'];
                            dispRegFn([num2str(n_regs),') feedback: money obtained (amount) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2 % |money amount|
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG feedback ',RP_fbk_nm,' ',splitE_fbk_nm,': abs(money obtained)'];
                            dispRegFn([num2str(n_regs),') feedback: |money obtained| (amount) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                    end
                    
                    % effort performed
                    switch GLMprm.fbk.(task_id_nm).(RP_fbk_nm).(splitE_fbk_nm).E_made
                        case 1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG feedback ',RP_fbk_nm,' ',splitE_fbk_nm,': effort performed'];
                            dispRegFn([num2str(n_regs),') feedback: effort performed (levels) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                    end
                    
                    % trial number
                    switch GLMprm.fbk.(task_id_nm).(RP_fbk_nm).(splitE_fbk_nm).trialN
                        case 1
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG feedback ',RP_fbk_nm,' ',splitE_fbk_nm,': trial number'];
                            dispRegFn([num2str(n_regs),') feedback: trial number '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 2
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG feedback ',RP_fbk_nm,' ',splitE_fbk_nm,': trial number x (Echosen-Eunchosen) (E levels)'];
                            dispRegFn([num2str(n_regs),') feedback: trial number x (effort chosen-effort unchosen) (effort levels) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case 3
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG feedback ',RP_fbk_nm,' ',splitE_fbk_nm,': trial number x (EnonDef-Edef) (E levels)'];
                            dispRegFn([num2str(n_regs),') feedback: trial number x (effort non-default - effort default) (effort levels) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                    end
                    
                    % confidence
                    switch GLMprm.fbk.(task_id_nm).(RP_fbk_nm).(splitE_fbk_nm).confidence
                        case 1 % confidence rated by the subjects
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG feedback ',RP_fbk_nm,' ',splitE_fbk_nm,': confidence'];
                            dispRegFn([num2str(n_regs),') feedback: confidence (levels) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                        case {2,3,4} % confidence inferred by the model
                            n_regs = n_regs + 1;
                            reg_names{n_regs} = ['REG feedback ',RP_fbk_nm,' ',splitE_fbk_nm,': confidence'];
                            dispRegFn([num2str(n_regs),') feedback: confidence (inferred by the model) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs = n_regs + add_drv;
                    end
            end % DCM mode
        end % loop effort levels
    end % loop reward/punishment
end % choice onset

%% movement regressors
n_mvmt = 6;
for iMvmt = 1:n_mvmt
    n_regs = n_regs + 1;
    reg_names{n_regs} = 'movement';
end
dispRegFn([num2str(n_regs-n_mvmt+1),'-',num2str(n_regs),') movement regressors '], dispRegs);
% note: if temporal derivative has been added, movement regressors
% are not doubled (ie number stays stable independent of add_drv value)

%% run constants
for iRun = 1:n_runs
    run_nb_nm = ['r',num2str(runs.runsToKeep(iRun))];
    n_regs = n_regs + 1;
    reg_names{n_regs} = [run_nb_nm,' constant'];
end % run loop
dispRegFn([num2str(n_regs-n_runs+1),'-',num2str(n_regs),') run constant regressors '], dispRegs);

end % function

function[] = dispRegFn(text, dispRegs)
switch dispRegs
    case 1
        disp(text);
end
end