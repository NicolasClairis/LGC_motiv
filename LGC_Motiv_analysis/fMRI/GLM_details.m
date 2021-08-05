function [reg_names, n_regs] = GLM_details(GLM)
% [Ep_regs, Em_regs] = GLM_details(GLM)
%GLM_details will provide the number of regressors, the text with the
%details of the GLM, the order of the regressors in the current GLM for
%each task.
%
% INPUTS
% GLM: number of the GLM
%
% OUTPUTS
% reg_names: structure with physical (Ep) and mental (Em) regressors for
% each run
%
% n_regs: structure with number of regressors for each physical (Ep) and
% each mental (Em) run

%% load GLM parameters
GLMprm = which_GLM(GLM);
% disp GLM number
disp(['GLM ',num2str(GLM)]);

%% initialize the variables of interest
[reg_names.Ep, reg_names.Em] = deal( {} );
n_regs.Ep = 0;
n_regs.Em = 0;

%% main parameters
generalPrm = GLMprm.gal;
disp('** general parameters of the GLM **');

% grey matter mask
grey_mask = generalPrm.grey_mask;
switch grey_mask
    case 0
        disp('all voxels included (no grey matter filter)');
    case 1
        disp('Use of 1st level probability grey mask for each subject.');
    case 2
        disp('Use of 1st level probability grey mask for each subject.');
end

% temporal/spatial derivative
% check if derivative has been included => if so each regressor will be
% modelled as baseline + temporal and/or spatial derivative
add_drv = generalPrm.add_drv;
switch add_drv
    case 0
        disp('no derivative included');
    case 1
        disp('Use of temporal derivative => be aware that each regressor will be doubled in the matrix list.');
    case 2
        disp('Use of temporal and spatial derivative => be aware that each regressor will be tripled in the matrix list');
end

% zscore variables per run?
z_perRun = generalPrm.zPerRun;
switch z_perRun
    case 1
        disp('all variables are zscored per run');
end

% orthogonalize variables or not?
orth_vars = generalPrm.orth_vars;
switch orth_vars
    case 0
        disp('Variables are not orthogonalized');
    case 1
        disp('Variables are orthogonalized.');
end

%% loop through both tasks
Epm = {'Ep','Em'};
for iEpm = 1:2
    task_id_nm = Epm{iEpm};
    switch task_id_nm
        case 'Ep'
            disp('** Physical effort design **');
        case 'Em'
            disp('** Mental effort design **');
    end
    
    %% Ep: fixation cross (across the whole task)
    if ~strcmp(GLMprm.model_onset.(task_id_nm).cross,'none')
        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
        reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'ONSET fixation cross';
        disp([num2str(n_regs.(task_id_nm)),') ONSET cross - ',GLMprm.model_onset.(task_id_nm).cross,' ']);
        % if derivative added => add derivatives
        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
    end % cross
    
    %% Ep: choice period
    if ~strcmp(GLMprm.model_onset.(task_id_nm).choice,'none')
        % check if trials are split or not
        if GLMprm.choice.(task_id_nm).RPpool == 1 % pool reward and punishment trials
            n_RP = 1;
            RP = {'RP'};
        elseif GLMprm.choice.(task_id_nm).RPpool == 0 % split reward and punishment
            n_RP = 2;
            RP = {'R','P'};
        end % RP pool
        
        % loop through conditions for choice period
        for iRP = 1:n_RP
            RP_nm = RP{iRP};
            
            % choice onset
            n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
            reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['ONSET choice display options - ',RP_nm];
            disp([num2str(n_regs.(task_id_nm)),') ONSET choice display options ',RP_nm,' - ',GLMprm.model_onset.(task_id_nm).choice,' ']);
            % if derivative added => add derivatives
            n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
            
            %choice regressors
            
            % Reward > Punishment
            if GLMprm.choice.(task_id_nm).(RP_nm).R_vs_P == 1
                n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG choice - R min P';
                disp([num2str(n_regs.(task_id_nm)),') choice - Reward>Punishment ']);
                % if derivative added => add derivatives
                n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
            end
            
            % money left
            switch GLMprm.choice.(task_id_nm).(RP_nm).money_left
                case 1
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                    reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG choice - money left';
                    disp([num2str(n_regs.(task_id_nm)),') choice - money left ']);
                    % if derivative added => add derivatives
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
            end
            
            % money right
            switch GLMprm.choice.(task_id_nm).(RP_nm).money_right
                case 1
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                    reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG choice - money right';
                    disp([num2str(n_regs.(task_id_nm)),') choice - money right ']);
                    % if derivative added => add derivatives
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
            end
            
            % money chosen option
            switch GLMprm.choice.(task_id_nm).(RP_nm).money_chosen
                case 1
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                    reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG choice - money chosen';
                    disp([num2str(n_regs.(task_id_nm)),') choice - money chosen ']);
                    % if derivative added => add derivatives
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                case 2
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                    reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG choice - money chosen';
                    disp([num2str(n_regs.(task_id_nm)),') choice - |money chosen| ']);
                    % if derivative added => add derivatives
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
            end
            
            % money unchosen option
            switch GLMprm.choice.(task_id_nm).(RP_nm).money_unchosen
                case 1
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                    reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG choice - money unchosen';
                    disp([num2str(n_regs.(task_id_nm)),') choice - money unchosen ']);
                    % if derivative added => add derivatives
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
            end
            
            % money chosen - money unchosen
            switch GLMprm.choice.(task_id_nm).(RP_nm).money_ch_min_unch
                case 1
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                    reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG choice - money ch-unch';
                    disp([num2str(n_regs.(task_id_nm)),') choice - money ch-unch ']);
                    % if derivative added => add derivatives
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
            end
            
            % money sum of both options
            switch GLMprm.choice.(task_id_nm).(RP_nm).money_sum
                case 1
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                    reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG choice - money sum';
                    disp([num2str(n_regs.(task_id_nm)),') choice - money sum ']);
                    % if derivative added => add derivatives
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
            end
            
            % effort left
            switch GLMprm.choice.(task_id_nm).(RP_nm).E_left
                case 1
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                    reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG choice - effort left';
                    disp([num2str(n_regs.(task_id_nm)),') choice - effort left (level) ']);
                    % if derivative added => add derivatives
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                case 2
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                    reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG choice - effort left';
                    disp([num2str(n_regs.(task_id_nm)),') choice - effort left (durations) ']);
                    % if derivative added => add derivatives
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
            end
            
            % effort right
            switch GLMprm.choice.(task_id_nm).(RP_nm).E_right
                case 1
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                    reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG choice - effort right';
                    disp([num2str(n_regs.(task_id_nm)),') choice - effort right (level) ']);
                    % if derivative added => add derivatives
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                case 2
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                    reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG choice - effort right';
                    disp([num2str(n_regs.(task_id_nm)),') choice - effort right (durations) ']);
                    % if derivative added => add derivatives
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
            end
            
            % effort chosen option
            switch GLMprm.choice.(task_id_nm).(RP_nm).E_chosen
                case 1
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                    reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG choice - effort chosen';
                    disp([num2str(n_regs.(task_id_nm)),') choice - effort chosen (levels) ']);
                    % if derivative added => add derivatives
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                case 2
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                    reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG choice - effort chosen';
                    disp([num2str(n_regs.(task_id_nm)),') choice - effort chosen (durations) ']);
                    % if derivative added => add derivatives
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
            end
            
            % effort unchosen option
            switch GLMprm.choice.(task_id_nm).(RP_nm).money_unchosen
                case 1
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                    reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG choice - effort unchosen';
                    disp([num2str(n_regs.(task_id_nm)),') choice - effort unchosen (levels) ']);
                    % if derivative added => add derivatives
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                case 2
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                    reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG choice - effort unchosen';
                    disp([num2str(n_regs.(task_id_nm)),') choice - effort unchosen (durations) ']);
                    % if derivative added => add derivatives
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
            end
            
            % effort chosen - unchosen option
            switch GLMprm.choice.(task_id_nm).(RP_nm).E_ch_min_unch
                case 1
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                    reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG choice - effort ch-unch';
                    disp([num2str(n_regs.(task_id_nm)),') choice - effort ch-unch (levels) ']);
                    % if derivative added => add derivatives
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                case 2
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                    reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG choice - effort ch-unch';
                    disp([num2str(n_regs.(task_id_nm)),') choice - effort ch-unch (durations) ']);
                    % if derivative added => add derivatives
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
            end
            
            % effort sum of both options
            switch GLMprm.choice.(task_id_nm).(RP_nm).E_sum
                case 1
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                    reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG choice - effort sum';
                    disp([num2str(n_regs.(task_id_nm)),') choice - effort sum (levels) ']);
                    % if derivative added => add derivatives
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                case 2
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                    reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG choice - effort sum';
                    disp([num2str(n_regs.(task_id_nm)),') choice - effort sum (durations) ']);
                    % if derivative added => add derivatives
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
            end
            
            % confidence
            switch GLMprm.choice.(task_id_nm).(RP_nm).confidence
                case 1
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                    reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG choice - confidence';
                    disp([num2str(n_regs.(task_id_nm)),') choice - confidence (levels) ']);
                    % if derivative added => add derivatives
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
            end
            
            % RT
            switch GLMprm.choice.(task_id_nm).(RP_nm).RT
                case 1
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                    reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG choice - RT';
                    disp([num2str(n_regs.(task_id_nm)),') choice - RT (raw) ']);
                    % if derivative added => add derivatives
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                case 2
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                    reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG choice - RT';
                    disp([num2str(n_regs.(task_id_nm)),') choice - RT (zscored per run) ']);
                    % if derivative added => add derivatives
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                case 3
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                    reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG choice - RT';
                    disp([num2str(n_regs.(task_id_nm)),') choice - RT (zscored per subject ie across all runs) ']);
                    % if derivative added => add derivatives
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
            end
            
        end % loop reward/punishment
    end % choice onset
    
    %% Ep: chosen period
    if ~strcmp(GLMprm.model_onset.(task_id_nm).chosen,'none')
        % check if trials are split or not
        if GLMprm.chosen.(task_id_nm).RPpool == 1 % pool reward and punishment trials
            n_RP = 1;
            RP = {'RP'};
        elseif GLMprm.chosen.(task_id_nm).RPpool == 0 % split reward and punishment
            n_RP = 2;
            RP = {'R','P'};
        end % RP pool
        
        % loop through conditions for chosen period
        for iRP = 1:n_RP
            RP_nm = RP{iRP};
            
            % chosen onset
            n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
            reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['ONSET chosen option display - ',RP_nm];
            disp([num2str(n_regs.(task_id_nm)),') ONSET chosen option display ',RP_nm,' - ',GLMprm.model_onset.(task_id_nm).chosen,' ']);
            % if derivative added => add derivatives
            n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
            
            % chosen regressors
            % money chosen
            switch GLMprm.chosen.(task_id_nm).(RP_nm).money_chosen
                case 1
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                    reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG chosen - money chosen';
                    disp([num2str(n_regs.(task_id_nm)),') chosen - money chosen (levels) ']);
                    % if derivative added => add derivatives
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                case 2
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                    reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG chosen - money chosen';
                    disp([num2str(n_regs.(task_id_nm)),') chosen - |money chosen| (levels) ']);
                    % if derivative added => add derivatives
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
            end
            
            % effort chosen
            switch GLMprm.chosen.(task_id_nm).(RP_nm).E_chosen
                case 1
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                    reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG chosen - effort chosen';
                    disp([num2str(n_regs.(task_id_nm)),') chosen - effort chosen (levels) ']);
                    % if derivative added => add derivatives
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
            end
            
            % money unchosen
            switch GLMprm.chosen.(task_id_nm).(RP_nm).money_unchosen
                case 1
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                    reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG chosen - money unchosen';
                    disp([num2str(n_regs.(task_id_nm)),') chosen - money unchosen (levels) ']);
                    % if derivative added => add derivatives
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
            end
            
            % effort unchosen
            switch GLMprm.chosen.(task_id_nm).(RP_nm).E_unchosen
                case 1
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                    reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG chosen - effort unchosen';
                    disp([num2str(n_regs.(task_id_nm)),') chosen - effort unchosen (levels) ']);
                    % if derivative added => add derivatives
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
            end
            
            % confidence
            switch GLMprm.chosen.(task_id_nm).(RP_nm).confidence
                case 1
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                    reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG chosen - confidence';
                    disp([num2str(n_regs.(task_id_nm)),') chosen - confidence (levels) ']);
                    % if derivative added => add derivatives
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
            end
            
        end % loop reward/punishment
    end % chosen onset
    
    %% Ep: effort performance period
    if ~strcmp(GLMprm.model_onset.(task_id_nm).Eperf,'none')
        % check if trials are split or not
        if GLMprm.Eperf.(task_id_nm).RPpool == 1 % pool reward and punishment trials
            n_RP = 1;
            RP = {'RP'};
        elseif GLMprm.Eperf.(task_id_nm).RPpool == 0 % split reward and punishment
            n_RP = 2;
            RP = {'R','P'};
        end % RP pool
        
        % loop through conditions for choice period
        for iRP = 1:n_RP
            RP_nm = RP{iRP};
            
            % effort period onset
            n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
            reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['ONSET effort period - ',RP_nm];
            disp([num2str(n_regs.(task_id_nm)),') ONSET effort period ',RP_nm,' - ',GLMprm.model_onset.(task_id_nm).Eperf,' ']);
            % if derivative added => add derivatives
            n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
            
            % effort period regressors
            % money chosen
            switch GLMprm.Eperf.(task_id_nm).(RP_nm).money_chosen
                case 1
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                    reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG effort - money chosen';
                    disp([num2str(n_regs.(task_id_nm)),') effort period - money chosen (levels) ']);
                    % if derivative added => add derivatives
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                case 2
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                    reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG effort - money chosen';
                    disp([num2str(n_regs.(task_id_nm)),') effort period - |money chosen| (levels) ']);
                    % if derivative added => add derivatives
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
            end
            
            % effort chosen
            switch GLMprm.Eperf.(task_id_nm).(RP_nm).E_chosen
                case 1
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                    reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG effort - effort chosen';
                    disp([num2str(n_regs.(task_id_nm)),') effort period - effort chosen (levels) ']);
                    % if derivative added => add derivatives
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
            end
            
            % force peak
            switch task_id_nm
                case 'Ep'
                    switch GLMprm.Eperf.(task_id_nm).(RP_nm).F_peak
                        case 1
                            n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                            reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG effort - force peak';
                            disp([num2str(n_regs.(task_id_nm)),') effort period - force peak ']);
                            % if derivative added => add derivatives
                            n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    end
                    
                    % force integral
                    switch GLMprm.Eperf.(task_id_nm).(RP_nm).F_integral
                        case 1
                            n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                            reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG effort - force integral';
                            disp([num2str(n_regs.(task_id_nm)),') effort period - force integral ']);
                            % if derivative added => add derivatives
                            n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    end
                    
                case 'Em'
                    % average RT for all numbers of each trial
                    switch GLMprm.Eperf.(task_id_nm).(RP_nm).RT_avg
                        case 1
                            n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                            reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG effort - average RT effort';
                            disp([num2str(n_regs.(task_id_nm)),') effort period - average RT effort ']);
                            % if derivative added => add derivatives
                            n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    end
                    
                    % number of errors
                    switch GLMprm.Eperf.(task_id_nm).(RP_nm).n_errors
                        case 1
                            n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                            reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG effort - number of errors';
                            disp([num2str(n_regs.(task_id_nm)),') effort period - number of errors ']);
                            % if derivative added => add derivatives
                            n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    end
            end
            
            % RT first answer
            switch GLMprm.Eperf.(task_id_nm).(RP_nm).RT_1stAnswer
                case 1
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                    reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG effort - RT 1st answer';
                    disp([num2str(n_regs.(task_id_nm)),') effort period - RT 1st answer (raw) ']);
                    % if derivative added => add derivatives
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
            end
            
        end % loop reward/punishment
    end % choice onset
    
    %% Ep: feedback period
    if ~strcmp(GLMprm.model_onset.(task_id_nm).fbk,'none')
        % check if trials are split or not
        if GLMprm.fbk.(task_id_nm).RPpool == 1 % pool reward and punishment trials
            n_RP = 1;
            RP = {'RP'};
        elseif GLMprm.fbk.(task_id_nm).RPpool == 0 % split reward and punishment
            n_RP = 2;
            RP = {'R','P'};
        end % RP pool
        
        % loop through conditions for choice period
        for iRP = 1:n_RP
            RP_nm = RP{iRP};
            
            % feedback onset
            n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
            reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['ONSET feedback - ',RP_nm];
            disp([num2str(n_regs.(task_id_nm)),') ONSET feedback ',RP_nm,' - ',GLMprm.model_onset.(task_id_nm).fbk,' ']);
            % if derivative added => add derivatives
            n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
            
            % feedback regressors
            % win vs loss
            switch GLMprm.fbk.(task_id_nm).(RP_nm).win_vs_loss
                case 1
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                    reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG feedback - win min loss';
                    disp([num2str(n_regs.(task_id_nm)),') feedback - win min loss ']);
                    % if derivative added => add derivatives
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
            end
            % money obtained
            switch GLMprm.fbk.(task_id_nm).(RP_nm).money_obtained
                case 1
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                    reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG feedback - money obtained';
                    disp([num2str(n_regs.(task_id_nm)),') feedback - money obtained (levels) ']);
                    % if derivative added => add derivatives
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
            end
            
            % effort performed
            switch GLMprm.fbk.(task_id_nm).(RP_nm).E_made
                case 1
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                    reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG feedback - effort performed';
                    disp([num2str(n_regs.(task_id_nm)),') feedback - effort performed (levels) ']);
                    % if derivative added => add derivatives
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
            end
            
            % confidence
            switch GLMprm.fbk.(task_id_nm).(RP_nm).confidence
                case 1
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                    reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG feedback - confidence';
                    disp([num2str(n_regs.(task_id_nm)),') feedback - confidence (levels) ']);
                    % if derivative added => add derivatives
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
            end
            
        end % loop reward/punishment
    end % choice onset
    
    %% Ep: movement regressors
    n_mvmt = 6;
    for iMvmt = 1:n_mvmt
        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
        reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'movement';
    end
    disp([num2str(n_regs.(task_id_nm)-n_mvmt+1),'-',num2str(n_regs.(task_id_nm)),') movement regressors ']);
    warning('need to check whether including the derivative also impacts the movement regressors or not. My guess is no but recheck please');
    
end % physical/mental loop
end % function