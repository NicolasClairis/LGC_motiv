function [reg_names, n_regs] = GLM_details(GLM, dispRegs)
% [Ep_regs, Em_regs] = GLM_details(GLM, dispRegs)
%GLM_details will provide the number of regressors, the text with the
%details of the GLM, the order of the regressors in the current GLM for
%each task.
%
% INPUTS
% GLM: number of the GLM
%
% dispRegs: display information on regressors (1) or not (0)
%
% OUTPUTS
% reg_names: structure with physical (Ep) and mental (Em) regressors for
% each run
%
% n_regs: structure with number of regressors for each physical (Ep) and
% each mental (Em) run

%% display parameters by default if dispRegs left empty
if ~exist('dispRegs','var') ||...
        isempty(dispRegs) ||...
        ~ismember(dispRegs,[0,1])
    dispRegs = 1;
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
[reg_names.Ep, reg_names.Em] = deal( {} );
n_regs.Ep = 0;
n_regs.Em = 0;

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

%% loop through both tasks
Epm = {'Ep','Em'};
for iEpm = 1:length(Epm)
    task_id_nm = Epm{iEpm};
    switch task_id_nm
        case 'Ep'
            dispRegFn('** Physical effort design **',dispRegs);
        case 'Em'
            dispRegFn('** Mental effort design **',dispRegs);
    end
    
    %% all fixation crosses (pool of pre-choice and pre-effort cross)
    if ~strcmp(GLMprm.model_onset.(task_id_nm).allCrosses,'none')
        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
        reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'ONSET fixation cross';
        dispRegFn([num2str(n_regs.(task_id_nm)),') ONSET fixation cross: ',GLMprm.model_onset.(task_id_nm).allCrosses,' '],dispRegs);
        % if derivative added => add derivatives
        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
    end
    
    %% pre-choice fixation cross
    if ~strcmp(GLMprm.model_onset.(task_id_nm).preChoiceCross,'none')
        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
        reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'ONSET preChoice white fixation cross';
        dispRegFn([num2str(n_regs.(task_id_nm)),') ONSET preChoice white cross: ',GLMprm.model_onset.(task_id_nm).preChoiceCross,' '],dispRegs);
        % if derivative added => add derivatives
        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
        
        % RT (first regressor)
        switch GLMprm.preChoiceCross.(task_id_nm).RT
            case 4
                n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG preChoice Cross RP E: RT';
                dispRegFn([num2str(n_regs.(task_id_nm)),') preChoice Cross: RT (raw) '],dispRegs);
                % if derivative added => add derivatives
                n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
            case 5
                n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG preChoice Cross RP E: RT';
                dispRegFn([num2str(n_regs.(task_id_nm)),') preChoice Cross: RT (zscored per run) '],dispRegs);
                % if derivative added => add derivatives
                n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
            case 6
                n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG preChoice Cross RP E: RT';
                dispRegFn([num2str(n_regs.(task_id_nm)),') preChoice Cross: RT (zscored per subject ie across all runs) '],dispRegs);
                % if derivative added => add derivatives
                n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
        end
        
        % binary variable indicating when choice = high effort option
        if GLMprm.preChoiceCross.(task_id_nm).choiceHighE == 1
            n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
            reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG preChoice Cross RP E: choice = highE';
            dispRegFn([num2str(n_regs.(task_id_nm)),') preChoice Cross: choice hE '],dispRegs);
            % if derivative added => add derivatives
            n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
        end
        
        % effort chosen
        if GLMprm.preChoiceCross.(task_id_nm).E_chosen == 1
            n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
            reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG preChoice Cross RP E: effort chosen';
            dispRegFn([num2str(n_regs.(task_id_nm)),') preChoice Cross: effort chosen (levels) '],dispRegs);
            % if derivative added => add derivatives
            n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
        end
        
        % RT (last regressor)
        switch GLMprm.preChoiceCross.(task_id_nm).RT
            case 1
                n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG preChoice Cross RP E: RT';
                dispRegFn([num2str(n_regs.(task_id_nm)),') preChoice Cross: RT (raw) '],dispRegs);
                % if derivative added => add derivatives
                n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
            case 2
                n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG preChoice Cross RP E: RT';
                dispRegFn([num2str(n_regs.(task_id_nm)),') preChoice Cross: RT (zscored per run) '],dispRegs);
                % if derivative added => add derivatives
                n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
            case 3
                n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'REG preChoice Cross RP E: RT';
                dispRegFn([num2str(n_regs.(task_id_nm)),') preChoice Cross: RT (zscored per subject ie across all runs) '],dispRegs);
                % if derivative added => add derivatives
                n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
        end
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
                
                %% choice onset
                n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['ONSET choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm];
                dispRegFn([num2str(n_regs.(task_id_nm)),') ONSET choice display options ',...
                    RP_dispChoice_nm,' ',splitE_dispChoice_nm,': ',GLMprm.model_onset.(task_id_nm).choice,' '],dispRegs);
                % if derivative added => add derivatives
                n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                
                %% choice regressors
                
                % RT (first regressor)
                switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).RT
                    case 4
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',...
                            RP_dispChoice_nm,' ',splitE_dispChoice_nm,': RT'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: RT (raw) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 5
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',...
                            RP_dispChoice_nm,' ',splitE_dispChoice_nm,': RT'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: RT (zscored per run) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 6
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',...
                            RP_dispChoice_nm,' ',splitE_dispChoice_nm,': RT'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: RT (zscored per subject ie across all runs) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
                
                % Reward > Punishment
                if GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).R_vs_P == 1
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                    reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': R-P'];
                    dispRegFn([num2str(n_regs.(task_id_nm)),') choice: Reward>Punishment '],dispRegs);
                    % if derivative added => add derivatives
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
                
                % binary variable indicating when choice = high effort option
                if GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).choiceHighE == 1
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                    reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': choice = highE'];
                    dispRegFn([num2str(n_regs.(task_id_nm)),') choice: choice hE '],dispRegs);
                    % if derivative added => add derivatives
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
                
                % reward variable option
                switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).R_varOption
                    case 0
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': R amount highE'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: R amount hE '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': R level highE'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: R level hE '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    otherwise
                        error('not ready yet');
                end
                
                % reward chosen option
                switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).R_chosen
                    case 0
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': R amount chosen'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: R amount chosen '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': R level chosen'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: R level chosen '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    otherwise
                        error('not ready yet');
                end
                
                % punishment variable option
                switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).P_varOption
                    case 0
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': P amount highE'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: P amount hE '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': P level highE'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: P level hE '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    otherwise
                        error('not ready yet');
                end
                
                % punishment chosen option
                switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).P_chosen
                    case 0
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': P amount chosen'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: P amount chosen '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': P level chosen'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: P level chosen '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    otherwise
                        error('not ready yet');
                end
                
                % money left
                switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_left
                    case 0
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money left'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: money left (amounts)'],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money left'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: |money left| (amounts)'],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 3
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money left'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: money left (levels)'],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 4
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money left'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: |money left| (levels)'],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    otherwise
                        error('not ready yet');
                end
                
                % money right
                switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_right
                    case 0
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money right'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: money right (amounts)'],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money right'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: |money right| (amounts)'],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 3
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money right'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: money right (levels)'],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 4
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money right'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: |money right| (levels)'],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    otherwise
                        error('not ready yet');
                end
                
                % money chosen option
                switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_chosen
                    case 0
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money chosen'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: money chosen (amounts)'],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money chosen'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: |money chosen| (amounts)'],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 3
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money chosen'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: money chosen (levels)'],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 4
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money chosen'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: |money chosen| (levels)'],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    otherwise
                        error('not ready yet');
                end
                
                % money unchosen option
                switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_unchosen
                    case 0
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money unchosen'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: money unchosen (amounts)'],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money unchosen'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: |money unchosen| (amounts)'],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 3
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money unchosen'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: money unchosen (levels)'],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 4
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money unchosen'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: |money unchosen| (levels)'],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    otherwise
                        error('not ready yet');
                end
                
                % money associated to the option which varies (the non-default
                % option)
                switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_varOption
                    case 0
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money non-default option'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: money non-default option (amount) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money non-default option'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: |money| non-default option (amount) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 3
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money non-default option'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: money non-default option (level) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 4
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money non-default option'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: |money| non-default option (level) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    otherwise
                        error('not ready yet');
                end
                
                % money chosen - money unchosen
                switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_ch_min_unch
                    case 0
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money ch-unch'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: money ch-unch '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    otherwise
                        error('not ready yet');
                end
                
                % money chosen - money default
                switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_ch_min_fixOption
                    case 0
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money ch-def'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: money chosen-default (amount) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money ch-def'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: money chosen-default (level) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    otherwise
                        error('not ready yet');
                end
                
                % money sum of both options
                switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_sum
                    case 0
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money sum'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: money sum '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    otherwise
                        error('not ready yet');
                end
                
                % effort left
                switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).E_left
                    case 0
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': effort left'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: effort left (level) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': effort left'];
                        switch task_id_nm
                            case 'Ep'
                                dispRegFn([num2str(n_regs.(task_id_nm)),') choice: effort left (durations) '],dispRegs);
                            case 'Em'
                                dispRegFn([num2str(n_regs.(task_id_nm)),') choice: effort left (nb answers to give) '],dispRegs);
                        end
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    otherwise
                        error('not ready yet');
                end
                
                % effort right
                switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).E_right
                    case 0
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': effort right'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: effort right (level) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': effort right'];
                        switch task_id_nm
                            case 'Ep'
                                dispRegFn([num2str(n_regs.(task_id_nm)),') choice: effort right (durations) '],dispRegs);
                            case 'Em'
                                dispRegFn([num2str(n_regs.(task_id_nm)),') choice: effort right (nb answers to give) '],dispRegs);
                        end
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    otherwise
                        error('not ready yet');
                end
                
                % effort chosen option
                switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).E_chosen
                    case 0
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': effort chosen'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: effort chosen (levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': effort chosen'];
                        switch task_id_nm
                            case 'Ep'
                                dispRegFn([num2str(n_regs.(task_id_nm)),') choice: effort chosen (durations) '],dispRegs);
                            case 'Em'
                                dispRegFn([num2str(n_regs.(task_id_nm)),') choice: effort chosen (nb answers to give) '],dispRegs);
                        end
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    otherwise
                        error('not ready yet');
                end
                
                % effort unchosen option
                switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).E_unchosen
                    case 0
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': effort unchosen'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: effort unchosen (levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': effort unchosen'];
                        switch task_id_nm
                            case 'Ep'
                                dispRegFn([num2str(n_regs.(task_id_nm)),') choice: effort unchosen (durations) '],dispRegs);
                            case 'Em'
                                dispRegFn([num2str(n_regs.(task_id_nm)),') choice: effort unchosen (nb answers to give) '],dispRegs);
                        end
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    otherwise
                        error('not ready yet');
                end
                
                % effort non-default option
                switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).E_varOption
                    case 0
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': effort non-default option'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: effort non-default option (levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': effort non-default option'];
                        switch task_id_nm
                            case 'Ep'
                                dispRegFn([num2str(n_regs.(task_id_nm)),') choice: effort non-default option (durations) '],dispRegs);
                            case 'Em'
                                dispRegFn([num2str(n_regs.(task_id_nm)),') choice: effort non-default option (nb answers to give) '],dispRegs);
                        end
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    otherwise
                        error('not ready yet');
                end
                
                % effort chosen - unchosen option
                switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).E_ch_min_unch
                    case 0
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': effort ch-unch'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: effort ch-unch (levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': effort ch-unch'];
                        switch task_id_nm
                            case 'Ep'
                                dispRegFn([num2str(n_regs.(task_id_nm)),') choice: effort ch-unch (durations) '],dispRegs);
                            case 'Em'
                                dispRegFn([num2str(n_regs.(task_id_nm)),') choice: effort ch-unch (nb answers to give) '],dispRegs);
                        end
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    otherwise
                        error('not ready yet');
                end
                
                % effort chosen - fixed low effort option
                switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).E_ch_min_fixOption
                    case 0
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': effort ch-def'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: effort ch-def (levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    otherwise
                        error('not ready yet');
                end
                
                % effort sum of both options
                switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).E_sum
                    case 0
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': effort sum'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: effort sum (levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': effort sum'];
                        switch task_id_nm
                            case 'Ep'
                                dispRegFn([num2str(n_regs.(task_id_nm)),') choice: effort sum (durations) '],dispRegs);
                            case 'Em'
                                dispRegFn([num2str(n_regs.(task_id_nm)),') choice: effort sum (nb answers to give) '],dispRegs);
                        end
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    otherwise
                        error('not ready yet');
                end
                
                % (money)*(effort) for high effort option
                switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_level_x_E_varOption
                    case 0
                    case 1 % (money level)*(effort level) high E option
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money x effort non-default'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: money x effort non-default option '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    otherwise
                        error('not ready yet');
                end
                
                % (money)*(effort) for chosen option
                switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).money_level_x_E_chosen
                    case 0
                    case 1 % (money level)*(effort level) chosen option
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': money x effort chosen'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: money x effort chosen '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    otherwise
                        error('not ready yet');
                end
                
                % (reward)*(effort) for high effort option
                switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).R_level_x_E_varOption
                    case 0
                    case 1 % (R level)*(effort level) high E option
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': R x E non-default'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: R x E non-default option '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    otherwise
                        error('not ready yet');
                end
                
                % (reward)*(effort) for chosen option
                switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).R_level_x_E_chosen
                    case 0
                    case 1 % (R level)*(effort level) chosen option
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': R x E chosen'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: R x E chosen '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    otherwise
                        error('not ready yet');
                end
                
                % (punishment)*(effort) for high effort option
                switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).P_level_x_E_varOption
                    case 0
                    case 1 % (P level)*(effort level) high E option
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': P x E non-default'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: P x E non-default option '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    otherwise
                        error('not ready yet');
                end
                
                % (punishment)*(effort) for chosen option
                switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).P_level_x_E_chosen
                    case 0
                    case 1 % (P level)*(effort level) chosen option
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': P x E chosen'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: P x E chosen '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    otherwise
                        error('not ready yet');
                end
                
                % net value chosen option
                switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).NV_chosen
                    case 0
                    case {1,3} % NV(chosen)-NV(unchosen)
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': NVch-NVunch'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: NVch-NVunch '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2 % p(chosen)
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': p(chosen)'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: p(chosen) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    otherwise
                        error('not ready yet');
                end
                
                % net value variable option
                switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).NV_varOption
                    case 0
                    case 1 % NV(high E) - NV(low E)
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': net value high E - low E'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: net value high E - low E '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2 % |NV(high E) - NV(low E)|
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': |net value high E - low E|'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: |net value high E - low E| '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 3 % p(choice = high E)
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': p(choice = high E)'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: p(choice = high E) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 4 % NV(high E) - NV(low E) + bias
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': net value high E - low E + bias'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: net value high E - low E + bias '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 5 % |NV(high E) - NV(low E) + bias|
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': |net value high E - low E + bias|'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: |net value high E - low E + bias| '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    otherwise
                        error('not ready yet');
                end
                
                % net value variable option bis
                switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).NV_varOption_bis
                    case 0
                    case 1 % NV(high E) - NV(low E)
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': net value high E - low E'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: net value high E - low E '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2 % |NV(high E) - NV(low E)|
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': |net value high E - low E|'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: |net value high E - low E| '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 3 % p(choice = high E)
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': p(choice = high E)'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: p(choice = high E) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 4 % NV(high E) - NV(low E) + bias
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': net value high E - low E + bias'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: net value high E - low E + bias '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 5 % |NV(high E) - NV(low E) + bias|
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': |net value high E - low E + bias|'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: |net value high E - low E + bias| '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    otherwise
                        error('not ready yet');
                end
                
                if strcmp(task_id_nm, 'Ep') % physical effort only
                    % area under the curve of the force that is gonna be produced
                    switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).F_integral
                        case 0
                        case {1,3}
                            n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                            reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': effort integral'];
                            dispRegFn([num2str(n_regs.(task_id_nm)),') choice: effort integral '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                        case {2,4}
                            n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                            reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': effort integral overshoot'];
                            dispRegFn([num2str(n_regs.(task_id_nm)),') choice: effort integral overshoot '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    otherwise
                        error('not ready yet');
                    end
                    
                    % fatigue
                    switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).fatigue
                        case 0
                        case 1
                            n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                            reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': fatigue'];
                            dispRegFn([num2str(n_regs.(task_id_nm)),') choice: fatigue '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % Ech*fatigue
                    switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).Ech_x_fatigue
                        case 0
                        case 1
                            n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                            reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': Ech_x_fatigue'];
                            dispRegFn([num2str(n_regs.(task_id_nm)),') choice: Ech_x_fatigue '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                end % physical effort filter
                
                if strcmp(task_id_nm,'Em') % mental effort only
                    % efficacy of the next trial
                    switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).efficacy
                        case 0
                        case {1,2,3,4}
                            n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                            reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': efficacy'];
                            dispRegFn([num2str(n_regs.(task_id_nm)),') choice: efficacy '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % efficacy during the previous trial
                    switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).prevEfficacy
                        case 0
                        case {1,2,3,4}
                            n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                            reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': previous efficacy'];
                            dispRegFn([num2str(n_regs.(task_id_nm)),') choice: previous efficacy '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                    
                    % (effort chosen)*(efficacy during the previous trial)
                    switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).Ech_x_prevEfficacy
                        case 0
                        case {1,2,3,4}
                            n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                            reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': (Ech)x(previous efficacy)'];
                            dispRegFn([num2str(n_regs.(task_id_nm)),') choice: (Ech)x(previous efficacy) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                        otherwise
                            error('not ready yet');
                    end
                end % mental effort filter
                
                % trial number
                switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).trialN
                        case 0
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': trial number'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: trial number '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': (trial number)x(Echosen-Eunchosen) (E levels)'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: (trial number)x(effort chosen-effort unchosen) (effort levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 3
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': trial number x (EnonDef-Edef) (E levels)'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: (trial number)x(effort non-default - effort default) (effort levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 4
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': (trial number)x(EnonDef) (E levels)'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: (trial number)x(effort non-default) (effort levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    otherwise
                        error('not ready yet');
                end
                
                % confidence
                switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).confidence
                        case 0
                    case 1 % confidence ratings 0/1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': confidence'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: confidence (levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case {2,3,4} % confidence inferred by the model
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': confidence'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: confidence (inferred by the model) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    otherwise
                        error('not ready yet');
                end
                
                % RT (last regressor)
                switch GLMprm.choice.(task_id_nm).(RP_dispChoice_nm).(splitE_dispChoice_nm).RT
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': RT'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: RT (raw) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': RT'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: RT (zscored per run) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 3
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG choice ',RP_dispChoice_nm,' ',splitE_dispChoice_nm,': RT'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') choice: RT (zscored per subject ie across all runs) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
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
                
                %% chosen onset
                n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['ONSET chosen ',...
                    RP_dispChosen_nm,' ',splitE_dispChosen_nm];
                dispRegFn([num2str(n_regs.(task_id_nm)),') ONSET chosen option display ',...
                    RP_dispChosen_nm,' ',splitE_dispChosen_nm,': ',GLMprm.model_onset.(task_id_nm).chosen,' '],dispRegs);
                % if derivative added => add derivatives
                n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                
                %% chosen regressors
                
                % RT (first regressor)
                switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).RT
                    case 4
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': RT'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: RT (raw) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 5
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': RT'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: RT (zscored per run) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 6
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': RT'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: RT (zscored per subject ie across all runs) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
                
                % reward/punishment trial
                switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).R_vs_P
                        case 0
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': R-P'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: R vs P '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
                
                % binary variable indicating when choice = high effort option
                if GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).choiceHighE == 1
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                    reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': choice = highE'];
                    dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: choice hE '],dispRegs);
                    % if derivative added => add derivatives
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
                
                % reward variable option
                switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).R_varOption
                        case 0
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': R amount highE'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: R amount hE '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': R level highE'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: R level hE '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    otherwise
                        error('not ready yet');
                end
                
                % reward chosen option
                switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).R_chosen
                    case 0
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': R amount chosen'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: R amount chosen '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': R level chosen'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: R level chosen '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    otherwise
                        error('not ready yet');
                end
                
                % punishment variable option
                switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).P_varOption
                        case 0
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': P amount highE'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: P amount hE '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': P level highE'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: P level hE '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    otherwise
                        error('not ready yet');
                end
                
                % punishment chosen option
                switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).P_chosen
                        case 0
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': P amount chosen'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: P amount chosen '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': P level chosen'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: P level chosen '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    otherwise
                        error('not ready yet');
                end
                
                % money left
                switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_left
                        case 0
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money left'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: money left (amounts) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money left'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: |money left| (amounts) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 3
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money left'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: money left (levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 4
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money left'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: |money left| (levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
                
                % money right
                switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_right
                        case 0
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money right'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: money right (amounts) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money right'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: |money right| (amounts) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 3
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money right'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: money right (levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 4
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money right'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: |money right| (levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
                
                % money chosen
                switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_chosen
                        case 0
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money chosen'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: money chosen (amounts) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money chosen'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: |money chosen| (amounts) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 3
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money chosen'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: money chosen (levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 4
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money chosen'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: |money chosen| (levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
                
                % money unchosen
                switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_unchosen
                        case 0
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money unchosen'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: money unchosen (amounts) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money unchosen'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: |money unchosen| (amounts) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 3
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money unchosen'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: money unchosen (levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 4
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money unchosen'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: |money unchosen| (amounts) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
                
                % money associated to the option which varies (the non-default
                % option)
                switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_varOption
                        case 0
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money non-default option'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: money non-default option (amount) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money non-default option'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: |money| non-default option (amount) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 3
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money non-default option'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: money non-default option (level) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 4
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money non-default option'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: |money| non-default option (level) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
                
                % money chosen - money unchosen
                switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_ch_min_unch
                        case 0
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money ch-unch'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: money ch-unch '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
                
                % money chosen - money default
                switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_ch_min_fixOption
                        case 0
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money ch-def'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: money chosen-default (amount) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money ch-def'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: money chosen-default (level) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
                
                % money sum of both options
                switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_sum
                        case 0
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money sum'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: money sum '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
                
                % effort chosen
                switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).E_chosen
                        case 0
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': effort chosen'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: effort chosen (levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
                
                % effort unchosen
                switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).E_unchosen
                        case 0
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': effort unchosen'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: effort unchosen (levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
                
                % effort non-default option
                switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).E_varOption
                        case 0
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': effort non-default option'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: effort non-default option (levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
                
                % effort chosen - effort unchosen
                switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).E_ch_min_unch
                        case 0
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': effort ch-unch'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: effort ch-unch (levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': effort ch-unch'];
                        switch task_id_nm
                            case 'Ep'
                                dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: effort ch-unch (durations) '],dispRegs);
                            case 'Em'
                                dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: effort ch-unch (nb answers to give) '],dispRegs);
                        end
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
                
                % effort chosen - fixed low effort option
                switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).E_ch_min_fixOption
                        case 0
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': effort ch-def'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: effort ch-def (levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
                
                % effort sum of both options
                switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).E_sum
                        case 0
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': effort sum'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: effort sum (levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': effort sum'];
                        switch task_id_nm
                            case 'Ep'
                                dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: effort sum (durations) '],dispRegs);
                            case 'Em'
                                dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: effort sum (nb answers to give) '],dispRegs);
                        end
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
                
                % (money)*(effort) for high effort option
                switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_level_x_E_varOption
                    case 0
                    case 1 % (money level)*(effort level) high E option
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money x effort non-default'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: money x effort non-default option '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    otherwise
                        error('not ready yet');
                end
                
                % (money)*(effort) for chosen option
                switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).money_level_x_E_chosen
                    case 0
                    case 1 % (money level)*(effort level) chosen option
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': money x effort chosen'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: money x effort chosen '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    otherwise
                        error('not ready yet');
                end
                
                % (reward)*(effort) for high effort option
                switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).R_level_x_E_varOption
                    case 0
                    case 1 % (R level)*(effort level) high E option
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': R x E non-default'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: R x E non-default option '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    otherwise
                        error('not ready yet');
                end
                
                % (reward)*(effort) for chosen option
                switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).R_level_x_E_chosen
                    case 0
                    case 1 % (R level)*(effort level) chosen option
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': R x E chosen'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: R x E chosen '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    otherwise
                        error('not ready yet');
                end
                
                % (punishment)*(effort) for high effort option
                switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).P_level_x_E_varOption
                    case 0
                    case 1 % (P level)*(effort level) high E option
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': P x E non-default'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: P x E non-default option '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    otherwise
                        error('not ready yet');
                end
                
                % (punishment)*(effort) for chosen option
                switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).P_level_x_E_chosen
                    case 0
                    case 1 % (P level)*(effort level) chosen option
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': P x E chosen'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: P x E chosen '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    otherwise
                        error('not ready yet');
                end
                
                % net value chosen option
                switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).NV_chosen
                    case 0
                    case {1,3}
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': NVch-NVunch'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: NVch-NVunch '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2 % p(chosen)
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': p(chosen)'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: p(chosen) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
                
                % net value variable option
                switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).NV_varOption
                    case 0
                    case 1 % NV(high E) - NV(low E)
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': net value high E - low E'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: net value high E - low E '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2 % |NV(high E) - NV(low E)|
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': |net value high E - low E|'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: |net value high E - low E| '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 3 % p(choice = high E)
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': p(choice = high E)'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: p(choice = high E) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 4 % NV(high E) - NV(low E) + bias
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': net value high E - low E + bias'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: net value high E - low E + bias '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 5 % |NV(high E) - NV(low E) + bias|
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': |net value high E - low E + bias|'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: |net value high E - low E + bias| '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
                
                % net value variable option bis
                switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).NV_varOption_bis
                    case 0
                    case 1 % NV(high E) - NV(low E)
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': net value high E - low E'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: net value high E - low E '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2 % |NV(high E) - NV(low E)|
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': |net value high E - low E|'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: |net value high E - low E| '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 3 % p(choice = high E)
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': p(choice = high E)'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: p(choice = high E) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 4 % NV(high E) - NV(low E) + bias
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': net value high E - low E + bias'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: net value high E - low E + bias '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 5 % |NV(high E) - NV(low E) + bias|
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': |net value high E - low E + bias|'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: |net value high E - low E + bias| '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
                
                if strcmp(task_id_nm, 'Ep') % physical effort only
                    % area under the curve of the force that is gonna be produced
                    switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).F_integral
                        case 0
                        case {1,3}
                            n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                            reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': force integral'];
                            dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: effort integral '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                        case {2,4}
                            n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                            reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': force integral overshoot'];
                            dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: effort integral overshoot '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    end
                    
                    % fatigue
                    switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).fatigue
                        case 0
                        case 1
                            n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                            reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': fatigue'];
                            dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: fatigue '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    end
                    
                    % Ech*fatigue
                    switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).Ech_x_fatigue
                        case 0
                        case 1
                            n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                            reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': Ech_x_fatigue'];
                            dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: Ech_x_fatigue '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    end
                end % physical effort filter
                
                if strcmp(task_id_nm,'Em') % mental effort only
                    % efficacy of the next trial
                    switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).efficacy
                        case 0
                        case {1,2,3,4}
                            n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                            reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': efficacy'];
                            dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: efficacy '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    end
                    
                    % efficacy during the previous trial
                    switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).prevEfficacy
                        case 0
                        case {1,2,3,4}
                            n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                            reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': previous efficacy'];
                            dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: previous efficacy '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    end
                    
                    % (effort chosen)*(efficacy during the previous trial)
                    switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).Ech_x_prevEfficacy
                        case 0
                        case {1,2,3,4}
                            n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                            reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': (Ech)x(previous efficacy)'];
                            dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: (Ech)x(previous efficacy) '],dispRegs);
                            % if derivative added => add derivatives
                            n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    end
                end % mental effort filter
                
                % trial number
                switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).trialN
                        case 0
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': trial number'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: trial number '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': (trial number)x(Echosen-Eunchosen) (E levels)'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: (trial number)x(effort chosen-effort unchosen) (effort levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 3
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': (trial number)x(EnonDef-Edef) (E levels)'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: (trial number)x(effort non-default - effort default) (effort levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 4
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': (trial number)x(EnonDef) (E levels)'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: (trial number)x(effort non-default) (effort levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
                
                % confidence
                switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).confidence
                        case 0
                    case 1 % confidence rating by the subjects
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': confidence'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: confidence (levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case {2,3,4} % confidence inferred by the model
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': confidence'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: confidence (inferred by the model) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
                
                % RT (last regressor)
                switch GLMprm.chosen.(task_id_nm).(RP_dispChosen_nm).(splitE_dispChosen_nm).RT
                        case 0
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': RT'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: RT (raw) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': RT'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: RT (zscored per run) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 3
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG chosen ',RP_dispChosen_nm,' ',splitE_dispChosen_nm,': RT'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') chosen: RT (zscored per subject ie across all runs) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
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
                
                %% effort period onset
                n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['ONSET preEffort black cross ',RP_preEcross_nm,' ',splitE_preEcross_nm];
                dispRegFn([num2str(n_regs.(task_id_nm)),') ONSET preEffort black cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,': ',...
                    GLMprm.model_onset.(task_id_nm).preEffortCross,' '],dispRegs);
                % if derivative added => add derivatives
                n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                
                %% effort period regressors
                % binary variable indicating when choice = high effort option
                if GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).choiceHighE == 1
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                    reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,': choice = highE'];
                    dispRegFn([num2str(n_regs.(task_id_nm)),') pre-effort cross: choice hE '],dispRegs);
                    % if derivative added => add derivatives
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
                
                % money chosen
                switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).money_chosen
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,': money chosen'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') pre-effort cross: money chosen (amounts) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,': money chosen'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') pre-effort cross: |money chosen| (amounts) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 3
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,': money chosen'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') pre-effort cross: money chosen (levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 4
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,': money chosen'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') pre-effort cross: |money chosen| (levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
                
                % effort chosen
                switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).E_chosen
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,': effort chosen'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') pre-effort cross: effort chosen (levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
                
                % force peak
                switch task_id_nm
                    case 'Ep'
                        switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).F_peak
                            case {1,2}
                                n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                                reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,': force peak'];
                                dispRegFn([num2str(n_regs.(task_id_nm)),') pre-effort cross: force peak '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                        end
                        
                        % force integral
                        switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).F_integral
                            case {1,3}
                                n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                                reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,': force integral'];
                                dispRegFn([num2str(n_regs.(task_id_nm)),') pre-effort cross: force integral '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                            case {2,4}
                                error('please update');
                        end
                        
                    case 'Em'
                        % average RT for all numbers of each trial
                        switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).RT_avg
                            case 1
                                n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                                reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,': average RT effort'];
                                dispRegFn([num2str(n_regs.(task_id_nm)),') pre-effort cross: average RT effort '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                        end
                        
                        % number of correct answers
                        switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).n_correct
                            case 1
                                n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                                reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,': number of correct answers'];
                                dispRegFn([num2str(n_regs.(task_id_nm)),') pre-effort cross: number of correct answers '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                        end
                        
                        % number of errors
                        switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).n_errors
                            case 1
                                n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                                reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,': number of errors'];
                                dispRegFn([num2str(n_regs.(task_id_nm)),') pre-effort cross: number of errors '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                        end
                end
                
                % net value chosen option
                switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).NV_chosen
                    case 0
                    case {1,3}
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,': NVch-NVunch'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') pre-effort cross: NVch-NVunch '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG preEffort cross ',RP_preEcross_nm,' ',splitE_preEcross_nm,': p(chosen)'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') pre-effort cross: p(chosen) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    otherwise
                        error('problem with pre-E cross NV_chosen value');
                end
                
                % net value non-default option
                switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).NV_varOption
                    case 1 % NV(high E) - NV(low E)
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,': net value high E - low E'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') pre-effort cross: net value high E - low E '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2 % |NV(high E) - NV(low E)|
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,': |net value high E - low E|'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') pre-effort cross: |net value high E - low E| '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 3 % p(choice = high E)
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,': p(choice = high E)'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') pre-effort cross: p(choice = high E) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 4 % NV(high E) - NV(low E) + bias
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,': net value high E - low E + bias'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') pre-effort cross: net value high E - low E + bias '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 5 % |NV(high E) - NV(low E) + bias|
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,': |net value high E - low E + bias|'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') pre-effort cross: |net value high E - low E + bias| '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
                
                % net value non-default option bis
                switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).NV_varOption_bis
                    case 1 % NV(high E) - NV(low E)
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,': net value high E - low E'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') pre-effort cross: net value high E - low E '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2 % |NV(high E) - NV(low E)|
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,': |net value high E - low E|'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') pre-effort cross: |net value high E - low E| '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 3 % p(choice = high E)
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,': p(choice = high E)'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') pre-effort cross: p(choice = high E) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 4 % NV(high E) - NV(low E) + bias
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,': net value high E - low E + bias'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') pre-effort cross: net value high E - low E + bias '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 5 % |NV(high E) - NV(low E) + bias|
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,': |net value high E - low E + bias|'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') pre-effort cross: |net value high E - low E + bias| '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
                
                % RT first answer
                switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).RT_1stAnswer
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,': RT 1st answer'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') pre-effort cross: RT 1st answer (raw) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
                
                % trial number
                switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).trialN
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,': trial number'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') pre-effort cross: trial number '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,': (trial number)x(Echosen-Eunchosen) (E levels)'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') pre-effort cross: (trial number)x(effort chosen-effort unchosen) (effort levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 3
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,': (trial number)x(EnonDef-Edef) (E levels)'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') pre-effort cross: (trial number)x(effort non-default - effort default) (effort levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 4
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,': (trial number)x(EnonDef) (E levels)'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') pre-effort cross: (trial number)x(effort non-default) (effort levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
                
                % confidence
                switch GLMprm.preEffortCross.(task_id_nm).(RP_preEcross_nm).(splitE_preEcross_nm).confidence
                    case 1 % confidence ratings 0/1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,': confidence'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') pre-effort cross: confidence (levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case {2,3,4} % confidence inferred by the model
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_preEcross_nm,' ',splitE_preEcross_nm,': confidence'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') pre-effort cross: confidence (inferred by the model) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
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
                
                %% effort period onset
                n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['ONSET effort ',RP_Eperf_nm,' ',splitE_Eperf_nm];
                dispRegFn([num2str(n_regs.(task_id_nm)),') ONSET effort period ',RP_Eperf_nm,' ',splitE_Eperf_nm,': ',...
                    GLMprm.model_onset.(task_id_nm).Eperf,' '],dispRegs);
                % if derivative added => add derivatives
                n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                
                %% effort period regressors
                % binary variable indicating when choice = high effort option
                if GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).choiceHighE == 1
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                    reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': choice = highE'];
                    dispRegFn([num2str(n_regs.(task_id_nm)),') effort period: choice hE '],dispRegs);
                    % if derivative added => add derivatives
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
                
                % money chosen
                switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).money_chosen
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': money chosen'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') effort period: money chosen (amounts) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': money chosen'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') effort period: |money chosen| (amounts) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 3
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': money chosen'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') effort period: money chosen (levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 4
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': money chosen'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') effort period: |money chosen| (levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
                
                % effort chosen
                switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).E_chosen
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': effort chosen'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') effort period: effort chosen (levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
                
                % force peak
                switch task_id_nm
                    case 'Ep'
                        % force peak
                        switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).F_peak
                            case {1,2}
                                n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                                reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': force peak'];
                                dispRegFn([num2str(n_regs.(task_id_nm)),') effort period: force peak '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                        end
                        
                        % force integral
                        switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).F_integral
                            case {1,3}
                                n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                                reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': force integral'];
                                dispRegFn([num2str(n_regs.(task_id_nm)),') effort period: force integral '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                            case {2,4}
                                n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                                reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': force integral overshoot'];
                                dispRegFn([num2str(n_regs.(task_id_nm)),') effort period: force integral overshoot '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                        end
                        
                    case 'Em'
                        % efficacy
                        switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).efficacy
                            case {1,2,3,4}
                                n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                                reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': efficacy'];
                                dispRegFn([num2str(n_regs.(task_id_nm)),') effort period: efficacy '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                        end
                        
                        % average RT for all numbers of each trial
                        switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).RT_avg
                            case 1
                                n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                                reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': average RT effort'];
                                dispRegFn([num2str(n_regs.(task_id_nm)),') effort period: average RT effort '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                        end
                        
                        % number of correct answers
                        switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).n_correct
                            case 1
                                n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                                reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': number of correct answers'];
                                dispRegFn([num2str(n_regs.(task_id_nm)),') effort period: number of correct answers '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                        end
                        
                        % number of errors
                        switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).n_errors
                            case 1
                                n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                                reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': number of errors'];
                                dispRegFn([num2str(n_regs.(task_id_nm)),') effort period: number of errors '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                        end
                end
                
                % net value chosen option
                switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).NV_chosen
                    case {1,3}
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': NVch-NVunch'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') effort period: NVch-NVunch '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2 % p(chosen)
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': p(chosen)'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') effort period: p(chosen) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
                
                % net value non-default option
                switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).NV_varOption
                    case 1 % NV(high E) - NV(low E)
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': net value high E - low E'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') effort period: net value high E - low E '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2 % |NV(high E) - NV(low E)|
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': |net value high E - low E|'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') effort period: |net value high E - low E| '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 3 % p(choice = high E)
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': p(choice = high E)'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') effort period: p(choice = high E) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 4 % NV(high E) - NV(low E) + bias
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': net value high E - low E + bias'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') effort period: net value high E - low E + bias '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 5 % |NV(high E) - NV(low E) + bias|
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': |net value high E - low E + bias|'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') effort period: |net value high E - low E + bias| '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
                
                % net value non-default option bis
                switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).NV_varOption_bis
                    case 1 % NV(high E) - NV(low E)
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': net value high E - low E'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') effort period: net value high E - low E '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2 % |NV(high E) - NV(low E)|
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': |net value high E - low E|'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') effort period: |net value high E - low E| '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 3 % p(choice = high E)
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': p(choice = high E)'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') effort period: p(choice = high E) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 4 % NV(high E) - NV(low E) + bias
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': net value high E - low E + bias'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') effort period: net value high E - low E + bias '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 5 % |NV(high E) - NV(low E) + bias|
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': |net value high E - low E + bias|'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') effort period: |net value high E - low E + bias| '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
                
                % RT first answer
                switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).RT_1stAnswer
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': effort latency'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') effort period: effort latency (raw) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
                
                % fatigue
                switch task_id_nm
                    case 'Ep'
                        % fatigue
                        switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).fatigue
                            case 1
                                n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                                reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': fatigue'];
                                dispRegFn([num2str(n_regs.(task_id_nm)),') effort period: fatigue '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                        end
                        
                        % Ech*fatigue
                        switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).Ech_x_fatigue
                            case 1
                                n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                                reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': Ech_x_fatigue'];
                                dispRegFn([num2str(n_regs.(task_id_nm)),') effort period: Ech_x_fatigue '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                        end
                        
                    case 'Em'
                        % efficacy during the previous trial (learning boost)
                        switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).prevEfficacy
                            case {1,2,3,4}
                                n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                                reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': previous efficacy'];
                                dispRegFn([num2str(n_regs.(task_id_nm)),') effort period: previous efficacy '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                        end
                        
                        % (effort chosen)*(efficacy during the previous trial)
                        switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).Ech_x_prevEfficacy
                            case {1,2,3,4}
                                n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                                reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': (Ech)x(previous efficacy)'];
                                dispRegFn([num2str(n_regs.(task_id_nm)),') effort period: (Ech)x(previous efficacy) '],dispRegs);
                                % if derivative added => add derivatives
                                n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                        end
                end
                
                % trial number
                switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).trialN
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': trial number'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') effort period: trial number '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': (trial number)x(Echosen-Eunchosen) (E levels)'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') effort period: (trial number)x(effort chosen-effort unchosen) (effort levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 3
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': (trial number)x(EnonDef-Edef) (E levels)'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') effort period: (trial number)x(effort non-default - effort default) (effort levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 4
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': (trial number)x(EnonDef) (E levels)'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') effort period: (trial number)x(effort non-default) (effort levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
                
                % confidence
                switch GLMprm.Eperf.(task_id_nm).(RP_Eperf_nm).(splitE_Eperf_nm).confidence
                    case 1 % confidence ratings 0/1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': confidence'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') effort period: confidence (levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case {2,3,4} % confidence inferred by the model
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG effort ',RP_Eperf_nm,' ',splitE_Eperf_nm,': confidence'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') effort period: confidence (inferred by the model) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
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
                
                %% feedback onset
                n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['ONSET feedback ',RP_fbk_nm,' ',splitE_fbk_nm];
                dispRegFn([num2str(n_regs.(task_id_nm)),') ONSET feedback: ',GLMprm.model_onset.(task_id_nm).fbk,' '],dispRegs);
                % if derivative added => add derivatives
                n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                
                %% feedback regressors
                % win vs loss
                switch GLMprm.fbk.(task_id_nm).(RP_fbk_nm).(splitE_fbk_nm).win_vs_loss
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG feedback ',RP_fbk_nm,' ',splitE_fbk_nm,': win-loss'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') feedback: win-loss '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
                
                % binary variable indicating when choice = high effort option
                if GLMprm.fbk.(task_id_nm).(RP_fbk_nm).(splitE_fbk_nm).choiceHighE == 1
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                    reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG feedback ',RP_fbk_nm,' ',splitE_fbk_nm,': choice = highE'];
                    dispRegFn([num2str(n_regs.(task_id_nm)),') feedback: choice hE '],dispRegs);
                    % if derivative added => add derivatives
                    n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
                
                % money obtained
                switch GLMprm.fbk.(task_id_nm).(RP_fbk_nm).(splitE_fbk_nm).money_obtained
                    case 1 % money amount
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG feedback ',RP_fbk_nm,' ',splitE_fbk_nm,': money obtained'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') feedback: money obtained (amount) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2 % |money amount|
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG feedback ',RP_fbk_nm,' ',splitE_fbk_nm,': |money obtained|'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') feedback: |money obtained| (amount) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
                
                % effort performed
                switch GLMprm.fbk.(task_id_nm).(RP_fbk_nm).(splitE_fbk_nm).E_made
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG feedback ',RP_fbk_nm,' ',splitE_fbk_nm,': effort performed'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') feedback: effort performed (levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
                
                % trial number
                switch GLMprm.fbk.(task_id_nm).(RP_fbk_nm).(splitE_fbk_nm).trialN
                    case 1
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG feedback ',RP_fbk_nm,' ',splitE_fbk_nm,': trial number'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') feedback: trial number '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 2
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG feedback ',RP_fbk_nm,' ',splitE_fbk_nm,': trial number x (Echosen-Eunchosen) (E levels)'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') feedback: trial number x (effort chosen-effort unchosen) (effort levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case 3
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG feedback ',RP_fbk_nm,' ',splitE_fbk_nm,': trial number x (EnonDef-Edef) (E levels)'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') feedback: trial number x (effort non-default - effort default) (effort levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
                
                % confidence
                switch GLMprm.fbk.(task_id_nm).(RP_fbk_nm).(splitE_fbk_nm).confidence
                    case 1 % confidence rated by the subjects
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG feedback ',RP_fbk_nm,' ',splitE_fbk_nm,': confidence'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') feedback: confidence (levels) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                    case {2,3,4} % confidence inferred by the model
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
                        reg_names.(task_id_nm){n_regs.(task_id_nm)} = ['REG feedback ',RP_fbk_nm,' ',splitE_fbk_nm,': confidence'];
                        dispRegFn([num2str(n_regs.(task_id_nm)),') feedback: confidence (inferred by the model) '],dispRegs);
                        % if derivative added => add derivatives
                        n_regs.(task_id_nm) = n_regs.(task_id_nm) + add_drv;
                end
            end % loop effort levels
        end % loop reward/punishment
    end % choice onset
    
    %% movement regressors
    n_mvmt = 6;
    for iMvmt = 1:n_mvmt
        n_regs.(task_id_nm) = n_regs.(task_id_nm) + 1;
        reg_names.(task_id_nm){n_regs.(task_id_nm)} = 'movement';
    end
    dispRegFn([num2str(n_regs.(task_id_nm)-n_mvmt+1),'-',num2str(n_regs.(task_id_nm)),') movement regressors '], dispRegs);
    % note: if temporal derivative has been added, movement regressors
    % are not doubled (ie number stays stable independent of add_drv value)
    
end % physical/mental loop
end % function

function[] = dispRegFn(text, dispRegs)
switch dispRegs
    case 1
        disp(text);
end
end