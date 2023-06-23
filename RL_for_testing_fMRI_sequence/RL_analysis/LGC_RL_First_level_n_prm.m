function [ n_prm, prm_idx ] = LGC_RL_First_level_n_prm( GLMprm )
%[ n_prm, prm_idx ] = LGC_RL_First_level_n_prm( GLMprm )
% estimate number and identity of regressors for each run for each task for
% each subject
%
% INPUTS
% GLMprm: structure with GLM parameters
%
% OUTPUTS
% n_prm: structure containing the number of parameters for each run for
% each task for the current subject
%
% prm_idx: structure with corresponding index for each regressor
%

%% extract main relevant parameters for the script to work

% need to know add derivative to know number of betas/regressor
add_drv = GLMprm.gal.add_drv;
% infer number of parameters/regressor accordingly
switch add_drv
    case 0
        nBperR = 1;
    case 1 % normal + temporal derivative
        nBperR = 2;
    case 2 % normal + temporal derivative + spatial derivative
        nBperR = 3;
end

% reinforcement-learning task
RLprm = GLMprm;

% total number of runs
n_total_RL_runs     = 3;
n_total_runs = n_total_RL_runs;

%% movement parameters
n_mvmt_perRun = 6*ones(1,n_total_runs);

%% parameters of interest for output
n_prm.all = 0 ;
task_names = {'RL'};
n_tasks = length(task_names);
task_nm = task_names{1};

% set initial number of parameters per task to zero
for iRun = 1:n_total_runs
    run_nm = ['run',num2str(iRun)];
    n_prm.(task_nm).(run_nm) = 0;
end

%% check which onsets and modulators are used in the current GLM and assign
% them the corresponding index in the list of parameters
bf_idx = 0; % number of regressors before the current one
prm_idx = struct;

for iRun = 1:n_total_runs
    run_nm = ['run',num2str(iRun)];
    
    
    %% RL
    task_nm = 'RL';
    task_id = 'L_';
    
    %% regrouping all sessions together
        %% stim display
    o_stim = RLprm.o_stim;
    mod_stim = RLprm.mod_stim;
    switch o_stim
        case 1 % all pairs pooled
            % onset stim
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_stim'], bf_idx);
            
            % modulators
            % RT
            if mod_stim.RT ~= 0 && ismember(mod_stim.RT,9:16)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT'], bf_idx);
            end
            % trial number
            if ismember(mod_stim.trialN,[1,2,3])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN'], bf_idx);
            end
            
            % SV = pA*QA+pB*QB
            if ismember(mod_stim.SV,1:5)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_SV'], bf_idx);
            end
            
            % dQ/dQ
            if mod_stim.dQ ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_dQ'], bf_idx);
            end
            
            % p(choice=best option)
            if mod_stim.pBest ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_pBest'], bf_idx);
            end
            
            % trial number
            if ismember(mod_stim.trialN,[4,5,6])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN'], bf_idx);
            end
            
            % ROI activity
            if mod_stim.ROI_activity_yn ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_GLM',mod_stim.ROI_activity_GLM,...
                    '_',mod_stim.ROI_activity_period,'_period_',mod_stim.ROI_activity_ROI_nm,'_activity'], bf_idx);
            end
            
            % RT
            if mod_stim.RT ~= 0 && ismember(mod_stim.RT,1:8)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT'], bf_idx);
            end
            
            % SV = pA*QA+pB*QB
            if ismember(mod_stim.SV,6:10)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_SV'], bf_idx);
            end
            
            % trial number
            if ismember(mod_stim.trialN,[7,8,9])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN'], bf_idx);
            end
            
        case 2 % separate depending on pair type gain/neutral/loss
            
            %% gain pair
            trialType_nm = '_gainPair';
            % onset stim
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_stim',trialType_nm], bf_idx);
            
            % modulators
            % RT
            if mod_stim.RT ~= 0 && ismember(mod_stim.RT,9:16)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_stim.trialN,[1,2,3])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx);
            end
            
            % SV = pA*QA+pB*QB
            if ismember(mod_stim.SV,1:5)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_SV',trialType_nm], bf_idx);
            end
            
            % dQ
            if mod_stim.dQ ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_dQ',trialType_nm], bf_idx);
            end
            
            % p(choice=best option)
            if mod_stim.pBest ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_pBest',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_stim.trialN,[4,5,6])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx);
            end
            
            % ROI activity
            if mod_stim.ROI_activity_yn ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_GLM',mod_stim.ROI_activity_GLM,...
                    '_',mod_stim.ROI_activity_period,'_period_',mod_stim.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx);
            end
            
            % RT
            if mod_stim.RT ~= 0 && ismember(mod_stim.RT,1:8)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx);
            end
            
            % SV = pA*QA+pB*QB
            if ismember(mod_stim.SV,6:10)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_SV',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_stim.trialN,[7,8,9])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx);
            end
            
            %% neutral pair
            trialType_nm = '_ntalPair';
            % onset stim
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_stim',trialType_nm], bf_idx);
            
            % modulators
            % RT
            if mod_stim.RT ~= 0 && ismember(mod_stim.RT,9:16)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_stim.trialN,[1,2,3])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx);
            end
            
            % no SV or dQ for neutral pair (always equal to zero)
            
            % trial number
            if ismember(mod_stim.trialN,[4,5,6])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx);
            end
            
            % ROI activity
            if mod_stim.ROI_activity_yn ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_GLM',mod_stim.ROI_activity_GLM,...
                    '_',mod_stim.ROI_activity_period,'_period_',mod_stim.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx);
            end
            
            % RT
            if mod_stim.RT ~= 0 && ismember(mod_stim.RT,1:8)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_stim.trialN,[7,8,9])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx);
            end
            
            %% loss pair
            trialType_nm = '_lossPair';
            % onset stim
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_stim',trialType_nm], bf_idx);
            
            % modulators
            % RT
            if mod_stim.RT ~= 0 && ismember(mod_stim.RT,9:16)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx);
            end
            % trial number
            if ismember(mod_stim.trialN,[1,2,3])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx);
            end
            
            % SV = pA*QA+pB*QB
            if ismember(mod_stim.SV,1:5)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_SV',trialType_nm], bf_idx);
            end
            
            % dQ
            if mod_stim.dQ ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_dQ',trialType_nm], bf_idx);
            end
            
            % p(choice = best option)
            if mod_stim.pBest ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_pBest',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_stim.trialN,[4,5,6])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx);
            end
            
            % ROI activity
            if mod_stim.ROI_activity_yn ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_GLM',mod_stim.ROI_activity_GLM,...
                    '_',mod_stim.ROI_activity_period,'_period_',mod_stim.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx);
            end
            
            % RT
            if mod_stim.RT ~= 0 && ismember(mod_stim.RT,1:8)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx);
            end
            
            % SV = pA*QA+pB*QB
            if ismember(mod_stim.SV,6:10)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_SV',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_stim.trialN,[7,8,9])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx);
            end
            
        case 3 % gain + loss pooled, neutral apart
            %% gain + loss pair
            trialType_nm = '_GL_Pairs';
            % onset stim
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_stim',trialType_nm], bf_idx);
            
            % modulators
            % RT
            if mod_stim.RT ~= 0 && ismember(mod_stim.RT,9:16)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_stim.trialN,[1,2,3])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx);
            end
            
            % SV = pA*QA+pB*QB
            if ismember(mod_stim.SV,1:5)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_SV',trialType_nm], bf_idx);
            end
            
            % dQ
            if mod_stim.dQ ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_dQ',trialType_nm], bf_idx);
            end
            
            % p(choice=best option)
            if mod_stim.pBest ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_pBest',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_stim.trialN,[4,5,6])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx);
            end
            
            % ROI activity
            if mod_stim.ROI_activity_yn ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_GLM',mod_stim.ROI_activity_GLM,...
                    '_',mod_stim.ROI_activity_period,'_period_',mod_stim.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx);
            end
            
            % RT
            if mod_stim.RT ~= 0 && ismember(mod_stim.RT,1:8)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx);
            end
            
            % SV = pA*QA+pB*QB
            if ismember(mod_stim.SV,6:10)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_SV',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_stim.trialN,[7,8,9])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx);
            end
            
            %% neutral pair
            trialType_nm = '_ntalPair';
            % onset stim
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_stim',trialType_nm], bf_idx);
            
            % modulators
            % RT
            if mod_stim.RT ~= 0 && ismember(mod_stim.RT,9:16)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_stim.trialN,[1,2,3])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx);
            end
            
            % no SV or dQ for neutral pair (always equal to zero)
            
            % trial number
            if ismember(mod_stim.trialN,[4,5,6])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx);
            end
            
            % ROI activity
            if mod_stim.ROI_activity_yn ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_GLM',mod_stim.ROI_activity_GLM,...
                    '_',mod_stim.ROI_activity_period,'_period_',mod_stim.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx);
            end
            
            % RT
            if mod_stim.RT ~= 0 && ismember(mod_stim.RT,1:8)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_stim.trialN,[7,8,9])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx);
            end
        case 4
            error('not ready yet');
        case 5 % split gain/neutral/loss and first/second half
            %% gain pair first trials
            trialType_nm = '_gainPair_first';
            % onset stim
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_stim',trialType_nm], bf_idx);
            
            % modulators
            % RT
            if mod_stim.RT ~= 0 && ismember(mod_stim.RT,9:16)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_stim.trialN,[1,2,3])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx);
            end
            
            % SV = pA*QA+pB*QB
            if ismember(mod_stim.SV,1:5)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_SV',trialType_nm], bf_idx);
            end
            
            % dQ
            if mod_stim.dQ ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_dQ',trialType_nm], bf_idx);
            end
            
            % p(choice=best option)
            if mod_stim.pBest ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_pBest',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_stim.trialN,[4,5,6])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx);
            end
            
            % ROI activity
            if mod_stim.ROI_activity_yn ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_GLM',mod_stim.ROI_activity_GLM,...
                    '_',mod_stim.ROI_activity_period,'_period_',mod_stim.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx);
            end
            
            % RT
            if mod_stim.RT ~= 0 && ismember(mod_stim.RT,1:8)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx);
            end
            
            % SV = pA*QA+pB*QB
            if ismember(mod_stim.SV,6:10)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_SV',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_stim.trialN,[7,8,9])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx);
            end
            
            %% neutral pair first trials
            trialType_nm = '_ntalPair_first';
            % onset stim
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_stim',trialType_nm], bf_idx);
            
            % modulators
            % RT
            if mod_stim.RT ~= 0 && ismember(mod_stim.RT,9:16)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_stim.trialN,[1,2,3])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx);
            end
            
            % no SV, no dQ for neutral pair
            
            % trial number
            if ismember(mod_stim.trialN,[4,5,6])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx);
            end
            
            % ROI activity
            if mod_stim.ROI_activity_yn ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_GLM',mod_stim.ROI_activity_GLM,...
                    '_',mod_stim.ROI_activity_period,'_period_',mod_stim.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx);
            end
            
            % RT
            if mod_stim.RT ~= 0 && ismember(mod_stim.RT,1:8)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_stim.trialN,[7,8,9])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx);
            end
            
            
            %% loss pair first trials
            trialType_nm = '_lossPair_first';
            % onset stim
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_stim',trialType_nm], bf_idx);
            
            % modulators
            % RT
            if mod_stim.RT ~= 0 && ismember(mod_stim.RT,9:16)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_stim.trialN,[1,2,3])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx);
            end
            
            % SV = pA*QA+pB*QB
            if ismember(mod_stim.SV,1:5)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_SV',trialType_nm], bf_idx);
            end
            
            % dQ
            if mod_stim.dQ ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_dQ',trialType_nm], bf_idx);
            end
            
            % p(choice = best option)
            if mod_stim.pBest ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_pBest',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_stim.trialN,[4,5,6])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx);
            end
            
            % ROI activity
            if mod_stim.ROI_activity_yn ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_GLM',mod_stim.ROI_activity_GLM,...
                    '_',mod_stim.ROI_activity_period,'_period_',mod_stim.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx);
            end
            
            % RT
            if mod_stim.RT ~= 0 && ismember(mod_stim.RT,1:8)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx);
            end
            
            % SV = pA*QA+pB*QB
            if ismember(mod_stim.SV,6:10)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_SV',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_stim.trialN,[7,8,9])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx);
            end
            
            %% gain pair last trials
            trialType_nm = '_gainPair_last';
            % onset stim
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_stim',trialType_nm], bf_idx);
            
            % modulators
            % RT
            if mod_stim.RT ~= 0 && ismember(mod_stim.RT,9:16)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_stim.trialN,[1,2,3])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx);
            end
            
            % SV = pA*QA+pB*QB
            if ismember(mod_stim.SV,1:5)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_SV',trialType_nm], bf_idx);
            end
            
            % dQ
            if mod_stim.dQ ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_dQ',trialType_nm], bf_idx);
            end
            
            % p(choice=best option)
            if mod_stim.pBest ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_pBest',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_stim.trialN,[4,5,6])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx);
            end
            
            % ROI activity
            if mod_stim.ROI_activity_yn ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_GLM',mod_stim.ROI_activity_GLM,...
                    '_',mod_stim.ROI_activity_period,'_period_',mod_stim.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx);
            end
            
            % RT
            if mod_stim.RT ~= 0 && ismember(mod_stim.RT,1:8)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx);
            end
            
            % SV = pA*QA+pB*QB
            if ismember(mod_stim.SV,6:10)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_SV',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_stim.trialN,[7,8,9])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx);
            end
            
            %% neutral pair last trials
            trialType_nm = '_ntalPair_last';
            % onset stim
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_stim',trialType_nm], bf_idx);
            
            % modulators
            % RT
            if mod_stim.RT ~= 0 && ismember(mod_stim.RT,9:16)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_stim.trialN,[1,2,3])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx);
            end
            
            % no SV, no dQ for neutral pair
            
            % trial number
            if ismember(mod_stim.trialN,[4,5,6])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx);
            end
            
            % ROI activity
            if mod_stim.ROI_activity_yn ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_GLM',mod_stim.ROI_activity_GLM,...
                    '_',mod_stim.ROI_activity_period,'_period_',mod_stim.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx);
            end
            
            % RT
            if mod_stim.RT ~= 0 && ismember(mod_stim.RT,1:8)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_stim.trialN,[7,8,9])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx);
            end
            
            %% loss pair last trials
            trialType_nm = '_lossPair_last';
            % onset stim
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_stim',trialType_nm], bf_idx);
            
            % modulators
            % RT
            if mod_stim.RT ~= 0 && ismember(mod_stim.RT,9:16)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_stim.trialN,[1,2,3])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx);
            end
            
            % SV = pA*QA+pB*QB
            if ismember(mod_stim.SV,1:5)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_SV',trialType_nm], bf_idx);
            end
            
            % dQ
            if mod_stim.dQ ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_dQ',trialType_nm], bf_idx);
            end
            
            % p(choice = best option)
            if mod_stim.pBest ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_pBest',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_stim.trialN,[4,5,6])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx);
            end
            
            % ROI activity
            if mod_stim.ROI_activity_yn ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_GLM',mod_stim.ROI_activity_GLM,...
                    '_',mod_stim.ROI_activity_period,'_period_',mod_stim.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx);
            end
            
            % RT
            if mod_stim.RT ~= 0 && ismember(mod_stim.RT,1:8)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx);
            end
            
            % SV = pA*QA+pB*QB
            if ismember(mod_stim.SV,6:10)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_SV',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_stim.trialN,[7,8,9])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx);
            end
        case 6 % pool gain+loss/neutral apart and split first/second half trials also
            %% gain + loss pair - first trials
            trialType_nm = '_GL_Pairs_first';
            % onset stim
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_stim',trialType_nm], bf_idx);
            
            % modulators
            % RT
            if mod_stim.RT ~= 0 && ismember(mod_stim.RT,9:16)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_stim.trialN,[1,2,3])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx);
            end
            
            % SV = pA*QA+pB*QB
            if ismember(mod_stim.SV,1:5)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_SV',trialType_nm], bf_idx);
            end
            
            % dQ
            if mod_stim.dQ ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_dQ',trialType_nm], bf_idx);
            end
            
            % p(choice=best option)
            if mod_stim.pBest ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_pBest',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_stim.trialN,[4,5,6])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx);
            end
            
            % ROI activity
            if mod_stim.ROI_activity_yn ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_GLM',mod_stim.ROI_activity_GLM,...
                    '_',mod_stim.ROI_activity_period,'_period_',mod_stim.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx);
            end
            
            % RT
            if mod_stim.RT ~= 0 && ismember(mod_stim.RT,1:8)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx);
            end
            
            % SV = pA*QA+pB*QB
            if ismember(mod_stim.SV,6:10)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_SV',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_stim.trialN,[7,8,9])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx);
            end
            
            %% neutral pair - first trials
            trialType_nm = '_ntalPair_first';
            % onset stim
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_stim',trialType_nm], bf_idx);
            
            % modulators
            % RT
            if mod_stim.RT ~= 0 && ismember(mod_stim.RT,9:16)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_stim.trialN,[1,2,3])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx);
            end
            
            % no SV or dQ for neutral pair (always equal to zero)
            
            % trial number
            if ismember(mod_stim.trialN,[4,5,6])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx);
            end
            
            % ROI activity
            if mod_stim.ROI_activity_yn ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_GLM',mod_stim.ROI_activity_GLM,...
                    '_',mod_stim.ROI_activity_period,'_period_',mod_stim.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx);
            end
            
            % RT
            if mod_stim.RT ~= 0 && ismember(mod_stim.RT,1:8)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_stim.trialN,[7,8,9])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx);
            end
            
            %% gain + loss pair - last trials
            trialType_nm = '_GL_Pairs_last';
            % onset stim
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_stim',trialType_nm], bf_idx);
            
            % modulators
            % RT
            if mod_stim.RT ~= 0 && ismember(mod_stim.RT,9:16)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_stim.trialN,[1,2,3])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx);
            end
            
            % SV = pA*QA+pB*QB
            if ismember(mod_stim.SV,1:5)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_SV',trialType_nm], bf_idx);
            end
            
            % dQ
            if mod_stim.dQ ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_dQ',trialType_nm], bf_idx);
            end
            
            % p(choice=best option)
            if mod_stim.pBest ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_pBest',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_stim.trialN,[4,5,6])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx);
            end
            
            % ROI activity
            if mod_stim.ROI_activity_yn ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_GLM',mod_stim.ROI_activity_GLM,...
                    '_',mod_stim.ROI_activity_period,'_period_',mod_stim.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx);
            end
            
            % RT
            if mod_stim.RT ~= 0 && ismember(mod_stim.RT,1:8)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx);
            end
            
            % SV = pA*QA+pB*QB
            if ismember(mod_stim.SV,6:10)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_SV',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_stim.trialN,[7,8,9])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx);
            end
            
            %% neutral pair - last trials
            trialType_nm = '_ntalPair_last';
            % onset stim
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_stim',trialType_nm], bf_idx);
            
            % modulators
            % RT
            if mod_stim.RT ~= 0 && ismember(mod_stim.RT,9:16)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_stim.trialN,[1,2,3])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx);
            end
            
            % no SV or dQ for neutral pair (always equal to zero)
            
            % trial number
            if ismember(mod_stim.trialN,[4,5,6])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx);
            end
            
            % ROI activity
            if mod_stim.ROI_activity_yn ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_GLM',mod_stim.ROI_activity_GLM,...
                    '_',mod_stim.ROI_activity_period,'_period_',mod_stim.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx);
            end
            
            % RT
            if mod_stim.RT ~= 0 && ismember(mod_stim.RT,1:8)
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_RT',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_stim.trialN,[7,8,9])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_stim_trialN',trialType_nm], bf_idx);
            end
    end % stimulus onset
    
    %% answer
    o_answer = RLprm.o_answer;
    %         mod_answer = RLprm.mod_answer;
    switch o_answer
        case 1
            % onset answer
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_answer'], bf_idx);
            
        case 2 % separate depending on pair type
            
            %% gain pair
            trialType_nm = '_gainPair';
            % onset answer
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_answer',trialType_nm], bf_idx);
            
            %% neutral pair
            trialType_nm = '_ntalPair';
            % onset answer
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_answer',trialType_nm], bf_idx);
            
            %% loss pair
            trialType_nm = '_lossPair';
            % onset answer
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_answer',trialType_nm], bf_idx);
            
    end % answer onset
    
    %% chosen option in red
    o_chosen = RLprm.o_chosen;
    mod_chosen = RLprm.mod_chosen;
    switch o_chosen
        case 1
            % onset stim
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_chosen'], bf_idx);
            
            % modulators
            % trial number
            if ismember(mod_chosen.trialN,[1,2,3])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN'], bf_idx);
            end
            
            % SV = pA*QA+pB*QB
            if mod_chosen.SV ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_SV'], bf_idx);
            end
            
            % dQ/dQ
            if mod_chosen.dQ ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_dQ'], bf_idx);
            end
            
            % p(choice=best option)
            if mod_chosen.pBest ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_pBest'], bf_idx);
            end
            
            % trial number
            if ismember(mod_chosen.trialN,[4,5,6])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN'], bf_idx);
            end
            
            % ROI activity
            if mod_chosen.ROI_activity_yn ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_GLM',mod_chosen.ROI_activity_GLM,...
                    '_',mod_chosen.ROI_activity_period,'_period_',mod_chosen.ROI_activity_ROI_nm,'_activity'], bf_idx);
            end
            
            % RT
            if mod_chosen.RT ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_RT'], bf_idx);
            end
            
            % trial number
            if ismember(mod_chosen.trialN,[7,8,9])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN'], bf_idx);
            end
            
        case 2 % separate depending on pair type
            
            %% gain pair
            trialType_nm = '_gainPair';
            % onset stim
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_chosen',trialType_nm], bf_idx);
            
            % modulators
            % trial number
            if ismember(mod_chosen.trialN,[1,2,3])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx);
            end
            
            % SV = pA*QA+pB*QB
            if mod_chosen.SV ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_SV',trialType_nm], bf_idx);
            end
            
            % dQ
            if mod_chosen.dQ ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_dQ',trialType_nm], bf_idx);
            end
            
            % p(choice=best option)
            if mod_chosen.pBest ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_pBest',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_chosen.trialN,[4,5,6])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx);
            end
            
            % ROI activity
            if mod_chosen.ROI_activity_yn ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_GLM',mod_chosen.ROI_activity_GLM,...
                    '_',mod_chosen.ROI_activity_period,'_period_',mod_chosen.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx);
            end
            
            % RT
            if mod_chosen.RT ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_RT',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_chosen.trialN,[7,8,9])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx);
            end
            
            %% neutral pair
            trialType_nm = '_ntalPair';
            % onset stim
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_chosen',trialType_nm], bf_idx);
            
            % modulators
            % trial number
            if ismember(mod_chosen.trialN,[1,2,3])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx);
            end
            
            % no SV or dQ for neutral pair (always equal to zero)
            
            % trial number
            if ismember(mod_chosen.trialN,[4,5,6])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx);
            end
            
            % ROI activity
            if mod_chosen.ROI_activity_yn ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_GLM',mod_chosen.ROI_activity_GLM,...
                    '_',mod_chosen.ROI_activity_period,'_period_',mod_chosen.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx);
            end
            
            % RT
            if mod_chosen.RT ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_RT',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_chosen.trialN,[7,8,9])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx);
            end
            
            %% loss pair
            trialType_nm = '_lossPair';
            % onset stim
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_chosen',trialType_nm], bf_idx);
            
            % modulators
            % trial number
            if ismember(mod_chosen.trialN,[1,2,3])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx);
            end
            
            % SV = pA*QA+pB*QB
            if mod_chosen.SV ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_SV',trialType_nm], bf_idx);
            end
            
            % dQ
            if mod_chosen.dQ ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_dQ',trialType_nm], bf_idx);
            end
            
            % p(choice = best option)
            if mod_chosen.pBest ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_pBest',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_chosen.trialN,[4,5,6])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx);
            end
            
            % ROI activity
            if mod_chosen.ROI_activity_yn ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_GLM',mod_chosen.ROI_activity_GLM,...
                    '_',mod_chosen.ROI_activity_period,'_period_',mod_chosen.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx);
            end
            
            % RT
            if mod_chosen.RT ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_RT',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_chosen.trialN,[7,8,9])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx);
            end
            
        case 3
            %% gain + loss pair
            trialType_nm = '_GL_Pairs';
            % onset stim
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_chosen',trialType_nm], bf_idx);
            
            % modulators
            % trial number
            if ismember(mod_chosen.trialN,[1,2,3])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx);
            end
            
            % SV = pA*QA+pB*QB
            if mod_chosen.SV ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_SV',trialType_nm], bf_idx);
            end
            
            % dQ
            if mod_chosen.dQ ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_dQ',trialType_nm], bf_idx);
            end
            
            % p(choice=best option)
            if mod_chosen.pBest ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_pBest',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_chosen.trialN,[4,5,6])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx);
            end
            
            % ROI activity
            if mod_chosen.ROI_activity_yn ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_GLM',mod_chosen.ROI_activity_GLM,...
                    '_',mod_chosen.ROI_activity_period,'_period_',mod_chosen.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx);
            end
            
            % RT
            if mod_chosen.RT ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_RT',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_chosen.trialN,[7,8,9])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx);
            end
            
            %% neutral pair
            trialType_nm = '_ntalPair';
            % onset stim
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_chosen',trialType_nm], bf_idx);
            
            % modulators
            % trial number
            if ismember(mod_chosen.trialN,[1,2,3])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx);
            end
            
            % no SV or dQ for neutral pair (always equal to zero)
            
            % trial number
            if ismember(mod_chosen.trialN,[4,5,6])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx);
            end
            
            % ROI activity
            if mod_chosen.ROI_activity_yn ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_GLM',mod_chosen.ROI_activity_GLM,...
                    '_',mod_chosen.ROI_activity_period,'_period_',mod_chosen.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx);
            end
            
            % RT
            if mod_chosen.RT ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_RT',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_chosen.trialN,[7,8,9])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx);
            end
        case 4
            error('not ready yet');
        case 5
            %% gain pair first trials
            trialType_nm = '_gainPair_first';
            % onset stim
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_chosen',trialType_nm], bf_idx);
            
            % modulators
            % trial number
            if ismember(mod_chosen.trialN,[1,2,3])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx);
            end
            
            % SV = pA*QA+pB*QB
            if mod_chosen.SV ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_SV',trialType_nm], bf_idx);
            end
            
            % dQ
            if mod_chosen.dQ ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_dQ',trialType_nm], bf_idx);
            end
            
            % p(choice=best option)
            if mod_chosen.pBest ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_pBest',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_chosen.trialN,[4,5,6])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx);
            end
            
            % ROI activity
            if mod_chosen.ROI_activity_yn ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_GLM',mod_chosen.ROI_activity_GLM,...
                    '_',mod_chosen.ROI_activity_period,'_period_',mod_chosen.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx);
            end
            
            % RT
            if mod_chosen.RT ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_RT',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_chosen.trialN,[7,8,9])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx);
            end
            
            %% neutral pair first trials
            trialType_nm = '_ntalPair_first';
            % onset stim
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_chosen',trialType_nm], bf_idx);
            
            % modulators
            % trial number
            if ismember(mod_chosen.trialN,[1,2,3])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx);
            end
            
            % no SV, no dQ for neutral pair
            
            % trial number
            if ismember(mod_chosen.trialN,[4,5,6])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx);
            end
            
            % ROI activity
            if mod_chosen.ROI_activity_yn ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_GLM',mod_chosen.ROI_activity_GLM,...
                    '_',mod_chosen.ROI_activity_period,'_period_',mod_chosen.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx);
            end
            
            % RT
            if mod_chosen.RT ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_RT',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_chosen.trialN,[7,8,9])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx);
            end
            
            
            %% loss pair first trials
            trialType_nm = '_lossPair_first';
            % onset stim
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_chosen',trialType_nm], bf_idx);
            
            % modulators
            % trial number
            if ismember(mod_chosen.trialN,[1,2,3])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx);
            end
            
            % SV = pA*QA+pB*QB
            if mod_chosen.SV ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_SV',trialType_nm], bf_idx);
            end
            
            % dQ
            if mod_chosen.dQ ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_dQ',trialType_nm], bf_idx);
            end
            
            % p(choice = best option)
            if mod_chosen.pBest ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_pBest',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_chosen.trialN,[4,5,6])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx);
            end
            
            % ROI activity
            if mod_chosen.ROI_activity_yn ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_GLM',mod_chosen.ROI_activity_GLM,...
                    '_',mod_chosen.ROI_activity_period,'_period_',mod_chosen.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx);
            end
            
            % RT
            if mod_chosen.RT ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_RT',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_chosen.trialN,[7,8,9])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx);
            end
            
            %% gain pair last trials
            trialType_nm = '_gainPair_last';
            % onset stim
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_chosen',trialType_nm], bf_idx);
            
            % modulators
            % trial number
            if ismember(mod_chosen.trialN,[1,2,3])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx);
            end
            
            % SV = pA*QA+pB*QB
            if mod_chosen.SV ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_SV',trialType_nm], bf_idx);
            end
            
            % dQ
            if mod_chosen.dQ ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_dQ',trialType_nm], bf_idx);
            end
            
            % p(choice=best option)
            if mod_chosen.pBest ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_pBest',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_chosen.trialN,[4,5,6])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx);
            end
            
            % ROI activity
            if mod_chosen.ROI_activity_yn ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_GLM',mod_chosen.ROI_activity_GLM,...
                    '_',mod_chosen.ROI_activity_period,'_period_',mod_chosen.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx);
            end
            
            % RT
            if mod_chosen.RT ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_RT',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_chosen.trialN,[7,8,9])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx);
            end
            
            %% neutral pair last trials
            trialType_nm = '_ntalPair_last';
            % onset stim
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_chosen',trialType_nm], bf_idx);
            
            % modulators
            % trial number
            if ismember(mod_chosen.trialN,[1,2,3])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx);
            end
            
            % no SV, no dQ for neutral pair
            
            % trial number
            if ismember(mod_chosen.trialN,[4,5,6])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx);
            end
            
            % ROI activity
            if mod_chosen.ROI_activity_yn ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_GLM',mod_chosen.ROI_activity_GLM,...
                    '_',mod_chosen.ROI_activity_period,'_period_',mod_chosen.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx);
            end
            
            % RT
            if mod_chosen.RT ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_RT',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_chosen.trialN,[7,8,9])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx);
            end
            
            %% loss pair last trials
            trialType_nm = '_lossPair_last';
            % onset stim
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_chosen',trialType_nm], bf_idx);
            
            % modulators
            % trial number
            if ismember(mod_chosen.trialN,[1,2,3])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx);
            end
            
            % SV = pA*QA+pB*QB
            if mod_chosen.SV ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_SV',trialType_nm], bf_idx);
            end
            
            % dQ
            if mod_chosen.dQ ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_dQ',trialType_nm], bf_idx);
            end
            
            % p(choice = best option)
            if mod_chosen.pBest ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_pBest',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_chosen.trialN,[4,5,6])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx);
            end
            
            % ROI activity
            if mod_chosen.ROI_activity_yn ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_GLM',mod_chosen.ROI_activity_GLM,...
                    '_',mod_chosen.ROI_activity_period,'_period_',mod_chosen.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx);
            end
            
            % RT
            if mod_chosen.RT ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_RT',trialType_nm], bf_idx);
            end
            
            % trial number
            if ismember(mod_chosen.trialN,[7,8,9])
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_chosen_trialN',trialType_nm], bf_idx);
            end
    end % chosen option in red
    
    %% feedback
    o_fbk = RLprm.o_fbk;
    mod_fbk = RLprm.mod_fbk;
    switch o_fbk
        case 1
            % onset feedback
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_fbk'], bf_idx);
            
            % modulators
            % trial number
            if mod_fbk.trialN ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN'], bf_idx);
            end
            
            % feedback
            if mod_fbk.fbk ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_fbk'], bf_idx);
            end
            
            % PE
            if mod_fbk.PE ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE'], bf_idx);
            end
            
            % PE
            if mod_fbk.PE_bis ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE_bis'], bf_idx);
            end
            
            % total gain
            if mod_fbk.totalGain ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_totalGain'], bf_idx);
            end
            
            % ROI activity
            if mod_fbk.ROI_activity_yn ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GLM',mod_fbk.ROI_activity_GLM,...
                    '_',mod_fbk.ROI_activity_period,'_period_',mod_fbk.ROI_activity_ROI_nm,'_activity'], bf_idx);
            end
            
        case 2 % separate depending on pair type (gain/neutral/loss)
            
            %% gain pair
            trialType_nm = '_gainPair';
            % onset feedback
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_fbk',trialType_nm], bf_idx);
            
            % modulators
            % trial number
            if mod_fbk.trialN ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN',trialType_nm], bf_idx);
            end
            
            % feedback
            if mod_fbk.fbk ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_fbk',trialType_nm], bf_idx);
            end
            
            % PE
            if mod_fbk.PE ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE',trialType_nm], bf_idx);
            end
            
            % PE bis
            if mod_fbk.PE_bis ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE_bis',trialType_nm], bf_idx);
            end
            
            % total gain
            if mod_fbk.totalGain ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_totalGain',trialType_nm], bf_idx);
            end
            
            % ROI activity
            if mod_fbk.ROI_activity_yn ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GLM',mod_fbk.ROI_activity_GLM,...
                    '_',mod_fbk.ROI_activity_period,'_period_',mod_fbk.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx);
            end
            
            %% neutral pair
            trialType_nm = '_ntalPair';
            % onset feedback
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_fbk',trialType_nm], bf_idx);
            
            % modulators
            % trial number
            if mod_fbk.trialN ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN',trialType_nm], bf_idx);
            end
            
            % total gain
            if mod_fbk.totalGain ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_totalGain',trialType_nm], bf_idx);
            end
            
            % ROI activity
            if mod_fbk.ROI_activity_yn ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GLM',mod_fbk.ROI_activity_GLM,...
                    '_',mod_fbk.ROI_activity_period,'_period_',mod_fbk.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx);
            end
            
            %% loss pair
            trialType_nm = '_lossPair';
            % onset feedback
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_fbk',trialType_nm], bf_idx);
            
            % modulators
            % trial number
            if mod_fbk.trialN ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN',trialType_nm], bf_idx);
            end
            
            % feedback
            if mod_fbk.fbk ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_fbk',trialType_nm], bf_idx);
            end
            
            % PE
            if mod_fbk.PE ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE',trialType_nm], bf_idx);
            end
            
            % PE bis
            if mod_fbk.PE_bis ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE_bis',trialType_nm], bf_idx);
            end
            
            % total gain
            if mod_fbk.totalGain ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_totalGain',trialType_nm], bf_idx);
            end
            
            % ROI activity
            if mod_fbk.ROI_activity_yn ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GLM',mod_fbk.ROI_activity_GLM,...
                    '_',mod_fbk.ROI_activity_period,'_period_',mod_fbk.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx);
            end
            
        case 3 % separate depending on feedback type (gain/neutral/loss)
            %% gain feedback
            trialType_nm = '_gainFbk';
            % onset feedback
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_fbk',trialType_nm], bf_idx);
            
            % modulators
            % trial number
            if mod_fbk.trialN ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN',trialType_nm], bf_idx);
            end
            
            % feedback
            if mod_fbk.fbk ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_fbk',trialType_nm], bf_idx);
            end
            
            % PE
            if mod_fbk.PE ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE',trialType_nm], bf_idx);
            end
            
            % PE bis
            if mod_fbk.PE_bis ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE_bis',trialType_nm], bf_idx);
            end
            
            % total gain
            if mod_fbk.totalGain ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_totalGain',trialType_nm], bf_idx);
            end
            
            % ROI activity
            if mod_fbk.ROI_activity_yn ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GLM',mod_fbk.ROI_activity_GLM,...
                    '_',mod_fbk.ROI_activity_period,'_period_',mod_fbk.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx);
            end
            
            %% neutral feedback
            trialType_nm = '_ntalFbk';
            % onset feedback
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_fbk',trialType_nm], bf_idx);
            
            % modulators
            % trial number
            if mod_fbk.trialN ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN',trialType_nm], bf_idx);
            end
            
            % feedback: none always 0 for neutral feedback...
            
            % PE
            if mod_fbk.PE ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE',trialType_nm], bf_idx);
            end
            
            % PE bis
            if mod_fbk.PE_bis ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE_bis',trialType_nm], bf_idx);
            end
            
            % total gain
            if mod_fbk.totalGain ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_totalGain',trialType_nm], bf_idx);
            end
            
            % ROI activity
            if mod_fbk.ROI_activity_yn ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GLM',mod_fbk.ROI_activity_GLM,...
                    '_',mod_fbk.ROI_activity_period,'_period_',mod_fbk.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx);
            end
            
            %% loss feedback
            trialType_nm = '_lossFbk';
            % onset feedback
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_fbk',trialType_nm], bf_idx);
            
            % modulators
            % trial number
            if mod_fbk.trialN ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN',trialType_nm], bf_idx);
            end
            
            % feedback
            if mod_fbk.fbk ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_fbk',trialType_nm], bf_idx);
            end
            
            % PE
            if mod_fbk.PE ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE',trialType_nm], bf_idx);
            end
            
            % PE bis
            if mod_fbk.PE_bis ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE_bis',trialType_nm], bf_idx);
            end
            
            % total gain
            if mod_fbk.totalGain ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_totalGain',trialType_nm], bf_idx);
            end
            
            % ROI activity
            if mod_fbk.ROI_activity_yn ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GLM',mod_fbk.ROI_activity_GLM,...
                    '_',mod_fbk.ROI_activity_period,'_period_',mod_fbk.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx);
            end
            
        case 4 % separate depending on feedback AND pair type
            %% gain pair - gain feedback
            trialType_nm = '_gainPair_gainFbk';
            % onset feedback
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_fbk',trialType_nm], bf_idx);
            
            % modulators
            % trial number
            if mod_fbk.trialN ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN',trialType_nm], bf_idx);
            end
            
            % feedback
            if mod_fbk.fbk ~= 0
                error('Pointless: feedback = stable for a given feedback type.');
            end
            
            % PE
            if mod_fbk.PE ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE',trialType_nm], bf_idx);
            end
            
            % PE bis
            if mod_fbk.PE_bis ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE_bis',trialType_nm], bf_idx);
            end
            
            % total gain
            if mod_fbk.totalGain ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_totalGain',trialType_nm], bf_idx);
            end
            
            % ROI activity
            if mod_fbk.ROI_activity_yn ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GLM',mod_fbk.ROI_activity_GLM,...
                    '_',mod_fbk.ROI_activity_period,'_period_',mod_fbk.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx);
            end
            
            %% gain pair - neutral feedback
            trialType_nm = '_gainPair_ntalFbk';
            % onset feedback
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_fbk',trialType_nm], bf_idx);
            
            % modulators
            % trial number
            if mod_fbk.trialN ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN',trialType_nm], bf_idx);
            end
            
            % feedback
            if mod_fbk.fbk ~= 0
                error('Pointless: feedback = stable for a given feedback type.');
            end
            
            % PE
            if mod_fbk.PE ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE',trialType_nm], bf_idx);
            end
            
            % PE bis
            if mod_fbk.PE_bis ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE_bis',trialType_nm], bf_idx);
            end
            
            % total gain
            if mod_fbk.totalGain ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_totalGain',trialType_nm], bf_idx);
            end
            
            % ROI activity
            if mod_fbk.ROI_activity_yn ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GLM',mod_fbk.ROI_activity_GLM,...
                    '_',mod_fbk.ROI_activity_period,'_period_',mod_fbk.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx);
            end
            
            %% neutral pair
            trialType_nm = '_ntalPair_ntalFbk';
            % onset feedback
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_fbk',trialType_nm], bf_idx);
            
            % modulators
            % trial number
            if mod_fbk.trialN ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN',trialType_nm], bf_idx);
            end
            
            % feedback
            if mod_fbk.fbk ~= 0
                error('Pointless: feedback = stable for a given feedback type.');
            end
            
            % PE: pointless should remain equal 0 always for neutral
            % pair
            
            % total gain
            if mod_fbk.totalGain ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_totalGain',trialType_nm], bf_idx);
            end
            
            % ROI activity
            if mod_fbk.ROI_activity_yn ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GLM',mod_fbk.ROI_activity_GLM,...
                    '_',mod_fbk.ROI_activity_period,'_period_',mod_fbk.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx);
            end
            
            %% loss pair
            trialType_nm = '_lossPair_ntalFbk';
            % onset feedback
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_fbk',trialType_nm], bf_idx);
            
            % modulators
            % trial number
            if mod_fbk.trialN ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN',trialType_nm], bf_idx);
            end
            
            % feedback
            if mod_fbk.fbk ~= 0
                error('Pointless: feedback = stable for a given feedback type.');
            end
            
            % PE
            if mod_fbk.PE ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE',trialType_nm], bf_idx);
            end
            
            % PE bis
            if mod_fbk.PE_bis ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE_bis',trialType_nm], bf_idx);
            end
            
            % total gain
            if mod_fbk.totalGain ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_totalGain',trialType_nm], bf_idx);
            end
            
            % ROI activity
            if mod_fbk.ROI_activity_yn ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GLM',mod_fbk.ROI_activity_GLM,...
                    '_',mod_fbk.ROI_activity_period,'_period_',mod_fbk.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx);
            end
            
            %% loss pair - loss feedback
            trialType_nm = '_lossPair_lossFbk';
            % onset feedback
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_fbk',trialType_nm], bf_idx);
            
            % modulators
            % trial number
            if mod_fbk.trialN ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN',trialType_nm], bf_idx);
            end
            
            % feedback
            if mod_fbk.fbk ~= 0
                error('Pointless: feedback = stable for a given feedback type.');
            end
            
            % PE
            if mod_fbk.PE ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE',trialType_nm], bf_idx);
            end
            
            % PE bis
            if mod_fbk.PE_bis ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE_bis',trialType_nm], bf_idx);
            end
            
            % total gain
            if mod_fbk.totalGain ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_totalGain',trialType_nm], bf_idx);
            end
            
            % ROI activity
            if mod_fbk.ROI_activity_yn ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GLM',mod_fbk.ROI_activity_GLM,...
                    '_',mod_fbk.ROI_activity_period,'_period_',mod_fbk.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx);
            end
            
        case 5 % pool gain+loss pairs, separate neutral apart
            
            %% gain +loss pair
            trialType_nm = '_GL_Pairs';
            % onset feedback
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_fbk',trialType_nm], bf_idx);
            
            % modulators
            % trial number
            if mod_fbk.trialN ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN',trialType_nm], bf_idx);
            end
            
            % feedback
            if mod_fbk.fbk ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_fbk',trialType_nm], bf_idx);
            end
            
            % PE
            if mod_fbk.PE ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE',trialType_nm], bf_idx);
            end
            
            % PE bis
            if mod_fbk.PE_bis ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE_bis',trialType_nm], bf_idx);
            end
            
            % total gain
            if mod_fbk.totalGain ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_totalGain',trialType_nm], bf_idx);
            end
            
            % ROI activity
            if mod_fbk.ROI_activity_yn ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GLM',mod_fbk.ROI_activity_GLM,...
                    '_',mod_fbk.ROI_activity_period,'_period_',mod_fbk.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx);
            end
            
            %% neutral pair
            trialType_nm = '_ntalPair';
            % onset feedback
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_fbk',trialType_nm], bf_idx);
            
            % modulators
            % trial number
            if mod_fbk.trialN ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN',trialType_nm], bf_idx);
            end
            
            % total gain
            if mod_fbk.totalGain ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_totalGain',trialType_nm], bf_idx);
            end
            
            % ROI activity
            if mod_fbk.ROI_activity_yn ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GLM',mod_fbk.ROI_activity_GLM,...
                    '_',mod_fbk.ROI_activity_period,'_period_',mod_fbk.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx);
            end
            
        case 6 % pool gain+loss pairs (but split first/second half trials) and keep neutral pair apart
            
            %% gain +loss pair (first trials)
            trialType_nm = '_GL_Pairs_first';
            % onset feedback
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_fbk',trialType_nm], bf_idx);
            
            % modulators
            % trial number
            if mod_fbk.trialN ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN',trialType_nm], bf_idx);
            end
            
            % feedback
            if mod_fbk.fbk ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_fbk',trialType_nm], bf_idx);
            end
            
            % PE
            if mod_fbk.PE ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE',trialType_nm], bf_idx);
            end
            
            % PE bis
            if mod_fbk.PE_bis ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE_bis',trialType_nm], bf_idx);
            end
            
            % total gain
            if mod_fbk.totalGain ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_totalGain',trialType_nm], bf_idx);
            end
            
            % ROI activity
            if mod_fbk.ROI_activity_yn ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GLM',mod_fbk.ROI_activity_GLM,...
                    '_',mod_fbk.ROI_activity_period,'_period_',mod_fbk.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx);
            end
            
            %% gain +loss pair (last trials)
            trialType_nm = '_GL_Pairs_last';
            % onset feedback
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_fbk',trialType_nm], bf_idx);
            
            % modulators
            % trial number
            if mod_fbk.trialN ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN',trialType_nm], bf_idx);
            end
            
            % feedback
            if mod_fbk.fbk ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_fbk',trialType_nm], bf_idx);
            end
            
            % PE
            if mod_fbk.PE ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE',trialType_nm], bf_idx);
            end
            
            % PE bis
            if mod_fbk.PE_bis ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_PE_bis',trialType_nm], bf_idx);
            end
            
            % total gain
            if mod_fbk.totalGain ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_totalGain',trialType_nm], bf_idx);
            end
            
            % ROI activity
            if mod_fbk.ROI_activity_yn ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GLM',mod_fbk.ROI_activity_GLM,...
                    '_',mod_fbk.ROI_activity_period,'_period_',mod_fbk.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx);
            end
            
            %% neutral pair
            trialType_nm = '_ntalPair';
            % onset feedback
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_fbk',trialType_nm], bf_idx);
            
            % modulators
            % trial number
            if mod_fbk.trialN ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_trialN',trialType_nm], bf_idx);
            end
            
            % total gain
            if mod_fbk.totalGain ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_totalGain',trialType_nm], bf_idx);
            end
            
            % ROI activity
            if mod_fbk.ROI_activity_yn ~= 0
                [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'mod_fbk_GLM',mod_fbk.ROI_activity_GLM,...
                    '_',mod_fbk.ROI_activity_period,'_period_',mod_fbk.ROI_activity_ROI_nm,'_activity',trialType_nm], bf_idx);
            end
    end
    
    %% cross
    o_cross = RLprm.o_cross;
    switch o_cross
        case 1
            % onset cross
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_cross'], bf_idx);
    end
    
    %% missed trials
    o_missed_trials_stim = RLprm.o_missed_trials_stim;
    switch o_missed_trials_stim
        case 1
            % onset answer
            [n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, [task_id,'o_missed_trials_stim'], bf_idx);
    end
    
    %% add movement parameters at the end of each run
    n_prm.all = n_prm.all + n_mvmt_perRun(iRun);
    n_prm.(task_nm).(run_nm) = n_prm.(task_nm).(run_nm) + n_mvmt_perRun(iRun);
    bf_idx = bf_idx + n_mvmt_perRun(iRun);
end % run loop

%% add one constant per task at the end
n_prm.all = n_prm.all + n_total_runs;

end % function end


function[n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, curr_reg_nm, bf_idx)
%[n_prm, prm_idx, bf_idx] = update_n_prm(n_prm, prm_idx, task_nm, run_nm, nBperR, curr_reg_nm, bf_idx)
% sub-function to update each regressor with one single line

%% update total number of parameters
n_prm.all = n_prm.all + nBperR; 
% update RL number of parameters for the current run
n_prm.(task_nm).(run_nm) = n_prm.(task_nm).(run_nm) + nBperR; 

%% extract index of the regressor
curr_reg_idx = bf_idx + 1;
% extract index for contrast across ALL runs
if ~isfield(prm_idx, curr_reg_nm)
    prm_idx.(curr_reg_nm) = curr_reg_idx; % extract index for this regressor
else % means that the current regressor has already been extracted in a previous run for example
    % => extract all the index when it is used inside the same variable
    prm_idx.(curr_reg_nm) = [prm_idx.(curr_reg_nm), curr_reg_idx];
end
% extract index for contrast specific to each run
prm_idx.([curr_reg_nm,'_',run_nm]) = curr_reg_idx; % extract index for this regressor

%% update index (need to be done after extracting current regressor index)
bf_idx = bf_idx + nBperR;

end