
%% install the VBA toolbox if not present
wherever = pwd;
if ~exist('RFT_GLM_contrast.m','file')
    VBA_folder = fullfile('C:', 'Users', 'nicolas.clairis', 'Documents', 'GitHub', 'VBA-toolbox');
    cd(VBA_folder);
    VBA_setup;
    cd(wherever);
end

%% define model parameters
mdlPrm.fatigueEvolOrObs = questdlg('Fatigue as hidden state or observation parameter?',...
    'Fatigue in f or g function',...
    'hidden state','observation parameter',...
    'observation parameter'); % the fatigue component can either be considered as a hidden state evolving 
% trial by trial, ie included in the evolution function; or it can be
% directly included in the observation function as kR and kE

mdlPrm.fatigueType = questdlg('How do you want to model fatigue?',...
    'Fatigue variable',...
    'trial number','sum of efforts','effort integral',...
    'trial number'); % consider trial number, sum of effort levels until current trial or integral of effort performed

mdlPrm.kR_tasksSplitOrPool = questdlg('Pool tasks together to estimate kR or check them separately?',...
    '1 kR per subject or per task per subject',...
    'pool','split',...
    'split');

%% prepare the model variables to use for VBA model estimation
nTrialsPerSession = 44;
nSessionsPerTask = 2; % number of runs per task
nTrialsPerTask = nTrialsPerSession*nSessionsPerTask;
task_names = {'mental','physical'};
nTasks = length(task_names); % physical/mental effort
nTotalTrials = nTrialsPerTask*nTasks;
nTotalSessions = nSessionsPerTask*nTasks;

% initialize parameters of interest
switch mdlPrm.kR_tasksSplitOrPool
    case 'pool'
        % one variable to explain for both tasks
        [choice_left,...
            trialN,...
            Rleft, Rright,...
            EleftPool_Em, ErightPool_Em,...
            EleftPool_Ep, ErightPool_Ep,...
            sumEprevPool_Em, intEprevPool_Em,...
            sumEprevPool_Ep, intEprevPool_Ep,...
            effort_type,...
            ElevelPrevTrialPool_Em, ElevelPrevTrialPool_Ep,...
            EintegralPrevTrialPool_Em, EintegralPrevTrialPool_Ep] = deal(NaN(1, nTotalTrials));
    case 'split'
        % model each task separately
        [choice_left_Ep,...
            trialN_Ep,...
            Rleft_Ep, Rright_Ep,...
            Eleft_Ep, Eright_Ep,...
            sumEprev_Ep, intEprev_Ep,...
            ElevelPrevTrial_Ep,...
            EintegralPrevTrial_Ep] = deal( NaN(1, nTrialsPerTask) );
        [choice_left_Em,...
            trialN_Em,...
            Rleft_Em, Rright_Em,...
            Eleft_Em, Eright_Em,...
            sumEprev_Em, intEprev_Em,...
            ElevelPrevTrial_Em,...
            EintegralPrevTrial_Em] = deal( NaN(1, nTrialsPerTask) );
end

% load the data
for iTask = 1:nTasks
    task_nm = task_names{iTask};
    
    for iSession = 1:nSessionPerTask
        
        
        %% load data
        loadStruct = load(['AB_s000_session0',num2str(iSession),'_',task_nm,'_task_messyAllStuff.mat']);
        
        % extract relevant parameters
        % variable to explain: y
        choice_tmp = loadStruct.choice;% choice=1*nTrials vector with the choices (-1 left, +1 right)
        % convert into binary variable
        choice_tmp(loadStruct.choice == -1) = true; % left choice
        choice_tmp(loadStruct.choice == 1) = false; % right choice
        
        % other variables of interest
        Rleft_tmp = loadStruct.choice_opt.R.left; % 1*nTrials vector with the reward level for left option
        Rright_tmp = loadStruct.choice_opt.R.right; % 1*nTrials vector with the reward level for right option
        Eleft_tmp = loadStruct.choice_opt.E.left; % 1*nTrials vector with the effort level for left option
        Eright_tmp = loadStruct.choice_opt.E.right; % 1*nTrials vector with the effort level for right option
        trialN_tmp = 1:nTrialsPerSession;
        Echosen_tmp = loadStruct.E_chosen; % 1*nTrials vector with the effort level chosen for each trial
        trialSuccess_tmp = loadStruct.trial_was_successfull; % 1*nTrials binary vector telling you if they reached
        perfSummary_tmp = loadStruct.perfSummary; % structure with details of performance for each trial
        [sumEprev_tmp, intEprev_tmp,...
            ElevelPrevTrial_tmp, EintegralPrevTrial_tmp] = deal(zeros(1,nTrialsPerSession));
        for iTrial = trialN_tmp(2:end) % first trial sum of effort = 0 and same for integral
            sumEprev_tmp(iTrial) = sumEprev_tmp(iTrial - 1) + Echosen_tmp(iTrial-1).*trialSuccess_tmp(iTrial-1);
            ElevelPrevTrial_tmp(iTrial) = Echosen_tmp(iTrial-1).*trialSuccess_tmp(iTrial-1);
            switch task_nm
                case 'mental'
                    intEprev_tmp(iTrial) = intEprev_tmp(iTrial - 1) + perfSummary_tmp{1,iTrial-1}.n_correctAnswersProvided;
                    EintegralPrevTrial_tmp(iTrial) = perfSummary_tmp{1,iTrial-1}.n_correctAnswersProvided;
                case 'physical'
                    error('define how to compute integral of effort for physical effort task');
            end
        end
        
        % store into variables for the model
        switch mdlPrm.kR_tasksSplitOrPool
            case 'pool'
                trials_idx = (1:nTrialsPerSession) + nTrialsPerSession*(iSession - 1) + nTrialsPerTask*(iTask - 1);
                choice_left(trials_idx) = choice_tmp;
                trialN(trials_idx) = trialN_tmp;
                Rleft(trials_idx) = Rleft_tmp;
                Rright(trials_idx) = Rright_tmp;
                switch task_nm
                    case 'mental'
                        EleftPool_Em(trials_idx)    = Eleft_tmp;
                        ErightPool_Em(trials_idx)   = Eright_tmp;
                        sumEprevPool_Em(trials_idx) = sumEprev_tmp;
                        intEprevPool_Em(trials_idx) = intEprev_tmp;
                        [EleftPool_Ep, ErightPool_Ep,...
                            sumEprevPool_Ep, intEprevPool_Ep,...
                            ElevelPrevTrialPool_Ep, EintegralPrevTrialPool_Ep] = deal(0);
                        effort_type(trials_idx) = 0;
                        ElevelPrevTrialPool_Em = ElevelPrevTrial_tmp;
                        EintegralPrevTrialPool_Em = EintegralPrevTrial_tmp;
                        
                    case 'physical'
                        EleftPool_Ep(trials_idx)    = Eleft_tmp;
                        ErightPool_Ep(trials_idx)   = Eright_tmp;
                        sumEprevPool_Ep(trials_idx) = sumEprev_tmp;
                        intEprevPool_Ep(trials_idx) = intEprev_tmp;
                        [EleftPool_Em, ErightPool_Em,...
                            sumEprevPool_Em, intEprevPool_Em,...
                            ElevelPrevTrialPool_Em, EintegralPrevTrialPool_Em] = deal(0);
                        effort_type(trials_idx) = 1;
                        ElevelPrevTrialPool_Ep(trials_idx) = ElevelPrevTrial_tmp;
                        EintegralPrevTrialPool_Ep(trials_idx) = EintegralPrevTrial_tmp;
                end
            case 'split'
                trials_idx = (1:nTrialsPerSession) + nTrialsPerSession*(iSession - 1);
                switch task_nm
                    case 'physical'
                        choice_left_Ep(trials_idx) = choice_tmp;
                        trialN_Ep(trials_idx) = trialN_tmp;
                        Rleft_Ep(trials_idx) = Rleft_tmp;
                        Rright_Ep(trials_idx) = Rright_tmp;
                        Eleft_Ep(trials_idx) = Eleft_tmp;
                        Eright_Ep(trials_idx) = Eright_tmp;
                        sumEprev_Ep(trials_idx) = sumEprev_tmp;
                        intEprev_Ep(trials_idx) = intEprev_tmp;
                        ElevelPrevTrial_Ep(trials_idx) = ElevelPrevTrial_tmp;
                        EintegralPrevTrial_Ep(trials_idx) = EintegralPrevTrial_tmp;
                    case 'mental'
                        choice_left_Em(trials_idx) = choice_tmp;
                        trialN_Em(trials_idx) = trialN_tmp;
                        Rleft_Em(trials_idx) = Rleft_tmp;
                        Rright_Em(trials_idx) = Rright_tmp;
                        Eleft_Em(trials_idx) = Eleft_tmp;
                        Eright_Em(trials_idx) = Eright_tmp;
                        sumEprev_Em(trials_idx) = sumEprev_tmp;
                        intEprev_Em(trials_idx) = intEprev_tmp;
                        ElevelPrevTrial_Em(trials_idx) = ElevelPrevTrial_tmp;
                        EintegralPrevTrial_Em(trials_idx) = EintegralPrevTrial_tmp;
                end
        end
        
    end % run loop
end % loop through physical and mental effort tasks

%% Define number of parameters to estimate
switch mdlPrm.fatigueEvolOrObs
    case 'hidden state'
        switch mdlPrm.kR_tasksSplitOrPool
            case 'pool'
                n_G_prm = 3; % kR, kEp, kEm
                % 1 fatigue component per task kF_Ep and kF_Em
                n_hiddenStates = 2; % kF physical and kF mental
                n_F_prm = 2; % fatigue component evolving through time
            case 'split'
                n_hiddenStates = 1;
                n_F_prm = 1; % fatigue component evolving through time
                n_G_prm = 2; % kR and kE for each task
        end
    case 'observation parameter' % no evolution function nor any hidden state
        n_hiddenStates = 0;
        n_F_prm = 0;
        switch mdlPrm.kR_tasksSplitOrPool
            case 'pool'
                n_G_prm = 5; % kR, kEp, kFEp, kEm, kFEm
            case 'split'
                n_G_prm = 3; % kR, kE and kF for each task
        end
end

switch mdlPrm.kR_tasksSplitOrPool
    case 'pool'
        nTrialsEstimation = nTotalTrials;
    case 'split'
        nTrialsEstimation = nTrialsPerTask;
end

%% Define VBA options

% enter the variables in VBA format
switch mdlPrm.kR_tasksSplitOrPool
    case 'pool'
        y = choice_left;
        u_t = [trialN; Rleft; Rright;...
            EleftPool_Em; ErightPool_Em;...
            EleftPool_Ep; ErightPool_Ep;...
            sumEprevPool_Em; intEprevPool_Em;...
            sumEprevPool_Ep; intEprevPool_Ep;...
            effort_type;...
            ElevelPrevTrialPool_Em; ElevelPrevTrialPool_Ep;...
            EintegralPrevTrialPool_Em; EintegralPrevTrialPool_Ep];
    case 'split'
        y_Ep = choice_left_Ep;
        y_Em = choice_left_Em;
        u_t_Ep = [trialN_Ep;...
            Rleft_Ep; Rright_Ep;...
            Eleft_Ep; Eright_Ep;...
            sumEprev_Ep; intEprev_Ep;...
            ElevelPrevTrial_Ep,...
            EintegralPrevTrial_Ep];
        u_t_Em = [trialN_Em;...
            Rleft_Em; Rright_Em;...
            Eleft_Em; Eright_Em;...
            sumEprev_Em; intEprev_Em;...
            ElevelPrevTrial_Em,...
            EintegralPrevTrial_Em];
end
% define number of parameters to estimate
dim = struct('n', n_hiddenStates,... % number of hidden states
    'n_t', nTrialsEstimation,... % number of trials across runs
    'n_theta',n_F_prm,... % number of evolution parameters
    'n_phi',n_G_prm); % number of observation parameters

options = struct;
options.sources.type = 1; % binary data
options.dim         = dim; % define number of parameters that need to be estimated by the model
options.GnFigs      = 0;
options.DisplayWin  = 1; % display figure during inversion (1) or not (0)
options.verbose     = 1; % display text during inversion (1) or not (0)

% store model parameters in the inputs of each function
if strcmp(mdlPrm.fatigueEvolOrObs,'hidden state')
    options.inF.mdlPrm = mdlPrm;
end
options.inG.mdlPrm = mdlPrm;

% define priors to use for your parameters
if strcmp(mdlPrm.fatigueEvolOrObs,'hidden state')
    % priors for fatigue parameter in evolution function
    options.priors.muTheta = zeros(n_F_prm,1);
    options.priors.SigmaTheta = eye(n_F_prm);
    % priors for initial fatigue state
    options.priors.muX0 = zeros(n_hiddenStates,1);
    options.priors.SigmaX0 = eye(n_hiddenStates);
end
% g parameters priors
options.priors.muPhi = zeros(n_G_prm, 1);
options.priors.SigmaPhi = eye(n_G_prm);


% use multisession to pool sessions together but considering that they can have
% different noise levels
switch mdlPrm.kR_tasksSplitOrPool
    case 'pool'
        options.multisession.split = repmat(nTrialsPerSession, 1, nTotalSessions);
    case 'split'
        options.multisession.split = repmat(nTrialsPerSession, 1, nSessionsPerTask);
end
% fix the parameters which should be equivalent across sessions and those
% which can vary independently in each session
switch mdlPrm.kR_tasksSplitOrPool
    case 'pool' % 
        
    case 'split' % physical and mental effort task estimated separately => fix the parameters within each task across sessions
        if strcmp(mdlPrm.fatigueEvolOrObs,'hidden state')
            options.multisession.fixed.theta = 1:n_F_prm;
        end
        options.multisession.fixed.phi = 1:n_G_prm;
end

%% Define functions
switch mdlPrm.fatigueEvolOrObs
    case 'hidden state'
        f_fname = @LGCM_f_evolution;
    case 'observation parameter'
        f_fname = []; % no evolution function if no hidden state considered
end
g_fname = @LGCM_g_observation;

%% Do Laplacian inversion of the model with the VBA
switch mdlPrm.kR_tasksSplitOrPool
    case 'pool'
        [posterior,out] = VBA_NLStateSpaceModel(choice_left, u_t, f_fname, g_fname, dim, options);
    case 'split'
        [posteriorEp,outEp] = VBA_NLStateSpaceModel(yEp, u_t_Ep, f_fname, g_fname, dim, options);
        [posteriorEm,outEm] = VBA_NLStateSpaceModel(yEm, u_t_Em, f_fname, g_fname, dim, options);
end