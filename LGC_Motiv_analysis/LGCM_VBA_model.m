
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

mdlPrm.kRP_tasksSplitOrPool = questdlg('Pool tasks together to estimate kR and kP or check them separately?',...
    '1 kR and 1 kP per subject or per task per subject',...
    'pool','split',...
    'split');

% extract model parameters
fatigueEvolOrObs    = mdlPrm.fatigueEvolOrObs;
fatigueType         = mdlPrm.fatigueType;
kRP_tasksSplitOrPool = mdlPrm.kRP_tasksSplitOrPool;

%% prepare the model variables to use for VBA model estimation
nTrialsPerSession = 44;
nSessionsPerTask = 2; % number of runs per task
nTrialsPerTask = nTrialsPerSession*nSessionsPerTask;
task_names = {'mental','physical'};
nTasks = length(task_names); % physical/mental effort
nTotalTrials = nTrialsPerTask*nTasks;
nTotalSessions = nSessionsPerTask*nTasks;

% initialize parameters of interest
switch kRP_tasksSplitOrPool
    case 'pool'
        % one variable to explain for both tasks
        [choice_left,...
            trialN,...
            Rleft, Rright,...
            R_or_P,...
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
            R_or_P_Ep,...
            Eleft_Ep, Eright_Ep,...
            sumEprev_Ep, intEprev_Ep,...
            ElevelPrevTrial_Ep,...
            EintegralPrevTrial_Ep] = deal( NaN(1, nTrialsPerTask) );
        [choice_left_Em,...
            trialN_Em,...
            Rleft_Em, Rright_Em,...
            R_or_P_Em,...
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
        choice_tmp(loadStruct.choice == 0) = NaN; % NO choice
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
        R_or_P_str_tmp = loadStruct.R_or_P; % load reward or punishment string
        % convert in a variable = to (0) for punishments and to (+1) for
        % rewards
        R_or_P_tmp = strcmp(R_or_P_str_tmp,'R');
        
        % store into variables for the model
        switch kRP_tasksSplitOrPool
            case 'pool'
                trials_idx = (1:nTrialsPerSession) + nTrialsPerSession*(iSession - 1) + nTrialsPerTask*(iTask - 1);
                choice_left(trials_idx) = choice_tmp;
                trialN(trials_idx) = trialN_tmp;
                Rleft(trials_idx) = Rleft_tmp;
                Rright(trials_idx) = Rright_tmp;
                R_or_P(trials_idx) = R_or_P_tmp;
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
                        R_or_P_Ep(trials_idx) = R_or_P_tmp;
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
                        R_or_P_Em(trials_idx) = R_or_P_tmp;
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
switch fatigueEvolOrObs
    case 'hidden state'
        switch kRP_tasksSplitOrPool
            case 'pool'
                G_Phi_prm_nm = {'kR','kP','kEm','kEp'};
                % 1 fatigue component per task kF_Ep and kF_Em
                F_Phi_prm_nm = {'kFEm','kFE'}; % Fatigue physical and Fatigue mental
            case 'split'
                G_Phi_prm_nm = {'kR','kP','kE'};
                % 1 fatigue component
                F_Phi_prm_nm = {'kF'}; % fatigue component evolving through time
        end
        n_G_prm = length(G_Phi_prm_nm);
        n_F_prm = length(F_Phi_prm_nm); 
        n_hiddenStates = n_F_prm; % (hidden states = fatigue components = number of F parameters)
    case 'observation parameter' % no evolution function nor any hidden state
        % no hidden states in this case (=no evolution function)
        F_Phi_prm_nm = {''};
        n_hiddenStates = 0;
        n_F_prm = 0;
        % observation parameters
        switch kRP_tasksSplitOrPool
            case 'pool'
                G_Phi_prm_nm = {'kR','kP','kEm','kEp','kFEm','kFEp'};
            case 'split'
                G_Phi_prm_nm = {'kR','kP','kE','kF'};
        end
        n_G_prm = length(G_Phi_prm_nm);
end

switch kRP_tasksSplitOrPool
    case 'pool'
        nTrialsEstimation = nTotalTrials;
    case 'split'
        nTrialsEstimation = nTrialsPerTask;
end

%% Define VBA options

% enter the variables in VBA format
switch kRP_tasksSplitOrPool
    case 'pool'
        u_t = [trialN;...
            Rleft; Rright;...
            EleftPool_Em; ErightPool_Em;...
            EleftPool_Ep; ErightPool_Ep;...
            sumEprevPool_Em; intEprevPool_Em;...
            sumEprevPool_Ep; intEprevPool_Ep;...
            effort_type;...
            ElevelPrevTrialPool_Em; ElevelPrevTrialPool_Ep;...
            EintegralPrevTrialPool_Em; EintegralPrevTrialPool_Ep;...
            R_or_P];
    case 'split'
        u_t_Ep = [trialN_Ep;...
            Rleft_Ep; Rright_Ep;...
            Eleft_Ep; Eright_Ep;...
            sumEprev_Ep; intEprev_Ep;...
            ElevelPrevTrial_Ep,...
            EintegralPrevTrial_Ep;...
            R_or_P_Ep];
        u_t_Em = [trialN_Em;...
            Rleft_Em; Rright_Em;...
            Eleft_Em; Eright_Em;...
            sumEprev_Em; intEprev_Em;...
            ElevelPrevTrial_Em,...
            EintegralPrevTrial_Em;...
            R_or_P_Em];
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
if strcmp(fatigueEvolOrObs,'hidden state')
    options.inF.mdlPrm = mdlPrm;
end
options.inG.mdlPrm = mdlPrm;

%% define priors to use for your parameters
if strcmp(fatigueEvolOrObs,'hidden state')
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


%% use multisession to pool sessions together but considering that they can have
% different noise levels
switch kRP_tasksSplitOrPool
    case 'pool'
        options.multisession.split = repmat(nTrialsPerSession, 1, nTotalSessions);
    case 'split'
        options.multisession.split = repmat(nTrialsPerSession, 1, nSessionsPerTask);
end
% fix the parameters which should be equivalent across sessions and those
% which can vary independently in each session
switch kRP_tasksSplitOrPool
    case 'pool' % need to see how to fix kR across sessions but not the parameters specific to each task
        if strcmp(fatigueEvolOrObs,'hidden state')
            options.multisession.fixed.theta = 1:n_F_prm;
        end
        options.multisession.fixed.phi = 1:n_G_prm;
    case 'split' % physical and mental effort task estimated separately => fix the parameters within each task across sessions
        if strcmp(fatigueEvolOrObs,'hidden state')
            options.multisession.fixed.theta = 1:n_F_prm;
        end
        options.multisession.fixed.phi = 1:n_G_prm;
end

%% need to define options separetely (for trials to ignore in particular)
% for each model when tasks are studied separately
switch kRP_tasksSplitOrPool
    case 'split'
        options_Em = options;
        options_Ep = options;
end

%% select trials to ignore

% ignore trials where no answer was provided for the estimation of the model
switch kRP_tasksSplitOrPool
    case 'split'
        % mental effort
        noChoiceTrials_Em = isnan(choice_left_Em);
        options_Em.isYout = noChoiceTrials_Em;
        % physical effort
        noChoiceTrials_Ep = isnan(choice_left_Ep);
        options_Ep.isYout = noChoiceTrials_Ep;
    case 'pool'
        noChoiceTrials = isnan(choice_left);
        options.isYout = noChoiceTrials;
end

% avoid updating hidden states (fatigue) for the first trial and for any
% trial where there was no effort performance (because choice too slow or
% performance failure)
switch kRP_tasksSplitOrPool
    case 'split'
        options_Em.skipf = zeros(nTrialsEstimation,1);
        options_Ep.skipf = zeros(nTrialsEstimation,1);
        % skip first trial
        options_Em.skipf(trialN_Em == 1) = 1;
        options_Ep.skipf(trialN_Ep == 1) = 1;
        % skip choice failure
        post_noChoiceTrials_Em = [false, noChoiceTrials_Em(1:(end-1))];
        options_Em.skipf(post_noChoiceTrials_Em) = 1;
        post_noChoiceTrials_Ep = [false, noChoiceTrials_Ep(1:(end-1))];
        options_Ep.skipf(post_noChoiceTrials_Ep) = 1;
    case 'pool'
        options.skipf = zeros(nTrialsEstimation,1);
        % skip first trial
        options.skipf(trialN == 1) = 1;
        % skip choice failure
        post_noChoiceTrials = [false, noChoiceTrials(1:(end-1))];
        options.skipf(post_noChoiceTrials) = 1;
end

%% Define functions
switch fatigueEvolOrObs
    case 'hidden state'
        f_fname = @LGCM_f_evolution;
    case 'observation parameter'
        f_fname = []; % no evolution function if no hidden state considered
end
g_fname = @LGCM_g_observation;

%% Do Laplacian inversion of the model with the VBA
switch kRP_tasksSplitOrPool
    case 'pool'
        [posterior, out] = VBA_NLStateSpaceModel(choice_left, u_t, f_fname, g_fname, dim, options);
    case 'split'
        [posterior_Ep, outEp] = VBA_NLStateSpaceModel(choice_left_Ep, u_t_Ep, f_fname, g_fname, dim, options_Em);
        [posterior_Em, outEm] = VBA_NLStateSpaceModel(choice_left_Em, u_t_Em, f_fname, g_fname, dim, options_Ep);
end

%% extract estimated parameters

% note: if you want to add some constraints on the parameters (ie > 0 or
% else) you might want to use functions like exp/log to transform them
switch kRP_tasksSplitOrPool
    case 'pool'
        switch fatigueEvolOrObs
            case 'hidden state'
                % extract mean
                kFEm.mean = posterior.muTheta(1);
                kFEp.mean = posterior.muTheta(2);
                kR.mean      = posterior.muPhi(1);
                kP.mean      = posterior.muPhi(2);
                kEm.mean      = posterior.muPhi(3);
                kEp.mean      = posterior.muPhi(4);
                % extract variance
                kFEm.sigma = posterior.SigmaTheta(1);
                kFEp.sigma = posterior.SigmaTheta(2);
                kR.sigma      = posterior.SigmaPhi(1);
                kP.sigma      = posterior.SigmaPhi(2);
                kEm.sigma      = posterior.SigmaPhi(3);
                kEp.sigma      = posterior.SigmaPhi(4);
            case 'observation parameter'
                % extract mean
                kR.mean      = posterior.muPhi(1);
                kP.mean      = posterior.muPhi(2);
                kEm.mean      = posterior.muPhi(3);
                kEp.mean      = posterior.muPhi(4);
                kFEm.mean      = posterior.muPhi(5);
                kFEp.mean      = posterior.muPhi(6);
                % extract variance
                kR.sigma      = posterior.SigmaPhi(1);
                kP.sigma      = posterior.SigmaPhi(2);
                kEm.sigma      = posterior.SigmaPhi(3);
                kEp.sigma      = posterior.SigmaPhi(4);
                kFEm.sigma      = posterior.SigmaPhi(5);
                kFEp.sigma      = posterior.SigmaPhi(6);
        end
    case 'split'
        switch fatigueEvolOrObs
            case 'hidden state'
                % extract mean
                kFEm.mean = posterior_Em.muTheta(1);
                kFEp.mean = posterior_Ep.muTheta(1);
                kR_Em.mean      = posterior_Em.muPhi(1);
                kR_Ep.mean      = posterior_Ep.muPhi(1);
                kP_Em.mean      = posterior_Em.muPhi(2);
                kP_Ep.mean      = posterior_Ep.muPhi(2);
                kEm.mean      = posterior_Em.muPhi(3);
                kEp.mean      = posterior_Ep.muPhi(3);
                % extract variance
                kFEm.sigma = posterior_Em.SigmaTheta(1);
                kFEp.sigma = posterior_Ep.SigmaTheta(1);
                kR_Em.sigma      = posterior_Em.SigmaPhi(1);
                kR_Ep.sigma      = posterior_Ep.SigmaPhi(1);
                kP_Em.sigma      = posterior_Em.SigmaPhi(2);
                kP_Ep.sigma      = posterior_Ep.SigmaPhi(2);
                kEm.sigma      = posterior_Em.SigmaPhi(3);
                kEp.sigma      = posterior_Ep.SigmaPhi(3);
            case 'observation parameter'
                % extract mean
                kR_Em.mean      = posterior_Em.muPhi(1);
                kR_Ep.mean      = posterior_Ep.muPhi(1);
                kP_Em.mean      = posterior_Em.muPhi(2);
                kP_Ep.mean      = posterior_Ep.muPhi(2);
                kEm.mean      = posterior_Em.muPhi(3);
                kEp.mean      = posterior_Ep.muPhi(3);
                kFEm.mean      = posterior_Em.muPhi(4);
                kFEp.mean      = posterior_Ep.muPhi(4);
                % extract variance
                kR_Em.sigma      = posterior_Em.SigmaPhi(1);
                kR_Ep.sigma      = posterior_Ep.SigmaPhi(1);
                kP_Em.sigma      = posterior_Em.SigmaPhi(2);
                kP_Ep.sigma      = posterior_Ep.SigmaPhi(2);
                kEm.sigma      = posterior_Em.SigmaPhi(3);
                kEp.sigma      = posterior_Ep.SigmaPhi(3);
                kFEm.sigma      = posterior_Em.SigmaPhi(4);
                kFEp.sigma      = posterior_Ep.SigmaPhi(4);
        end
end % pool both effort tasks or not?