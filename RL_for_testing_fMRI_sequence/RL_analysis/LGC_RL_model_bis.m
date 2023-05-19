function [  ] = LGC_RL_model_bis( RL_model_n )
%[  ] = LGC_RL_model_bis( RL_model_n )
%MS2_RL_model_bis similar to MS2_RL_Fabien_version.m will compute a model
%for the Q.learning in RL version but with a greater focus on uncertainty.
%
% INPUTS
% RL_model_n: model number
%
% See also MS2_RL_model_define.m for model details

%% install the VBA toolbox if not present
wherever = pwd;
if ~exist('RFT_GLM_contrast.m','file')
    VBA_folder = fullfile('C:', 'Users', 'Nicolas Clairis',...
        'Documents', 'GitHub', 'VBA-toolbox');
    cd(VBA_folder);
    VBA_setup;
    cd(wherever);
end

%% working directories
root = ['E:',filesep,'study2',filesep,...
    'pilots', filesep];
warning('check path here');

%% subject selection
warning('ADD LIST SUBJECTS HERE')
subject_id = {''};
NS = length(subject_id);

%% load model parameters
[RL_mdl_prm] = MS2_RL_model_define(RL_model_n);
alpha_prm       = RL_mdl_prm.alpha_prm;
sigmaQ_prm      = RL_mdl_prm.sigmaQ_prm;
RP_weights_prm  = RL_mdl_prm.RP_weights;
side_bias_prm   = RL_mdl_prm.side_bias;
opt_bias_prm    = RL_mdl_prm.opt_bias;

% number of hidden states
n_hiddenStates = 4; % QA gain gain/QB neutral gain/QE loss loss/QE loss neutral ignore neutral Q.values (no updating)

% number of F parameters
n_F_prm = 0;
switch alpha_prm
    case 'one_learningRate'
        n_F_prm = n_F_prm + 1; % alpha
    case 'GL_learningRates'
        n_F_prm = n_F_prm + 2; % alpha gain/alpha loss
end
switch RP_weights_prm
    case 1
        n_F_prm = n_F_prm + 2; % R/P modulators of reinforcement
end

% number of G parameters
n_G_prm = 1; % beta temperature
% G_Phi_prm_nm{n_G_prm} = 'inverse_temperature';
G_Phi_prm_nm{n_G_prm} = 'temperature';
switch side_bias_prm
    case 1
        n_G_prm = n_G_prm + 1;
        G_Phi_prm_nm{n_G_prm} = 'side_bias';
end

switch opt_bias_prm
    case 'ntal'
        n_G_prm = n_G_prm + 1;
        G_Phi_prm_nm{n_G_prm} = 'ntal_option_bias';
    case 'all_pairs'
       n_G_prm = n_G_prm + 1;
        G_Phi_prm_nm{n_G_prm} = 'gain_option_bias';
        n_G_prm = n_G_prm + 1;
        G_Phi_prm_nm{n_G_prm} = 'ntal_option_bias';
        n_G_prm = n_G_prm + 1;
        G_Phi_prm_nm{n_G_prm} = 'loss_option_bias';
end

%% main task parameters
n_trials_per_run    = 60;
n_RL_runs           = 3;
n_total_trials      = n_trials_per_run*n_RL_runs;
trialN_vector = 1:n_trials_per_run;
trialN = repmat(trialN_vector, 1, n_RL_runs);

%% prepare variables of interest
[alpha.mean.perSub.aRuns,...
    alpha_G.mean.perSub.aRuns,...
    alpha_L.mean.perSub.aRuns,...
    R_weight.mean.perSub.aRuns,...
    P_weight.mean.perSub.aRuns,...
    beta.mean.perSub.aRuns,...
    side_bias.mean.perSub.aRuns,...
    alpha.sigma.perSub.aRuns,...
    alpha_G.sigma.perSub.aRuns,...
    alpha_L.sigma.perSub.aRuns,...
    R_weight.sigma.perSub.aRuns,...
    P_weight.sigma.perSub.aRuns,...
    beta.sigma.perSub.aRuns,...
    side_bias.sigma.perSub.aRuns] = deal(NaN(1, NS));
[opt_bias.mean.G.perSub.perRun,...
    opt_bias.mean.N.perSub.perRun,...
    opt_bias.mean.L.perSub.perRun,...
    opt_bias.sigma.G.perSub.perRun,...
    opt_bias.sigma.N.perSub.perRun,...
    opt_bias.sigma.L.perSub.perRun] = deal( NaN(n_RL_runs, NS) );

for iRun = 1:n_RL_runs
    run_nm = ['run_',num2str(iRun)];
    [Qvalues.mean.gainPair_best.perSub.(run_nm),...
        Qvalues.mean.gainPair_bad.perSub.(run_nm),...
        Qvalues.mean.lossPair_best.perSub.(run_nm),...
        Qvalues.mean.lossPair_bad.perSub.(run_nm),...
        Qvalues.sigma.gainPair_best.perSub.(run_nm),...
        Qvalues.sigma.gainPair_bad.perSub.(run_nm),...
        Qvalues.sigma.lossPair_best.perSub.(run_nm),...
        Qvalues.sigma.lossPair_bad.perSub.(run_nm),...
        pChoice_best.perSub.(run_nm), pChoice_worse.perSub.(run_nm),...
        uncertainty_sigmaOptions.perSub.(run_nm),...
        uncertainty_overlap_distrib_options.perSub.(run_nm)] = deal(NaN(n_trials_per_run, NS));
end
model_quality = cell(1,NS);

%%

for iS = 1:NS
        
        %% subject identification
        sub_nm = subject_id{iS};
        if strcmp(sub_nm(3),'_')
            subid   = sub_nm(2);
        elseif ~strcmp(sub_nm(3),'_') && strcmp(sub_nm(4),'_')
            subid = sub_nm(2:3);
        end
        
        behavior_folder = [root,sub_nm,filesep,'behavior',filesep];
        onsets_folder = [root,sub_nm,filesep,'fMRI_analysis',filesep];
        
        %% vars of interest
        [choice_best,...
            goodSide,...
            pairValence,...
            lastPairValence,...
            lastChoice,...
            lastSide,...
            lastOutcome] = deal( NaN(1, n_total_trials) ); % trialN
        
        %% pool data across runs
        [runs_idx] = MS2_task_runs_extraction('RL', sub_nm);
        
        % loop through runs
        for iRun = 1:n_RL_runs
            jRun = runs_idx(iRun);
            run_nm = num2str(jRun);
            
            % load data
            loadStruct = load([behavior_folder,'global_sub_',subid,'_session_',run_nm,'_learning.mat']);
            feedback_tmp    = loadStruct.feedback; % -1/1 worse/best possible feedback, 0=missed trial
            gain_tmp        = loadStruct.gain; % -1/0/1 loss/neutral/gain
            npair_tmp       = loadStruct.npair; % 1/2/3 gain/neutral/loss
            valence_tmp = npair_tmp; % valence -1/0/1 loss/neutral/gain
            valence_tmp(npair_tmp == 3) = -1;
            valence_tmp(npair_tmp == 2) = 0;
            choice_LR_tmp   = loadStruct.choice; % -1/1 left/right
            bad_trials = choice_LR_tmp == 0;
            choice_LR_tmp(bad_trials) = NaN; % put NaN when no choice was made
            side_best_tmp   = loadStruct.side; % -1/1 left/right
            choice_best_tmp = choice_LR_tmp == side_best_tmp;
            choice_best_tmp = double(choice_best_tmp); % convert from logical to double
            choice_best_tmp(bad_trials) = NaN;
%             prev_choice_LR_tmp          = [NaN, choice_LR_tmp(1:(end - 1))];
            choice_best_prevTrial_tmp   = [NaN, choice_best_tmp(1:(end - 1))];
            valence_prevTrial_tmp       = [NaN, valence_tmp(1:(end - 1))];
            sideBest_prevTrial_tmp      = [NaN, side_best_tmp(1:(end - 1))];
            gain_prevTrial_tmp          = [NaN, gain_tmp(1:(end - 1))];
            
            good_run_trials = feedback_tmp ~= 0;
            
            trialN_tmp = trialN_vector(good_run_trials); % trial index for the current run after removing bad trials
            trialN_tmp_bis = trialN_tmp + n_trials_per_run*(iRun - 1); % trial index for the across runs variables after removing bad trials
            trialN_tmp_ter = trialN_vector + n_trials_per_run*(iRun - 1); % includes all trials of the run
            
            
            choice_best(trialN_tmp_bis)     = choice_best_tmp(trialN_tmp);
            goodSide(trialN_tmp_ter)        = side_best_tmp;
            pairValence(trialN_tmp_ter)     = valence_tmp;
            lastPairValence(trialN_tmp_ter) = valence_prevTrial_tmp;
            lastChoice(trialN_tmp_ter)      = choice_best_prevTrial_tmp;
            lastSide(trialN_tmp_ter)        = sideBest_prevTrial_tmp;
            lastOutcome(trialN_tmp_ter)     = gain_prevTrial_tmp;
            
        end % run loop
        
        
        %% prepare the model to use for VBA
        y = choice_best;
        u_t = [trialN; pairValence; goodSide; lastPairValence; lastChoice; lastOutcome];
        dim = struct('n', n_hiddenStates,... % number of hidden states: 6 Q.values
            'n_t', n_total_trials,... % number of trials across runs
            'n_theta',n_F_prm,... % number of evolution parameters
            'n_phi',n_G_prm); % number of observation parameters
        
        
        options = struct;
        options.sources.type = 1; % binary data
        options.dim         = dim;
        options.GnFigs      = 0;
        options.DisplayWin  = 0; % display figure during inversion
        options.verbose     = 0; % display text during inversion
        
        % ignore NaN trials (=those where no (clear) answer was provided) for the estimation of the model
        NaN_trials = isnan(y);
        options.isYout = NaN_trials;
        % avoid updating hidden states for the first trial of every session
        % + for all trials coming after a NaN trial
        options.skipf = zeros(n_total_trials,1);
        options.skipf(trialN == 1) = 1;
        postNaN_trials = [0, NaN_trials(1:(end-1))];
        postNaN_trials = postNaN_trials == 1; % reconvert from double to logical
        options.skipf(postNaN_trials) = 1;
        
        % use multisession to pool 3 runs together
        options.multisession.split = repmat(n_trials_per_run, 1, n_RL_runs);
        % fix the parameters (only the hidden states should vary in each
        % session, in theory)
        options.multisession.fixed.theta = 1:n_F_prm; % fix the evolution parameters across sessions (learning rate, R/P weights)
        options.multisession.fixed.phi = 1:n_G_prm; % fix the observation parameters across sessions (temperature and side bias)
        
        % do not fix the option bias across sessions
        if ~strcmp(opt_bias_prm,'')
            % remove neutral option bias (should be independent for every
            % session)
            ntal_opt_bias_idx = find(strcmp(G_Phi_prm_nm, 'ntal_option_bias'));
            prm_to_fix_ntal_idx = options.multisession.fixed.phi == ntal_opt_bias_idx;
            options.multisession.fixed.phi(prm_to_fix_ntal_idx) = [];
            switch opt_bias_prm
                case 'all_pairs'
                    % gain
                    gain_opt_bias_idx = find(strcmp(G_Phi_prm_nm, 'gain_option_bias'));
                    prm_to_fix_gain_idx = options.multisession.fixed.phi == gain_opt_bias_idx;
                    options.multisession.fixed.phi(prm_to_fix_gain_idx) = [];
                    
                    % loss
                    loss_opt_bias_idx = find(strcmp(G_Phi_prm_nm, 'loss_option_bias'));
                    prm_to_fix_loss_idx = options.multisession.fixed.phi == loss_opt_bias_idx;
                    options.multisession.fixed.phi(prm_to_fix_loss_idx) = [];
            end
        end
        
        % information about the model used
        options.inF.RL_mdl_prm = RL_mdl_prm;
        options.inG.RL_mdl_prm = RL_mdl_prm;
        
        %% define priors
        % hidden values priors
        options.priors.muX0 = zeros(n_hiddenStates, 1);
        options.priors.SigmaX0 = eye(n_hiddenStates).*sigmaQ_prm;
        
        % f parameters priors = learning rate(s) + R/P weigths
        options.priors.muTheta = zeros(n_F_prm, 1);
        options.priors.SigmaTheta = eye(n_F_prm);
        
        % g parameters priors = beta inverse temperature + motor bias
        %         options.priors.muPhi = [1, zeros(n_G_prm-1, 1)];
        options.priors.muPhi = zeros(n_G_prm, 1);
        options.priors.SigmaPhi = eye(n_G_prm);
%         % increase a bit the sigma for the inverse temperature
%         options.priors.SigmaPhi(1,1) = 10;
        
        %% Define functions
        f_fname = @LGC_RL_f_evolution;
        g_fname = @LGC_RL_g_evolution;
        
        %% Do inversion
        [posterior,out] = VBA_NLStateSpaceModel(y, u_t, f_fname, g_fname, dim, options);
        
        %% extract variables and split per run
        %% extract evolution parameters
        % extract learning rate
        iTheta = 0;
        switch alpha_prm
            case 'one_learningRate'
                iTheta = iTheta + 1;
                alpha.mean.perSub.aRuns(iS) = posterior.muTheta(iTheta);
                alpha.sigma.perSub.aRuns(iS) = posterior.SigmaTheta(iTheta);
            case 'GL_learningRates'
                iTheta = iTheta + 1;
                alpha_G.mean.perSub.aRuns(iS) = posterior.muTheta(iTheta);
                alpha_G.sigma.perSub.aRuns(iS) = posterior.SigmaTheta(iTheta);
                iTheta = iTheta + 1;
                alpha_L.mean.perSub.aRuns(iS) = posterior.muTheta(iTheta);
                alpha_L.sigma.perSub.aRuns(iS) = posterior.SigmaTheta(iTheta);
        end
        
        % extract R/P weigths
        switch RP_weights_prm
            case 1
                iTheta = iTheta + 1;
                R_weight.mean.perSub.aRuns(iS) = posterior.muTheta(iTheta);
                R_weight.sigma.perSub.aRuns(iS) = posterior.SigmaTheta(iTheta);
                iTheta = iTheta + 1;
                P_weight.mean.perSub.aRuns(iS) = posterior.muTheta(iTheta);
                P_weight.sigma.perSub.aRuns(iS) = posterior.SigmaTheta(iTheta);
        end
        
        %% extract observation parameters
%         % beta inverse temperature
%         iPhi = 1;
%         beta.mean.perSub.aRuns(iS) = posterior.muPhi(iPhi);
%         beta.sigma.perSub.aRuns(iS) = posterior.SigmaPhi(iPhi, iPhi);
        % beta temperature
        iPhi = 1;
        beta.mean.perSub.aRuns(iS) = exp(posterior.muPhi(iPhi));
        beta.sigma.perSub.aRuns(iS) = exp(posterior.SigmaPhi(iPhi, iPhi));

        % side and motor bias
        switch side_bias_prm
            case 1
                iPhi = iPhi + 1;
                side_bias.mean.perSub.aRuns(iS) = posterior.muPhi(iPhi);
                side_bias.sigma.perSub.aRuns(iS) = posterior.SigmaPhi(iPhi, iPhi);
        end
        
        % option bias (not fixed across runs (would be nonsense since pictures are different) => must be extracted for each
        % run separately)
        for iRun = 1:n_RL_runs
            switch opt_bias_prm
                case 'ntal'
                    iPhi = iPhi + 1;
                    opt_bias.mean.N.perSub.perRun(iRun, iS) = posterior.muPhi(iPhi);
                    opt_bias.sigma.N.perSub.perRun(iRun, iS) = posterior.SigmaPhi(iPhi, iPhi);
                case 'all_pairs'
                    iPhi = iPhi + 1;
                    opt_bias.mean.G.perSub.perRun(iRun, iS) = posterior.muPhi(iPhi);
                    opt_bias.sigma.G.perSub.perRun(iRun, iS) = posterior.SigmaPhi(iPhi, iPhi);
                    iPhi = iPhi + 1;
                    opt_bias.mean.N.perSub.perRun(iRun, iS) = posterior.muPhi(iPhi);
                    opt_bias.sigma.N.perSub.perRun(iRun, iS) = posterior.SigmaPhi(iPhi, iPhi);
                    iPhi = iPhi + 1;
                    opt_bias.mean.L.perSub.perRun(iRun, iS) = posterior.muPhi(iPhi);
                    opt_bias.sigma.L.perSub.perRun(iRun, iS) = posterior.SigmaPhi(iPhi, iPhi);
            end
        end
        
        %% extract hidden states and their variance
        % extract Q.values for each option and split them per run
        for iRun = 1:n_RL_runs
            run_nm = ['run_',num2str(iRun)];
            
            Q_GP_best_idx   = 1 + n_hiddenStates*(iRun - 1);
            Q_GP_worse_idx  = 2 + n_hiddenStates*(iRun - 1);
            Q_LP_best_idx   = 3 + n_hiddenStates*(iRun - 1);
            Q_LP_worse_idx  = 4 + n_hiddenStates*(iRun - 1);
            
            run_idx = trialN_vector + n_trials_per_run*(iRun - 1);
            
            % Q.values mean
            Qvalues.mean.gainPair_best.perSub.(run_nm)(:,iS)  = posterior.muX(Q_GP_best_idx, run_idx);
            Qvalues.mean.gainPair_bad.perSub.(run_nm)(:,iS)   = posterior.muX(Q_GP_worse_idx, run_idx);
            Qvalues.mean.lossPair_best.perSub.(run_nm)(:,iS)  = posterior.muX(Q_LP_best_idx, run_idx);
            Qvalues.mean.lossPair_bad.perSub.(run_nm)(:,iS)   = posterior.muX(Q_LP_worse_idx, run_idx);
            
            % Q.values sigma
            sigma_Qvalues = VBA_getVar(posterior.SigmaX.current);
            Qvalues.sigma.gainPair_best.perSub.(run_nm)(:,iS)  = sigma_Qvalues(Q_GP_best_idx, run_idx);
            Qvalues.sigma.gainPair_bad.perSub.(run_nm)(:,iS)   = sigma_Qvalues(Q_GP_worse_idx, run_idx);
            Qvalues.sigma.lossPair_best.perSub.(run_nm)(:,iS)  = sigma_Qvalues(Q_LP_best_idx, run_idx);
            Qvalues.sigma.lossPair_bad.perSub.(run_nm)(:,iS)   = sigma_Qvalues(Q_LP_worse_idx, run_idx);
            
            % p(choice = best) predicted
            pChoice_best.perSub.(run_nm)(:,iS)  = out.suffStat.gx(run_idx);
            pChoice_worse.perSub.(run_nm)(:,iS) = 1 - out.suffStat.gx(run_idx);
        end
        
        %% extract log-likelihood and free energy of the model for model comparison
        model_quality{iS} = out;
        
        %% save individual data
        alpha_sub.mean          = alpha.mean.perSub.aRuns(iS);
        alpha_sub.sigma         = alpha.sigma.perSub.aRuns(iS);
        alpha_G_sub.mean        = alpha_G.mean.perSub.aRuns(iS);
        alpha_G_sub.sigma       = alpha_G.sigma.perSub.aRuns(iS);
        alpha_L_sub.mean        = alpha_L.mean.perSub.aRuns(iS);
        alpha_L_sub.sigma       = alpha_L.sigma.perSub.aRuns(iS);
        R_weight_sub.mean       = R_weight.mean.perSub.aRuns(iS);
        R_weight_sub.sigma      = R_weight.sigma.perSub.aRuns(iS);
        P_weight_sub.mean       = P_weight.mean.perSub.aRuns(iS);
        P_weight_sub.sigma      = P_weight.sigma.perSub.aRuns(iS);
        beta_sub.mean           = beta.mean.perSub.aRuns(iS);
        beta_sub.sigma          = beta.sigma.perSub.aRuns(iS);
        side_bias_sub.mean      = side_bias.mean.perSub.aRuns(iS);
        side_bias_sub.sigma     = side_bias.sigma.perSub.aRuns(iS);
        for iRun = 1:n_RL_runs
            run_nm = ['run_',num2str(iRun)];
            
            % option bias
            opt_bias_G_sub.mean.(run_nm)     = opt_bias.mean.G.perSub.perRun(iRun,iS);
            opt_bias_G_sub.sigma.(run_nm)    = opt_bias.sigma.G.perSub.perRun(iRun,iS);
            opt_bias_N_sub.mean.(run_nm)     = opt_bias.mean.N.perSub.perRun(iRun,iS);
            opt_bias_N_sub.sigma.(run_nm)    = opt_bias.sigma.N.perSub.perRun(iRun,iS);
            opt_bias_L_sub.mean.(run_nm)     = opt_bias.mean.L.perSub.perRun(iRun,iS);
            opt_bias_L_sub.sigma.(run_nm)    = opt_bias.sigma.L.perSub.perRun(iRun,iS);
            
            % Q.values
            Qval_sub.mean.GP_best.(run_nm)   = Qvalues.mean.gainPair_best.perSub.(run_nm)(:,iS);
            Qval_sub.sigma.GP_best.(run_nm)  = Qvalues.sigma.gainPair_best.perSub.(run_nm)(:,iS);
            Qval_sub.mean.GP_worse.(run_nm)  = Qvalues.mean.gainPair_bad.perSub.(run_nm)(:,iS);
            Qval_sub.sigma.GP_worse.(run_nm) = Qvalues.sigma.gainPair_bad.perSub.(run_nm)(:,iS);
            Qval_sub.mean.LP_best.(run_nm)   = Qvalues.mean.lossPair_best.perSub.(run_nm)(:,iS);
            Qval_sub.sigma.LP_best.(run_nm)  = Qvalues.sigma.lossPair_best.perSub.(run_nm)(:,iS);
            Qval_sub.mean.LP_worse.(run_nm)  = Qvalues.mean.lossPair_bad.perSub.(run_nm)(:,iS);
            Qval_sub.sigma.LP_worse.(run_nm) = Qvalues.sigma.lossPair_bad.perSub.(run_nm)(:,iS);
            
            % extract uncertainty variables
            for iT = 1:n_trials_per_run
                jT = iT + n_trials_per_run*(iRun - 1);
                
                if pairValence(jT) ~= 0
                    switch pairValence(jT)
                        case -1
                            mBest       = Qval_sub.mean.LP_best.(run_nm)(iT);
                            sigmaBest   = Qval_sub.sigma.LP_best.(run_nm)(iT);
                            mWorse      = Qval_sub.mean.LP_worse.(run_nm)(iT);
                            sigmaWorse  = Qval_sub.sigma.LP_worse.(run_nm)(iT);
                        case 1
                            mBest       = Qval_sub.mean.GP_best.(run_nm)(iT);
                            sigmaBest   = Qval_sub.sigma.GP_best.(run_nm)(iT);
                            mWorse      = Qval_sub.mean.GP_worse.(run_nm)(iT);
                            sigmaWorse  = Qval_sub.sigma.GP_worse.(run_nm)(iT);
                    end
                    
                    % uncertainty1: sigmaA + sigmaB per trial (depending on pair)
                    uncertainty_sub.sigma_options.(run_nm)(iT) = sigmaBest + sigmaWorse;
                    
                    % uncertainty2: overlap between sigmaA & sigmaB per trial
                    % (depending on pair)
                    if sigmaQ_prm == 0 % ignore this parameter in that case (creates weird bugs)
                        uncertainty_sub.overlap_distrib_options.(run_nm)(iT) = NaN;
                    else
                        uncertainty_sub.overlap_distrib_options.(run_nm)(iT) =...
                            AUC_gaussians_overlap( mBest, sigmaBest, mWorse, sigmaWorse );
                    end
                end
            end
            % save uncertainty variables for pool across subjects
            uncertainty_sigmaOptions.perSub.(run_nm)(:,iS)              = uncertainty_sub.sigma_options.(run_nm);
            uncertainty_overlap_distrib_options.perSub.(run_nm)(:,iS)   = uncertainty_sub.overlap_distrib_options.(run_nm);

            % uncertainty3: p(choice = best) per trial
            uncertainty_sub.pChoice_Worse.(run_nm) = 1 - pChoice_best.perSub.(run_nm)(:,iS);
        end
        save([onsets_folder, 'RLmodel_bis_model',num2str(RL_model_n),'_sub',subid,'.mat'],...
            'alpha_sub','alpha_G_sub','alpha_L_sub','R_weight_sub','P_weight_sub',...
            'beta_sub','side_bias_sub','opt_bias_G_sub','opt_bias_N_sub','opt_bias_L_sub',...
            'Qval_sub',...
            'uncertainty_sub',...
            'out');
end % subject loop

end % function