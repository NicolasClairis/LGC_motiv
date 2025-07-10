% function[r_identif, r_recover] = computational_mdl_simulations(mdl_n, priors_type)
% [r_identif, r_recover] = computational_mdl_simulations(mdl_n, priors_type)
% script for simulating computational model in order to extract confusion
% matrix (i.e. autocorrelation matrix for measuring the identifiability of
% the parameters) and a recoverability matrix (to test if simulated
% parameters can be well recovered).
% Note that this script is based on the mean posteriors actually extracted
% from the subjects of the study => some cases may not be ready yet
% depending on which subjects, sessions, model, you want to use and you
% may need to extract those values first to make the script work.
%
% Script unfolds in 4 steps:
% 1) Preparing inputs of the model based on the real data
% 2) simulating data based on these inputs and parameter priors
% 3) applying the model on the simulated data to derive posteriors
% 4) check the identifiability and recoverability of the derived posteriors
%
% INPUTS
% mdl_n: model number (will be asked if left empty)
%
% priors_type: use informed priors based on the actual data (1) or
% uninformed priors (0)
%
% OUTPUTS
% r_identif: structure with correlation coefficients (and associated
% p.values) for the auto-correlation between the recovered parameters (to
% perform confusion matrix, aka identifiability matrix)
%
% r_recover: structure with correlation coefficients (and associated
% p.values) for the correlation between the simulated and the recovered
% parameters

%% general parameters
fig_disp = true; % display figure
saveFolder = fullfile('P:','boulot','postdoc_CarmenSandi','results','behavior','simulations');

%% define main parameters
nTrialsPerRun = 54;
nRuns = 4;
nTotalTrials = nTrialsPerRun*nRuns;

%% 1) prepare model inputs based on the real data
disp('(1/4) **Load real data and average it**');
% extract average values for all the input variables (deltaR, deltaP,
% deltaE, etc.)
[~, mean_var, var_names, ~,~,m_AUC_perEch, m_Eff_perEch] = load_computational_mdl_input_across_subs();
% keep everything except for Fp/Lm where we will adjust the data to
% simulated choices
vars = mean_var;
vars(ismember(var_names,{'Fp','currEff','prevEff'}),:) = NaN;

%% 2) simulate data based on inputs, model selected and parameter priors
disp('(2/4) **Start performing simulations**');
%% define general simulation parameters
n_simulations = 10000; % number of simulations to perform for each model considered

% replicate var by number of simulations in 3rd dimension
vars = repmat(vars,1,1,n_simulations);

% noise parameters
% 1/sd noise for f_evolution and g_obs. Inf adds 0 noise
alpha = Inf;
sigma = Inf;

% uninformed priors by default if variable is left empty
if ~exist('priors_type','var') || isempty(priors_type) || ~ismember(priors_type,[0,1])
    priors_type = 0;
end

% define subjects on which the simulations will be based
if priors_type == 1
    [study_nm, condition, ~, subject_id, NS] = sub_id;
end

%% define which model you want to use
mdl_n = 5; % by default use model 5, but consider changing if you want to also test other models
[mdl_n, mdl_n_nm] = which_bayesian_mdl_n(mdl_n);

%% extract model characteristics
[mdl_prm] = computational_mdl_prm(mdl_n);
n_F_prm = mdl_prm.n_F_prm;
G_prm_names = mdl_prm.G_prm_names;
n_G_prm = mdl_prm.n_G_prm;
number_of_parameters = n_F_prm + n_G_prm;

% apply multisession (ie VBA will allow different noise values in the
% estimation of the parameters for each block considering that some blocks
% may be more trustful than others for the estimation of the parameters)
is_multisession = true;

% extract info regarding the way choices are to be modeled
binary_answers = mdl_prm.binary_answers;

%% define bayesian evolution and observation functions
% no evolution function with current models
if n_F_prm == 0
    x0 = [];
    n_hiddenStates = length(x0);
    f_evol_function = [];
    theta = [];
else
    error(['simulations not ready for model with ',num2str(n_F_prm),' yet.']);
end
g_obs_function = @g_observation_mdl;

% define priors for parameters of interest
switch priors_type
    case 0 % uniformed priors
        muPhi = zeros(1, n_G_prm);
        sigmaPhi = eye(n_G_prm);
    case 1 % informed priors => load data
        % extract mean and SD of parameters across subjects
        [muPhi, sigmaPhi_from_data] = avg_g_observation_mdl_prm(mdl_n, condition, NS);
end

% transform parameters into normally distributed gaussian values
[parameterPhi, parameterPhi_transfo,...
    posteriorPhi, posteriorPhi_transfo] = deal(NaN(n_G_prm, n_simulations));
for iPrm = 1:n_G_prm
    G_prm_nm = G_prm_names{iPrm};
    switch priors_type
        case 0
            % version with uninformed prior
            parameterPhi(iPrm,:) = randn(n_simulations,1).*(sigmaPhi(iPrm,iPrm)*3) + muPhi(iPrm);
        case 1
            % version with sd taken from distribution of our dataset
            parameterPhi(iPrm,:) = randn(nRuns,1).*(sigmaPhi_from_data(iPrm)) + muPhi(iPrm);
    end % prior type (uninformed or based on actual data)
    
    % transform corresponding parameter with positivity constraint (just
    % ignore it for the parameters which are not transformed)
    switch mdl_prm.pos.(G_prm_nm)
        case true % positivity constraint => adapt parameter accordingly
            parameterPhi_transfo(iPrm,:) = log(1 + exp(parameterPhi(iPrm,:)));
        case false % no constraint
            parameterPhi_transfo(iPrm,:) = parameterPhi(iPrm,:);
    end
end % parameter loop

%% perform simulations
[pChoice_HE_simu, choice_HE_simu] = deal(NaN(nTotalTrials, n_simulations));
for i_simu = 1:n_simulations
    % alternate between Ep/Em/Ep/Em order and Em/Ep/Em/Ep order for
    % simulations to mimick actual subjects
    switch mod(i_simu,2)
        case 0 % pair Ep/Em order
            mean_var(strcmp(var_names,'EpEm'),:) = [ones(1,nTrialsPerRun), zeros(1,nTrialsPerRun),ones(1,nTrialsPerRun), zeros(1,nTrialsPerRun)];
        otherwise % impair Em/Ep order
            mean_var(strcmp(var_names,'EpEm'),:) = [zeros(1,nTrialsPerRun), ones(1,nTrialsPerRun),zeros(1,nTrialsPerRun), ones(1,nTrialsPerRun)];
    end
    
    %% load input for model
    options = struct;
    switch binary_answers
        case false
            options.sources.type = 0; % 4-levels choices (0/0.25/0.75/1: low E + high Conf/low E + low Conf/high E + low Conf/high E + high Conf)
        case true
            options.sources.type = 1; % binary choices (0/1: low/high E choice)
    end
    options.DisplayWin  = 0; % display figure during inversion
    options.verbose     = 0; % display text during inversion (1) or not (0)
    % define number of parameters to estimate
    options.dim = struct('n', n_hiddenStates,... % number of hidden states
        'n_t', nTotalTrials,... % number of trials across runs
        'n_theta',n_F_prm,... % number of evolution parameters
        'n_phi',n_G_prm); % number of observation parameters
    options.dim.n_t = 1; % 1 trial/simulation
    % information about model parameters (positivity constraint, etc.)
    options.inG.mdl_prm = mdl_prm;
    options.inG.var_names = var_names;
    
    %% Launch VBA simulations
    for iRun = 1:nRuns
        % initialize fatigue/previous efficiency at 0
        Ep_sum_AUC = 0;
        Em_prevEff = 0;
        
        % loop over trials with one simulation/trial
        for iPrm = 1:nTrialsPerRun
            jTrial = iPrm + nTrialsPerRun*(iRun - 1);
            is_Ep_trial = vars(strcmp(var_names,'EpEm'),jTrial,i_simu);
            if  is_Ep_trial == 1 % Ep trial
                vars(strcmp(var_names,'Fp'),jTrial,i_simu) = Ep_sum_AUC;
                vars(strcmp(var_names,'prevEff'),jTrial,i_simu) = 0;
            else % Em trial
                vars(strcmp(var_names,'Fp'),jTrial,i_simu) = 0;
                vars(strcmp(var_names,'prevEff'),jTrial,i_simu) = Em_prevEff;
            end
            
            % launch the simulation trial/trial which allows to extract
            % physical fatigue and mental performance adjusted to simulated
            % choice instead of computing it a priori. The values are
            % nevertheless based on the average Fp/Lm values across
            % subjects for each condition
            [pChoice_HE_simu(jTrial,i_simu)] = VBA_simulate(1, f_evol_function, g_obs_function,...
                theta, parameterPhi(:,i_simu), vars(:,jTrial,i_simu), alpha, sigma, options, x0);
            
            % extract level of effort chosen
            if pChoice_HE_simu(jTrial, i_simu) < 0.5
                Ech_tmp = 0;
            else
                Ech_tmp = vars(strcmp(var_names,'dE'),jTrial,i_simu);
            end
            
            % update accumulated effort/efficiency parameter depending on
            % simulated choice and effort level
            switch is_Ep_trial
                case 1 % Ep trial
                    Ep_sum_AUC = Ep_sum_AUC + m_AUC_perEch(Ech_tmp+1)/1000;
                otherwise % Em trial
                    Em_prevEff = m_Eff_perEch(Ech_tmp+1);
            end
        end % trial loop
    end % run loop
    
    %% transform continuous p(high E choice) variable into 4-levels:
    % 0/0.25/0.75/1 to fit the real data
    for iPrm = 1:nTotalTrials
        if 0 <= pChoice_HE_simu(iPrm,i_simu) && pChoice_HE_simu(iPrm,i_simu) < 0.25
            choice_HE_simu(iPrm,i_simu) = 0;
        elseif 0.25 <= pChoice_HE_simu(iPrm,i_simu) && pChoice_HE_simu(iPrm,i_simu) < 0.5
            choice_HE_simu(iPrm,i_simu) = 0.25;
        elseif 0.5 <= pChoice_HE_simu(iPrm,i_simu) && pChoice_HE_simu(iPrm,i_simu) < 0.75
            choice_HE_simu(iPrm,i_simu) = 0.75;
        elseif 0.75 <= pChoice_HE_simu(iPrm,i_simu) && pChoice_HE_simu(iPrm,i_simu) <= 1
            choice_HE_simu(iPrm,i_simu) = 1;
        end
    end % trial loop
    
    %% indicate where we are
    disp(['Simulation ',num2str(i_simu),'/',num2str(n_simulations),' - done']);
end % simulation loop

%% 3) Compute the model on the simulated data to obtain the posteriors
% Match as much as possible the approach to what is done in
% computational_mdl.m
disp('(3/4) **Start model estimation on simulated data**');
for i_simu = 1:n_simulations
    % for posterior calculation, need priors on our k sensitivities.
    % add it only now as it might impact the simu otherwise.
    options.priors.muPhi = zeros(n_G_prm, 1); % moyenne des priors sur les paramètres
    options.priors.SigmaPhi = eye(n_G_prm)*100; % matrice de covariance des paramètres
    options.dim.n_t = nTotalTrials;
    options.MinIter = 3; % minimum number of VB iterations {1} by default
    options.TolFun  = 1e-7 ; % minimum absolute increase of the free energy {2e-2} by default
    % multisession
    switch is_multisession
        case true
            options.multisession.split = repmat(nTrialsPerRun,1,nRuns); % split in 4 equal sessions
            % fix the parameters to be equal across sessions
            options.multisession.fixed.phi = 1:n_G_prm;
        otherwise
            error(['is_multisession = ',num2str(is_multisession),' not ready yet.']);
    end
    % information about model parameters (positivity constraint, etc.)
    options.inG.mdl_prm = mdl_prm;
    options.inG.var_names = var_names;
    
    % compute posterior
    [posterior, out] = VBA_NLStateSpaceModel(pChoice_HE_simu(:,i_simu)', vars(:,:,i_simu), f_evol_function, g_obs_function, options.dim, options);
    
    %% save the data
    posteriorPhi(:, i_simu) = posterior.muPhi;
    for iPrm = 1:n_G_prm
        G_prm_nm = G_prm_names{iPrm};
        switch mdl_prm.pos.(G_prm_nm)
            case true % positivity constraint => adapt parameter accordingly
                posteriorPhi_transfo(iPrm,i_simu) = log(1 + exp(posteriorPhi(iPrm,i_simu)));
            case false % no constraint
                posteriorPhi_transfo(iPrm,i_simu) = posteriorPhi(iPrm,i_simu);
        end
    end % loop on parameters
    
    %% indicate where we are
    disp(['Model estimation on simulation ',num2str(i_simu),'/',num2str(n_simulations),' - done']);
end % simulation loop


%% 4) compare obtained posteriors to priors (recoverability matrix) and posteriors to themselves (identifiability/confusion matrix)
disp('(4/4) **Start computing identifiability and recoverability matrices**');
% recoverability matrix
[r_recover.r_corr, r_recover.pval] = deal(NaN(1,n_G_prm));
r_recover.var_nm = cell(1,n_G_prm);
for iPrm = 1:n_G_prm
    var_nm = var_names{iPrm};
    [r_recover.r_corr(iPrm), r_recover.pval(iPrm)] = corr(parameterPhi_transfo(iPrm,:)', posteriorPhi_transfo(iPrm,:)');
    r_recover.var_nm{iPrm} = var_nm;
end % loop over parameters
% extract min recovery between simulated and recovered paramters to assess
% recoverability
r_recover.min = min(r_recover.r_corr,[],2,'omitnan');

% identifiability matrix
[r_identif.r_corr, r_identif.pval, corrmat, autocorrmat] = deal(NaN(n_G_prm));
r_identif.var_nm = cell(n_G_prm,n_G_prm);
for iPrm = 1:n_G_prm
    var_nm1 = var_names{iPrm};
    for jPrm = 1:n_G_prm
        var_nm2 = var_names{jPrm};
        [r_identif.r_corr(iPrm, jPrm), r_identif.pval(iPrm, jPrm)] = corr(posteriorPhi_transfo(iPrm,:)', posteriorPhi_transfo(jPrm,:)');
        r_identif.var_nm{iPrm, jPrm} = [var_nm1,'_vs_',var_nm2];
        
        % old extraction by Arthur still useful for figures
        corrmat_real_tmp = corrcoef(parameterPhi_transfo(iPrm,:)', posteriorPhi_transfo(jPrm,:)');
        autocorrmat_tmp = corrcoef(posteriorPhi_transfo(iPrm,:)', posteriorPhi_transfo(jPrm,:)');
        corrmat(iPrm,jPrm) = corrmat_real_tmp(1,2);
        autocorrmat(iPrm,jPrm) = autocorrmat_tmp(1,2);
    end % loop over parameters
end % loop over parameters
% extract max confusion correlation between independent variables to assess
% identifiability
max_perLine = NaN(n_G_prm-1,1);
for iG_prm = 2:n_G_prm
   max_perLine(iG_prm - 1) = max(r_identif.r_corr(iG_prm,1:(iG_prm-1)),[],2,'omitnan'); 
end
r_identif.max_r_corr = max(max_perLine,[],1,'omitnan');
% save the data
save([saveFolder,filesep,...
    'Simulations_model_',mdl_n_nm,'_',num2str(n_simulations),'simus.mat']);
% indicate where we are to the user
disp('Computation finished, now script will eventually display correlation matrices');

%% figure display
if fig_disp == true
    fig;
    plot(mean(pChoice_HE_simu,2,'omitnan'));
    
    %% main figure with recoverability and identifiability matrices subplots
    fig;
    
    %% plot correlations between simulated and recovered parameters
    % (recoverability matrix)
    subplot(1,2,1);
    % correlation real-estimated plot
    imagesc(corrmat,[-1 1]);
    set(gca,'XTick',linspace(1,n_G_prm,n_G_prm));
    set(gca,'XTickLabel',G_prm_names);
    set(gca,'YTick',linspace(1,n_G_prm,n_G_prm));
    set(gca,'YTickLabel',G_prm_names);
    ax = gca;
    ax.FontSize = 16;
    xlabel('Simulated');
    ylabel('Recovered');
    title('Model recovery');
    hcb = colorbar;
    set(get(hcb,'label'),'string','Correlation coefficient (r)','Rotation',90.0,'FontSize',16);
    colormap(redblue(45));
    
    %% autocorrelation plot
    subplot(1,2,2);
    % put to 0 the center line for the figure for the paper
    autocorrmat(1,1) = 0;
    autocorrmat(2,2) = 0;
    autocorrmat(3,3) = 0;
    autocorrmat(4,4) = 0;
    autocorrmat(5,5) = 0;
    autocorrmat(6,6) = 0;
    autocorrmat(7,7) = 0;
    imagesc(autocorrmat,[-1 1]);
    set(gca,'XTick',linspace(1,n_G_prm,n_G_prm));
    set(gca,'XTickLabel',G_prm_names);
    set(gca,'YTick',linspace(1,n_G_prm,n_G_prm));
    set(gca,'YTickLabel',G_prm_names);
    ax = gca;
    ax.FontSize = 16;
    ylabel('Recovered');
    xlabel('Recovered');
    title('Model identifiability');
    hcb = colorbar;
    set(get(hcb,'label'),'string','Correlation coefficient (r)','Rotation',90.0,'FontSize',16);
    colormap(redblue(45));
    mean([corrmat(1,1),corrmat(2,2),corrmat(3,3),corrmat(4,4),corrmat(5,5),corrmat(6,6),corrmat(7,7)]);
    
end % show plot

% end % function