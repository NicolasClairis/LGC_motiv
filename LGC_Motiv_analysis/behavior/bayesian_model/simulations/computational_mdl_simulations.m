function[r_identif, r_recover, condition, NS] = computational_mdl_simulations(mdl_n, priors_type)
% [r_identif, r_recover, condition, NS] = computational_mdl_simulations(mdl_n, priors_type)
% script for simulating computational model in order to extract confusion
% matrix (i.e. autocorrelation matrix for measuring the identifiability of
% the parameters) and a recoverability matrix (to test if simulated
% parameters can be well recovered).
% Note that this script is based on the mean posteriors actually extracted
% from the subjects of the study => some cases may not be ready yet
% depending on which subjects, sessions, model, you want to use and you
% may need to extract those values first to make the script work.
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
%
% condition: indication of which condition was selected in input
%
% NS: number of subjects included

%% define subjects on which the simulations will be based
if priors_type == 1
    [study_nm, condition, ~, subject_id, NS] = sub_id;
end

%% define model
[mdl_n, mdl_n_nm] = which_bayesian_mdl_n(mdl_n);

%% extract model characteristics
[mdl_prm] = computational_mdl_prm(mdl_n);
n_F_prm = mdl_prm.n_F_prm;
n_G_prm = mdl_prm.n_G_prm;
number_of_parameters = n_F_prm + n_G_prm;

% check if there are F parameters for evolution function
if n_F_prm == 0
    % f function empty
    theta = [];
    f_fname = [];
    % no hidden states
    x0 = [];
    n_hiddenStates = length(x0);
else
    error(['simulations not ready for model with ',num2str(n_F_prm),' yet.']);
end

% define g observation function
g_fname = @g_observation_mdl;

fig_disp = true; % display figure
use_30ksim = true; % ?
prior_type = 'adjusted'; %'adjusted'/'basic'
nb_simulations = 5000; % number of simulations to perform for each model considered

% load variables used by the model
var_saved = getfield(load('var.mat'),'var_saved'); % average inputs utilized for actual participants
% create function here to recreate var here
ratio = getfield(load('M_ratio_for_Simu.mat'),'ratio'); % ?

% define main parameters: number of trials, number of conditions and
% parameters to simulate
n_trialsPerRun = 54;
n_runs = 4;
n_totalTrials = n_trialsPerRun*n_runs;

% 1/sd noise for f_evolution and g_obs. Inf adds 0 noise
alpha = Inf;
sigma = Inf;

% extract mean and SD of parameters across subjects
switch priors_type
    case 0 % uniformed priors
        muPhi = zeros(1, n_G_prm);
        sigmaPhi = eye(n_G_prm);
    case 1 % informed priors => load data
        [muPhi, sigmaPhi_from_data] = avg_g_observation_mdl_prm(mdl_n, condition, NS);
end

% create normally distributed gaussian values supposed to be our participants real sensitivities.
for i_param = 1:n_G_prm
    % version with sd taken from distribution of our dataset
    %     parameterPhi(i_param,:) = randn(nb_run,1).*(sigmaPhi_from_data(i_param)) + muPhi(i_param);
    
    % version with uninformed prior
    parameterPhi(i_param,:) = randn(nb_simulations,1).*(sigmaPhi(i_param,i_param)*3) + muPhi(i_param);
end

posterior_session_run = NaN(nb_simulations, n_G_prm);

%% perform simulations
for i_simu = 1:nb_simulations
    
    %% load input for model with var
    %     DeltaI, avg over participant deltaE         RP_idx            Ep_or_Em          trial_idx               E_chosen            SUm_Echosen
    var = [mean(var_saved(1,:,:),3,'omitnan'); var_saved(2,:,1); var_saved(3,:,1);...
        var_saved(4,:,1); var_saved(5,:,1); mean(var_saved(6,:,:),3,'omitnan');...
        mean(var_saved(7,:,:),3,'omitnan');...
        %       SumAUC+SumNb_answer avg     SumAUC+ratio avg
        mean(var_saved(8,:,:),3); mean(var_saved(9,:,:),3,'omitnan')];
    
    options.sources.type = 0; % binary data
    options.DisplayWin  = 0; % display figure during inversion
    options.verbose     = 0; % display text during inversion (1) or not (0)
    % define number of parameters to estimate
    options.dim = struct('n', n_hiddenStates,... % number of hidden states
        'n_t', n_totalTrials,... % number of trials across runs
        'n_theta',n_F_prm,... % number of evolution parameters
        'n_phi',n_G_prm); % number of observation parameters
    
    %% Launch VBA simulations
    for iRun = 1:n_runs
        sum_AUC = 0;
        m_ratio = 0;
        for trial_i = 1:54
            
            if  var(4,trial_i+54*(iRun-1)) == 1
                var_i = [var(1:8,trial_i+54*(iRun-1));sum_AUC];
            else
                var_i = [var(1:8,trial_i+54*(iRun-1));m_ratio];
            end
            save_var_i(trial_i+54*(iRun-1),i_simu) = var_i(9);
            options.dim.n_t = 1;
            [choiceLeftOptionSimu(trial_i+54*(iRun-1),i_simu), ~, ~, ~, ~, ~] = VBA_simulate(1, f_fname, g_fname,...
                theta, parameterPhi(:,i_simu), var_i, alpha, sigma, options, x0);
            p(trial_i+54*(iRun-1),i_simu) = choiceLeftOptionSimu(trial_i+54*(iRun-1),i_simu);
            % create a thousand values, using the reported proba for the number of 1 in the vector
            proba_vector = zeros(10000,1);
            proba_vector(1:round(choiceLeftOptionSimu(trial_i+54*(iRun-1),i_simu)*10000),1) = 1;
            % shuffle it, and use the first value as the choice picked by a participant
            random_idx = randperm(10000);
            proba_vector_random = proba_vector(random_idx);
            
            if  var(4,trial_i+54*(iRun-1)) == 1
                if proba_vector_random(1) == 0
                    sum_AUC = sum_AUC+44.5/1000;
                elseif proba_vector_random(1) == 1 && var(2,trial_i+54*(iRun-1)) == 1
                    sum_AUC = sum_AUC+100/1000;
                elseif proba_vector_random(1) == 1 && var(2,trial_i+54*(iRun-1)) == 2
                    sum_AUC = sum_AUC+186/1000;
                elseif proba_vector_random(1) == 1 && var(2,trial_i+54*(iRun-1)) == 3
                    sum_AUC = sum_AUC+263/1000;
                end
            else
                if trial_i <= 9
                    if proba_vector_random(1) == 0
                        m_ratio =   ratio(1,1);
                    elseif proba_vector_random(1) == 1 && var(2,trial_i+54*(iRun-1)) == 1
                        m_ratio =   ratio(1,2);
                    elseif proba_vector_random(1) == 1 && var(2,trial_i+54*(iRun-1)) == 2
                        m_ratio =   ratio(1,3);
                    elseif proba_vector_random(1) == 1 && var(2,trial_i+54*(iRun-1)) == 3
                        m_ratio =   ratio(1,4);
                    end
                elseif trial_i > 9 && trial_i <= 18
                    if proba_vector_random(1) == 0
                        m_ratio =   ratio(2,1);
                    elseif proba_vector_random(1) == 1 && var(2,trial_i+54*(iRun-1)) == 1
                        m_ratio =   ratio(2,2);
                    elseif proba_vector_random(1) == 1 && var(2,trial_i+54*(iRun-1)) == 2
                        m_ratio =   ratio(2,3);
                    elseif proba_vector_random(1) == 1 && var(2,trial_i+54*(iRun-1)) == 3
                        m_ratio =   ratio(2,4);
                    end
                elseif trial_i > 18 && trial_i <= 27
                    if proba_vector_random(1) == 0
                        m_ratio =   ratio(3,1);
                    elseif proba_vector_random(1) == 1 && var(2,trial_i+54*(iRun-1)) == 1
                        m_ratio =   ratio(3,2);
                    elseif proba_vector_random(1) == 1 && var(2,trial_i+54*(iRun-1)) == 2
                        m_ratio =   ratio(3,3);
                    elseif proba_vector_random(1) == 1 && var(2,trial_i+54*(iRun-1)) == 3
                        m_ratio =   ratio(3,4);
                    end
                elseif trial_i > 27 && trial_i <= 36
                    if proba_vector_random(1) == 0
                        m_ratio =   ratio(4,1);
                    elseif proba_vector_random(1) == 1 && var(2,trial_i+54*(iRun-1)) == 1
                        m_ratio =   ratio(4,2);
                    elseif proba_vector_random(1) == 1 && var(2,trial_i+54*(iRun-1)) == 2
                        m_ratio =   ratio(4,3);
                    elseif proba_vector_random(1) == 1 && var(2,trial_i+54*(iRun-1)) == 3
                        m_ratio =   ratio(4,4);
                    end
                elseif trial_i > 36 && trial_i <= 45
                    if proba_vector_random(1) == 0
                        m_ratio =   ratio(5,1);
                    elseif proba_vector_random(1) == 1 && var(2,trial_i+54*(iRun-1)) == 1
                        m_ratio =   ratio(5,2);
                    elseif proba_vector_random(1) == 1 && var(2,trial_i+54*(iRun-1)) == 2
                        m_ratio =   ratio(5,3);
                    elseif proba_vector_random(1) == 1 && var(2,trial_i+54*(iRun-1)) == 3
                        m_ratio =   ratio(5,4);
                    end
                elseif trial_i > 45 && trial_i <= 54
                    if proba_vector_random(1) == 0
                        m_ratio =   ratio(6,1);
                    elseif proba_vector_random(1) == 1 && var(2,trial_i+54*(iRun-1)) == 1
                        m_ratio =   ratio(6,2);
                    elseif proba_vector_random(1) == 1 && var(2,trial_i+54*(iRun-1)) == 2
                        m_ratio =   ratio(6,3);
                    elseif proba_vector_random(1) == 1 && var(2,trial_i+54*(iRun-1)) == 3
                        m_ratio =   ratio(6,4);
                    end
                end
            end
        end
    end
    
    %         0, 0.25 0.75 1
    for n = 1:size(choiceLeftOptionSimu,1)
        if 0 <= choiceLeftOptionSimu(n,i_simu) && 0.25 > choiceLeftOptionSimu(n,i_simu)
            choiceLeftOptionSimu(n,i_simu) = 0;
        elseif 0.25 <= choiceLeftOptionSimu(n,i_simu) && 0.5 > choiceLeftOptionSimu(n,i_simu)
            choiceLeftOptionSimu(n,i_simu) = 0.25;
        elseif 0.5 <= choiceLeftOptionSimu(n,i_simu) && 0.75 >= choiceLeftOptionSimu(n,i_simu)
            choiceLeftOptionSimu(n,i_simu) = 0.75;
        elseif 0.75 < choiceLeftOptionSimu(n,i_simu) && 1 >= choiceLeftOptionSimu(n,i_simu)
            choiceLeftOptionSimu(n,i_simu) = 1;
        end
    end
    %             figure()
    %              hist(choiceLeftOptionSimu)
    %% Launch posterior computation
    
    switch prior_type
        case 'adjusted'
            % for posterior calculation, need priors on our k sensitivities.
            % add it only now as it might impact the simu otherwise.
            options.priors.muPhi = zeros(n_G_prm, 1); %(moyenne des priors sur les paramètres)
            options.priors.SigmaPhi = zeros(n_G_prm, n_G_prm);
            options.priors.SigmaPhi(1,1) = sigmaPhi(1,1)*100;
            options.priors.SigmaPhi(2,2) = sigmaPhi(2,2)*100;
            options.priors.SigmaPhi(3,3) = sigmaPhi(3,3)*100;
            options.priors.SigmaPhi(4,4) = sigmaPhi(4,4)*100;
            options.priors.SigmaPhi(5,5) = sigmaPhi(5,5)*100;
            options.priors.SigmaPhi(6,6) = sigmaPhi(6,6)*100;
            options.priors.SigmaPhi(7,7) = sigmaPhi(7,7)*100;
            
        case 'basic'
            % for posterior calculation, need priors on our k sensitivities.
            % add it only now as it might impact the simu otherwise.
            options.priors.muPhi = zeros(n_G_prm, 1); %(moyenne des priors sur les paramètres)
            
            options.priors.SigmaPhi = eye(n_G_prm); % (matrice de covariance des paramètres)
    end
    options.dim.n_t = 216;
    % compute posterior
    [posterior,out] = VBA_NLStateSpaceModel(choiceLeftOptionSimu(:,i_simu)', var, f_fname, g_fname, options.dim, options);
    
    %% save the data
    run_nm = ['run_nb_',num2str(i_simu)];
    toBeSaved.(run_nm).n_trials = n_totalTrials;
    toBeSaved.(run_nm).theta = theta;
    toBeSaved.(run_nm).phi = parameterPhi;
    toBeSaved.(run_nm).u = var;
    toBeSaved.(run_nm).simulated_data_left_opt = choiceLeftOptionSimu;
    toBeSaved.(run_nm).posterior = posterior;
    toBeSaved.(run_nm).out = out;
    toBeSaved.(run_nm).options = options;
    
    
    close all
    
    posterior_session_run(i_simu,:) = posterior;
end % simulation loop


%%
fig;
plot(mean(p,2,'omitnan'))
% plot to check the decrease over time in physical trials, and increase in mental trials.
figure()
plot([nanmean(nanmean(p(1:10,:))),nanmean(nanmean(p(10:20,:))),nanmean(nanmean(p(20:30,:))),nanmean(nanmean(p(30:40,:))),nanmean(nanmean(p(40:54,:))),...
    nanmean(nanmean(p(55:65,:))),nanmean(nanmean(p(65:75,:))),nanmean(nanmean(p(75:85,:))),nanmean(nanmean(p(85:95,:))),nanmean(nanmean(p(95:108,:)))])


% extract posterior for each participant
for j = 1:nb_simulations
    data_tmp = posterior_session_run(j).muPhi';
    param_per_run(j,:) = data_tmp ;
end
param_simu.fake_parameterPhi = parameterPhi';
param_simu.true_parameterPhi = param_per_run;
%         simu_name = strcat('non_specific_param_simu',num2str(i_simu),'.mat');
%         save(simu_name,'param_simu')

%% prepare all the results from the simulations. model
if use_30ksim == true
    parameterPhi = [];
    param_per_run = [];
    for i = 1:6
        param_simu = getfield(load(strcat('non_specific_param_simu',num2str(i),'.mat')),'param_simu');
        parameterPhi = [parameterPhi; param_simu.fake_parameterPhi];
        param_per_run = [param_per_run; param_simu.true_parameterPhi];
    end
end

%% plot results
if fig_disp == true
    
    %% main figure with recoverability and identifiability matrices subplots
    fig;
    
    %% plot correlations between simulated and recovered parameters (recoverability matrix)
    subplot(1,2,1);
    for n = 1:number_of_parameters
        for m = 1:number_of_parameters
            corrmat_real_tmp = corrcoef(param_per_run(:,n),parameterPhi(:,m));
            autocorrmat_tmp = corrcoef(param_per_run(:,n),param_per_run(:,m));
            corrmat(n,m) = corrmat_real_tmp(1,2);
            autocorrmat(n,m) = autocorrmat_tmp(1,2);
        end
    end
    
    % correlation real-estimated plot
    imagesc(corrmat,[-1 1])
    set(gca,'XTick',linspace(1,n_G_prm,n_G_prm))
    set(gca,'XTickLabel',{'kR','kP','kEp','kEm','bias','kFp','kLm'})
    set(gca,'YTick',linspace(1,n_G_prm,n_G_prm))
    set(gca,'YTickLabel',{'kR','kP','kEp','kEm','bias','kFp','kLm'})
    ax = gca;
    ax.FontSize = 16;
    ylabel('Fitted on real data')
    xlabel('Fitted on simulated data')
    title('Model recovery')
    hcb = colorbar;
    set(get(hcb,'label'),'string','Correlation coefficient (R)','Rotation',90.0,'FontSize',16);
    colormap(redblue(45))
    
    %% autocorrelation plot
    subplot(1,2,2)
    % put to 0 the center line for the figure for the paper
    autocorrmat(1,1) =0;autocorrmat(2,2) =0;autocorrmat(3,3) =0;autocorrmat(4,4) =0;autocorrmat(5,5) =0;autocorrmat(6,6) =0;autocorrmat(7,7) =0;
    imagesc(autocorrmat,[-1 1])
    set(gca,'XTick',linspace(1,n_G_prm,n_G_prm))
    set(gca,'XTickLabel',{'kR','kP','kEp','kEm','bias','kFp','kLm'})
    set(gca,'YTick',linspace(1,n_G_prm,n_G_prm))
    set(gca,'YTickLabel',{'kR','kP','kEp','kEm','bias','kFp','kLm'})
    ax = gca;
    ax.FontSize = 16;
    ylabel('Fitted on simulated data')
    xlabel('Fitted on simulated data')
    title('Confusion matrix')
    hcb = colorbar;
    set(get(hcb,'label'),'string','Correlation coefficient (R)','Rotation',90.0,'FontSize',16);
    colormap(redblue(45))
    mean([corrmat(1,1),corrmat(2,2),corrmat(3,3),corrmat(4,4),corrmat(5,5),corrmat(6,6),corrmat(7,7)])
    
end % show plot