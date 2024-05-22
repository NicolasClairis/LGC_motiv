% script for simulating first model. where f = [] and everything is in observation

%% clean workspace before startingclearvars;
clearvars;
close all;
clc;

set(0,'defaultfigurecolor',[1 1 1])
%% set run parameters
% WARNING, MODIFY G_OBSERVATION TOO IF YOU MODIFY THESE VARIABLES.
number_of_parameters = 7;
model_nb = 1;% define which model I am currently working on

show_plot = true;
use_30ksim = true;
number_of_different_sessions = 1; % 1-2-3-4
prior_type = 'adjusted'; %'adjusted' 'basic'
nb_run = 5000;
nb_iteration = 1;

%% variables initialization

% place to store results
results_folder = 'D:\Matlab codes\Data simulation\Simulated Data\';
addpath('D:\Matlab codes\Data simulation');
addpath('D:\Matlab codes\Main_analysis_all_dataset_codes');
addpath ('D:\Matlab codes\Data simulation\Extraction Sensitivities');
% load variables used by the model
load('var.mat')
load('M_ratio_for_Simu.mat')
% ratio = [mean(ratio,2) mean(ratio,2) mean(ratio,2) mean(ratio,2)];
% define main parameters: number of trials, number of conditions and
% parameters to simulate
if number_of_different_sessions == 1
    session = 4;
elseif number_of_different_sessions == 2
    session = [2,4];
elseif number_of_different_sessions == 3
    session = [2,4,6];
elseif number_of_different_sessions == 4
    session = [2,4,6,8];
end

n_trialsPerCondition = 27;
n_conditionsPerSession = 2; % reward AND punishment

% number of conditions
n_R_levels = 4;
n_E_levels = 4;

n_trialsPerSession = n_trialsPerCondition*n_conditionsPerSession;

f_fname = [];
g_fname = @g_observation55;

% % define nb of parameters
theta = [];
n_F_prm = length(theta);

% 1/sd noise for f_evolution and g_obs. Inf adds 0 noise
alpha = Inf;
sigma = Inf;

% in case we have f_evaluation, it would be our fatigue param but not here.
% (set initial physical and mental fatigue at zero (but could be higher depending on the baseline))
x0 = [];
n_hiddenStates = length(x0);

% want to have variability as ''real'' datam, use mean of data over participants
muPhi(1) = -0.68;
muPhi(2) = -1.54;
muPhi(3) = 1;
muPhi(4) = 4;
muPhi(5) = 2;
muPhi(6) = -2.63;
muPhi(7) = 4;

% version where you use parameters uninformed on the data, v1 with eff = 10x learning and fatigue
% muPhi(1) = log(exp(1)-1);
% muPhi(2) = log(exp(1)-1);
% muPhi(3) = log(exp(5/1.3)-1);
% muPhi(4) = log(exp(500/92)-1);
% muPhi(5) = 0;
% muPhi(6) = log(exp(5/6)-1);
% muPhi(7) = log(exp(50/92)-1);

% version where you use parameters uninformed on the data, v1 with eff = 2x learning and fatigue
% muPhi(1) = log(exp(1)-1);
% muPhi(2) = log(exp(1)-1);
% muPhi(3) = log(exp(2)-1);
% muPhi(4) = log(exp(25/3)-1);
% muPhi(5) = 0;
% muPhi(6) = log(exp(1)-1);
% muPhi(7) = log(exp(25/6)-1);

% version where you use parameters uninformed on the data, v1 with eff = 4x learning and fatigue
muPhi(1) = log(exp(1)-1);
muPhi(2) = log(exp(1)-1);
muPhi(3) = log(exp(20/7)-1);
muPhi(4) = log(exp(25/4)-1);
muPhi(5) = 0;
muPhi(6) = log(exp(5/7)-1);
muPhi(7) = log(exp(25/16)-1);


% redefine the length number of parameters in this case
n_G_prm = length(muPhi);

% define sigma for each param
sigmaPhi = zeros(n_G_prm);
sigmaPhi(1,1) = 1;
sigmaPhi(2,2) = 1;
sigmaPhi(3,3) = 1;
sigmaPhi(4,4) = 1;
sigmaPhi(5,5) = 1;
sigmaPhi(6,6) = 1;
sigmaPhi(7,7) = 1;

%for i_simu=1:6
    
    sigma_from_data = [3.2,2.2,2.5,3.6,2.5,1.84,5.25];
    % create normally distributed gaussian values supposed to be our participants real sensitivities.
    for i_param = 1:n_G_prm
        % version with sd taken from distribution of our dataset
        %     parameterPhi(i_param,:) = randn(nb_run,1).*(sigma_from_data(i_param)) + muPhi(i_param);
        
        % version with uninformed prior
        parameterPhi(i_param,:) = randn(nb_run,1).*(sigmaPhi(i_param,i_param)*3) + muPhi(i_param);
    end
    
    for i_iteration = 1:nb_iteration
        
        prefix_name = ['iteration_nb_',num2str(i_iteration),'_nb_param_',num2str(number_of_parameters),'_prior_',prior_type,'_model_nb_',num2str(model_nb)];
        
        %% compute choice_matrices, prepare the variables
        
        load('plan_B_bestMatrix_bis.mat')
        choice_opt = bestMatrix;
        
        n_sessions = 4;
        n_trials = n_trialsPerSession*n_sessions;
        
        % for each 'participant''
        for i_run = 1:nb_run
            i_run
            % sequence name
            sequence_name = [prefix_name,'_run_',num2str(i_run)];
            
            
            %     DeltaI, avg over participant deltaE         RP_idx            Ep_or_Em          trial_idx               E_chosen            SUm_Echosen
            var = [nanmean(var_saved(1,:,:),3);var_saved(2,:,1);var_saved(3,:,1);var_saved(4,:,1);var_saved(5,:,1);nanmean(var_saved(6,:,:),3);nanmean(var_saved(7,:,:),3);...
                %       SumAUC+SumNb_answer avg     SumAUC+ratio avg
                mean(var_saved(8,:,:),3);nanmean(var_saved(9,:,:),3)];
            
            options.sources.type = 0; % binary data
            options.DisplayWin  = 0; % display figure during inversion
            options.verbose     = 0; % display text during inversion (1) or not (0)
            % define number of parameters to estimate
            options.dim = struct('n', n_hiddenStates,... % number of hidden states
                'n_t', n_trials,... % number of trials across runs
                'n_theta',n_F_prm,... % number of evolution parameters
                'n_phi',n_G_prm); % number of observation parameters
            
            %% Launch VBA simulations
            
            for block_i = 1:4
                sum_AUC =0;
                m_ratio = 0;
                for trial_i = 1:54
                    
                    if  var(4,trial_i+54*(block_i-1)) == 1
                        var_i = [var(1:8,trial_i+54*(block_i-1));sum_AUC];
                    else
                        var_i = [var(1:8,trial_i+54*(block_i-1));m_ratio];
                    end
                    save_var_i(trial_i+54*(block_i-1),i_run) = var_i(9);
                    options.dim.n_t = 1;
                    [choiceLeftOptionSimu(trial_i+54*(block_i-1),i_run), ~, ~, ~, ~, ~] = VBA_simulate(1, f_fname, g_fname,...
                        theta, parameterPhi(:,i_run), var_i, alpha, sigma, options, x0);
                    p(trial_i+54*(block_i-1),i_run) = choiceLeftOptionSimu(trial_i+54*(block_i-1),i_run);
                    % create a thousand values, using the reported proba for the number of 1 in the vector
                    proba_vector = zeros(10000,1);
                    proba_vector(1:round(choiceLeftOptionSimu(trial_i+54*(block_i-1),i_run)*10000),1) = 1;
                    % shuffle it, and use the first value as the choice picked by a participant
                    random_idx = randperm(10000);
                    proba_vector_random = proba_vector(random_idx);
                    
                    if  var(4,trial_i+54*(block_i-1)) == 1
                        if proba_vector_random(1) == 0
                            sum_AUC = sum_AUC+44.5/1000;
                        elseif proba_vector_random(1) == 1 && var(2,trial_i+54*(block_i-1)) == 1
                            sum_AUC = sum_AUC+100/1000;
                        elseif proba_vector_random(1) == 1 && var(2,trial_i+54*(block_i-1)) == 2
                            sum_AUC = sum_AUC+186/1000;
                        elseif proba_vector_random(1) == 1 && var(2,trial_i+54*(block_i-1)) == 3
                            sum_AUC = sum_AUC+263/1000;
                        end
                    else
                        if trial_i <= 9
                            if proba_vector_random(1) == 0
                                m_ratio =   ratio(1,1);
                            elseif proba_vector_random(1) == 1 && var(2,trial_i+54*(block_i-1)) == 1
                                m_ratio =   ratio(1,2);
                            elseif proba_vector_random(1) == 1 && var(2,trial_i+54*(block_i-1)) == 2
                                m_ratio =   ratio(1,3);
                            elseif proba_vector_random(1) == 1 && var(2,trial_i+54*(block_i-1)) == 3
                                m_ratio =   ratio(1,4);
                            end
                        elseif trial_i > 9 && trial_i <= 18
                            if proba_vector_random(1) == 0
                                m_ratio =   ratio(2,1);
                            elseif proba_vector_random(1) == 1 && var(2,trial_i+54*(block_i-1)) == 1
                                m_ratio =   ratio(2,2);
                            elseif proba_vector_random(1) == 1 && var(2,trial_i+54*(block_i-1)) == 2
                                m_ratio =   ratio(2,3);
                            elseif proba_vector_random(1) == 1 && var(2,trial_i+54*(block_i-1)) == 3
                                m_ratio =   ratio(2,4);
                            end
                        elseif trial_i > 18 && trial_i <= 27
                            if proba_vector_random(1) == 0
                                m_ratio =   ratio(3,1);
                            elseif proba_vector_random(1) == 1 && var(2,trial_i+54*(block_i-1)) == 1
                                m_ratio =   ratio(3,2);
                            elseif proba_vector_random(1) == 1 && var(2,trial_i+54*(block_i-1)) == 2
                                m_ratio =   ratio(3,3);
                            elseif proba_vector_random(1) == 1 && var(2,trial_i+54*(block_i-1)) == 3
                                m_ratio =   ratio(3,4);
                            end
                        elseif trial_i > 27 && trial_i <= 36
                            if proba_vector_random(1) == 0
                                m_ratio =   ratio(4,1);
                            elseif proba_vector_random(1) == 1 && var(2,trial_i+54*(block_i-1)) == 1
                                m_ratio =   ratio(4,2);
                            elseif proba_vector_random(1) == 1 && var(2,trial_i+54*(block_i-1)) == 2
                                m_ratio =   ratio(4,3);
                            elseif proba_vector_random(1) == 1 && var(2,trial_i+54*(block_i-1)) == 3
                                m_ratio =   ratio(4,4);
                            end
                        elseif trial_i > 36 && trial_i <= 45
                            if proba_vector_random(1) == 0
                                m_ratio =   ratio(5,1);
                            elseif proba_vector_random(1) == 1 && var(2,trial_i+54*(block_i-1)) == 1
                                m_ratio =   ratio(5,2);
                            elseif proba_vector_random(1) == 1 && var(2,trial_i+54*(block_i-1)) == 2
                                m_ratio =   ratio(5,3);
                            elseif proba_vector_random(1) == 1 && var(2,trial_i+54*(block_i-1)) == 3
                                m_ratio =   ratio(5,4);
                            end
                        elseif trial_i > 45 && trial_i <= 54
                            if proba_vector_random(1) == 0
                                m_ratio =   ratio(6,1);
                            elseif proba_vector_random(1) == 1 && var(2,trial_i+54*(block_i-1)) == 1
                                m_ratio =   ratio(6,2);
                            elseif proba_vector_random(1) == 1 && var(2,trial_i+54*(block_i-1)) == 2
                                m_ratio =   ratio(6,3);
                            elseif proba_vector_random(1) == 1 && var(2,trial_i+54*(block_i-1)) == 3
                                m_ratio =   ratio(6,4);
                            end
                        end
                    end
                    %                 if trial_i ==27 && block_i ==2
                    %                 'breakpoint'
                    %                 [gx] = g_observation55( [], parameterPhi(:,i_run), var_i, [] )
                    %                 end
                    keep_ratio(trial_i+54*(block_i-1),i_run) = m_ratio;
                    keep_AUC(trial_i+54*(block_i-1),i_run) = sum_AUC;
                end
            end
            
            %         0, 0.25 0.75 1
            for n = 1:size(choiceLeftOptionSimu,1)
                if 0 <= choiceLeftOptionSimu(n,i_run) && 0.25 > choiceLeftOptionSimu(n,i_run)
                    choiceLeftOptionSimu(n,i_run) = 0;
                elseif 0.25 <= choiceLeftOptionSimu(n,i_run) && 0.5 > choiceLeftOptionSimu(n,i_run)
                    choiceLeftOptionSimu(n,i_run) = 0.25;
                elseif 0.5 <= choiceLeftOptionSimu(n,i_run) && 0.75 >= choiceLeftOptionSimu(n,i_run)
                    choiceLeftOptionSimu(n,i_run) = 0.75;
                elseif 0.75 < choiceLeftOptionSimu(n,i_run) && 1 >= choiceLeftOptionSimu(n,i_run)
                    choiceLeftOptionSimu(n,i_run) = 1;
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
            [posterior,out] = VBA_NLStateSpaceModel(choiceLeftOptionSimu(:,i_run)', var, f_fname, g_fname, options.dim, options);
            
            %% save the data
            run_nm = ['run_nb_',num2str(i_run)];
            toBeSaved.(run_nm).n_trials = n_trials;
            toBeSaved.(run_nm).theta = theta;
            toBeSaved.(run_nm).phi = parameterPhi;
            toBeSaved.(run_nm).u = var;
            toBeSaved.(run_nm).simulated_data_left_opt = choiceLeftOptionSimu;
            toBeSaved.(run_nm).posterior = posterior;
            toBeSaved.(run_nm).out = out;
            toBeSaved.(run_nm).options = options;
            
            
            close all
            
            posterior_session_run(i_run,:) = posterior;
        end
        %         save([results_folder,sequence_name,'simulation_data.mat'],'toBeSaved');
        
        figure()
        plot(nanmean(p,2))
        % plot to check the decrease over time in physical trials, and increase in mental trials.
        figure()
        plot([nanmean(nanmean(p(1:10,:))),nanmean(nanmean(p(10:20,:))),nanmean(nanmean(p(20:30,:))),nanmean(nanmean(p(30:40,:))),nanmean(nanmean(p(40:54,:))),...
            nanmean(nanmean(p(55:65,:))),nanmean(nanmean(p(65:75,:))),nanmean(nanmean(p(75:85,:))),nanmean(nanmean(p(85:95,:))),nanmean(nanmean(p(95:108,:)))])
        
        
        % extract posterior for each participant
        for j = 1:nb_run
            data_tmp = posterior_session_run(j).muPhi';
            param_per_run(j,:) = data_tmp ;
        end
        param_simu.fake_parameterPhi = parameterPhi';
        param_simu.true_parameterPhi = param_per_run;
%         simu_name = strcat('non_specific_param_simu',num2str(i_simu),'.mat');
%         save(simu_name,'param_simu')
    end
% end

%% prepare all the results from the simulations. model
if use_30ksim == true
    parameterPhi = [];
    param_per_run = [];
    for i = 1:6
        load(strcat('non_specific_param_simu',num2str(i),'.mat'))
        parameterPhi = [parameterPhi;param_simu.fake_parameterPhi];
        param_per_run = [param_per_run;param_simu.true_parameterPhi];
    end
end

%% plot results
if show_plot == true
    
    
    %         figure()
    %         x_param_idx = 1:n_G_prm;
    %
    %
    %         % compute SEM and SD
    %         err_sem = squeeze(std(param_per_run(:,:))./sqrt(size(param_per_run(:,:),1)));
    %         err = squeeze(std(param_per_run(:,:)));
    %         %     if wanted, compute correlation between ''real'' and posterior data
    %         for param_i = 1:n_G_prm
    %             corr_tmp = corrcoef(param_per_run(:,param_i),parameterPhi(param_i,:)');
    %             corr_mat(1,param_i) =corr_tmp(1,2);
    %         end
    %
    %         % plot the average over N simulation.
    %         bar(x_param_idx,mean(squeeze(param_per_run(:,:))))
    %         mean(squeeze(param_per_run(:,:)));
    %         hold on
    %         er = errorbar(x_param_idx,mean(squeeze(param_per_run(:,:))),err,err);
    %         er.Color = [0 0 0];
    %         er.LineStyle = 'none';
    %         er = errorbar(x_param_idx,mean(squeeze(param_per_run(:,:))),err_sem,err_sem);
    %         er.Color = [1 0 0];
    %         er.LineStyle = 'none';
    %         set(gca,'XTick',linspace(1,n_G_prm,n_G_prm))
    %         set(gca,'XTickLabel',{'kR','kP','kEp','kEm','biais','kFp','kLm'})
    %         title(['Mean ',num2str(i_run),' simulations with ',num2str(session),' session of efforts'])
    %         ylim
    %         xlim=get(gca,'xlim');
    %         hold on
    %         %         plot(xlim,[0.3 0.3],'m')
    %         %         plot(xlim,[0.6 0.6],'c')
    %         %         ylim([-1 1])
    figure()
    subplot(1,2,1)
    % plot correlations between ''real'' and posterior for each nb of sessions
    
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
    title(['Model recovery'])
    hcb = colorbar;
    set(get(hcb,'label'),'string','Correlation coefficient (R)','Rotation',90.0,'FontSize',16);
    colormap(redblue(45))
    
    % autocorrelation plot
    
    
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
    title(['Confusion matrix'])
    hcb = colorbar;
    set(get(hcb,'label'),'string','Correlation coefficient (R)','Rotation',90.0,'FontSize',16);
    colormap(redblue(45))
    mean([corrmat(1,1),corrmat(2,2),corrmat(3,3),corrmat(4,4),corrmat(5,5),corrmat(6,6),corrmat(7,7)])
    
end