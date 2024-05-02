function[CID_nb, mod, y_hat, NV, AIC, BIC, R2, free_energy_matrix, RMSE, MAE] = ModelSelectionAndFeatureExtraction_NC_modified()
% script to compute behavioral sensitivities and select the best one
%
% OUTPUTS
% CID_nb: list of subjects
%
% mod: structure with behavioral parameters
%
% y_hat: matrix with p(high E) for all subjects all trials all runs
%
% NV: matrix with net value for all subjects all trials all runs
%
% AIC, BIC, R2, free_energy_matrix: scores for each model
%
% RMSE, MAE: error for each model

%% clean workspace before starting

clc;
clearvars;
close all;
set(0,'defaultfigurecolor',[1 1 1])
%% prepare paths and working directory
% serverRoot = 'D:';
serverRoot = '\\sv-nas1.rcp.epfl.ch\Sandi-Lab\human_data_private\Summary\';
% main_folder                 = [pwd]; % you have to be sure that you are in the correct path when you launch the script
results_folder           = '\\sv-nas1.rcp.epfl.ch\Sandi-Lab\human_data_private\raw_data_subject\study1';
code_folder              = [serverRoot,'\matlab_codes\Data_simulation\'];

addpath(results_folder);
addpath(genpath(code_folder));
cd(results_folder)
addpath([serverRoot,'\matlab_codes\Main_analysis_all_dataset_codes']);

%% set run parameters

% first test phase
% f_fname = {[],[],[],[],[],[],[],@f_observation8,@f_observation8,@f_observation8,@f_observation8,@f_observation8,@f_observation8,@f_observation8,[],[],...
%     @f_observation8,@f_observation8,@f_observation8,@f_observation8,@f_observation8,@f_observation8,@f_observation8,@f_observation8};
% is_exponential = [false,false,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true,true];
% is_multisession = true;
% g_fname = {@g_observation1,@g_observation2,@g_observation3,@g_observation4,@g_observation5,@g_observation6,@g_observation7,...
%     @g_observation8,@g_observation9,@g_observation10,@g_observation11,@g_observation11,@g_observation12,@g_observation13,...
%     @g_observation14,@g_observation15,@g_observation16,@g_observation17,@g_observation18,@g_observation19,@g_observation20,@g_observation21,@g_observation22,@g_observation23};
% % WARNING TO CHECK ! THERE WERE NOT THE SAME NB OF PARAM (21-22 different type, but g_fname had 23
% % of them, so we need to check it.
% n_G_prm = [4,6,6,7,7,7,7,4,5,3,3,3,2,4,3,3,4,4,4,4,5,6,6,4];
% n_F_prm = [0,0,0,0,0,0,0,2,2,2,2,2,2,2,0,0,2,2,2,2,2,2,2,2];
% n_hiddenStates = [0,0,0,0,0,0,0,2,2,2,2,2,2,2,0,0,2,2,2,2,2,2,2,2];
%% result from first phase, basic model with adaptive fatigue best
% second test phase
% % groupe f functions to test
% f_fname = {@f_observation8,@f_observation9,@f_observation10,@f_observation8,@f_observation8...
%           ,@f_observation9,@f_observation10,@f_observation8,@f_observation8,@f_observation8...
%           ,@f_observation8,@f_observation8,@f_observation8,@f_observation8,@f_observation8,@f_observation8};
%       
% % apply multisession to each block. 
% is_multisession = true;
% 
% % groupe g functions to test
% g_fname = {@g_observation8,@g_observation8,@g_observation8,@g_observation24,@g_observation31...
%           ,@g_observation24,@g_observation24,@g_observation25,@g_observation26,@g_observation27...
%           ,@g_observation28,@g_observation29,@g_observation30,@g_observation24,@g_observation24,@g_observation24};
% 
% % Set the number of G parameters. For F param, they are the same number as hidden states in our case
% n_G_prm = [4,4,4,5,6,5,5,5,5,5,3,3,5,5,5,5];
% n_F_prm = [2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2];
% n_hiddenStates = [2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2];
%% results: model with bias on top of it, is best
% % third phase
% % groupe f functions to test
% f_fname = {@f_observation8,@f_observation8,[],[]};
%       
% % apply multisession to each block. 
% is_multisession = true;
% 
% % groupe g functions to test
% g_fname = {@g_observation8,@g_observation24,@g_observation32,@g_observation33};
% 
% % Set the number of G parameters. For F param, they are the same number as hidden states in our case
% n_G_prm = [4,5,5,7];
% n_F_prm = [2,2,2,0];
% n_hiddenStates = [2,2,2,0];
% n_trialsPerSession = 54;

%% results. Have all the parameters in G_obs change nothing. Now we have to define priors and the fct forcing positivity
% % second third phase
% % groupe f functions to test
% f_fname = {[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]};
%       
% % apply multisession to each block. 
% is_multisession = true;
% 
% % groupe g functions to test
% g_fname = {@g_observation37,@g_observation40,@g_observation41,@g_observation42,@g_observation43,@g_observation44,@g_observation45,...
%     @g_observation46,@g_observation47,@g_observation48,@g_observation49,@g_observation50,@g_observation51,@g_observation52,...
%     @g_observation53,@g_observation54};
% 
% % Set the number of G parameters. For F param, they are the same number as hidden states in our case
% n_G_prm = [7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7];
% n_F_prm = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
% n_hiddenStates = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
% n_trialsPerSession = 54;
%% 4th run, optimizing model
% 4 phase
% groupe f functions to test
f_fname = {[],[],[],[],[]};
%       
% % apply multisession to each block. 
is_multisession = true;
% 
% % groupe g functions to test current best @g_observation55
% test model : with kI (g_observation66), with kI+kE (g_observation1), with kI+kE+bias (g_observation64),with kI+kE+bias+kFlinear (g_observation65), with kI+kE+bias+kFadapted (g_observation55) 
g_fname = {@g_observation66,@g_observation1,@g_observation64,@g_observation65,@g_observation55};

% % Set the number of G parameters. For F param, they are the same number as hidden states in our case
n_G_prm = [2,4,5,7,7];
n_F_prm = [0,0,0,0,0];
n_hiddenStates = [0,0,0,0,0];
n_trialsPerSession = 54;

%% test run, trying out other models
% 
% f_fname = {[],[],[]};
% %       
% % % apply multisession to each block. 
% is_multisession = true;
% % 
% % % groupe g functions to test current best @g_observation55
% % test model : with kI, with kI +kE, with kI +kE+bias,with kI,kE+bias+kFlinear, with kI,kE+bias+kFadapted, 
% g_fname = {@g_observation64,@g_observation68,@g_observation55};
% 
% % % Set the number of G parameters. For F param, they are the same number as hidden states in our case
% n_G_prm = [5,6,7];
% n_F_prm = [0,0,0];
% n_hiddenStates = [0,0,0];
% n_trialsPerSession = 54;
 %% Just best model
% % % 4 phase
% % % groupe f functions to test
% f_fname = {[]};
% %       
% % % apply multisession to each block. 
% is_multisession = true;
% % 
% % % groupe g functions to test current best @g_observation55
% % test model : with kI, with kI +kE, with kI +kE+bias,with kI,kE+bias+kFlinear, with kI,kE+bias+kFadapted, 
% g_fname = {@g_observation55};
% 
% % % Set the number of G parameters. For F param, they are the same number as hidden states in our case
% n_G_prm = [7];
% n_F_prm = [0];
% n_hiddenStates = [0];
% n_trialsPerSession = 54;

nb_functions = length(g_fname);
% change the input date [-2,-1,1,2] to [0,0.25,0.75,1]
rescale_choice = true;

% do you want binary response from the model
binary_answers = false; % make response binary from model prediction
do_smooth = false;
do_save = false; % save the corresponding results
do_plot = true;
do_kEc_plot = true;

% include current efficiency or previous trial efficiency
include_currEff = false;
include_prevEff = true;

%% extract behavioral data and put it in a struct
% the list of saturators is updated by hand, as the end of the script gives the results
% remove saturators, removes any participants which modeling saturates on 1 run
remove_choice_saturators = true;
remove_choice_saturators_one_run = false;
% remove saturators, removes any participants which choices saturated on 1 run
remove_pred_saturators = false;
% Extract all information for modelling
[p_or_m,nb_participants,all_data,all_files,CID_nb] = extract_behavioral_data(remove_choice_saturators,remove_choice_saturators_one_run,remove_pred_saturators);

%% for each subject and each model, compute modelisation and compute error
RT_mental=NaN(54,72,2);
RT_phyiscal=NaN(54,72,2);
for i_sub = 1:length(all_files)
    i_sub % plot at what subject we currently are.
    [var,all_choices,all_choices_rescaled,deltaRP_idx,deltaRP,deltaE,tmp1,tmp2,all_MVC(i_sub,:),all_NMP(i_sub,:),physical_IP(i_sub)...
        ,mental_IP(i_sub),perf(i_sub),incentive_idx] = reshape_data_NC_modified(i_sub,all_data,rescale_choice,binary_answers,p_or_m, include_currEff, include_prevEff);
    % use the version rescaled as all_choices
    %     [var,~,all_choices,deltaRP_idx,deltaRP,deltaE,RT_physical,RT_mental,all_MVC(i_sub,:),all_NMP(i_sub,:),physical_IP(i_sub)...
    %         ,mental_IP(i_sub),perf(i_sub),incentive_idx] = reshape_data(i_sub,all_data,rescale_choice,binary_answers,p_or_m);
    if size(tmp1,1) == 2
        RT_physical(:,i_sub,:) = tmp1';
    elseif size(tmp1,1) == 1
        RT_physical(:,i_sub,1) = tmp1';
    end
    if size(tmp2,1) == 2
        RT_mental(:,i_sub,:) = tmp2';
    elseif size(tmp2,1) == 1
        RT_mental(:,i_sub,1) = tmp2;
    end
    
    

    % compute decision making reaction time and save it
    RTP(i_sub) = mean(mean(RT_physical(:,i_sub,:)));
    RTM(i_sub) = mean(mean(RT_mental(:,i_sub,:)));

    % save participants choices.
    sub_i_choices(:,i_sub) = all_choices;
    var_p_or_m(:,i_sub) = var(4,:);
    var_saved(:,:,i_sub) = var;
	% choose how many models to test
    for i_model=1:nb_functions
        
        if i_model == 6
                if var_p_or_m(1,i_sub) == 1
                    all_choices(1:54) = NaN;
                    all_choices(109:162) = NaN;
                else
                    all_choices(55:108) = NaN;
                    all_choices(163:216) = NaN;
                end
        end
        [options] = select_options(i_model,all_choices,n_hiddenStates,n_trialsPerSession,n_F_prm,n_G_prm,i_sub,is_multisession);
        
        % since NaN makes the toolbox crash, and that we ignore NaN choices using isYout in the
        % options, we redefine the no choice with a -1 idx, to take it into account later
        all_choices(isnan(all_choices)) = -1;
        % for i_sub 26, there is suddenly a nan popping out, take care ouf it for the function.
        var(9,isnan(var(9,:))) = 0;
        var(10,isnan(var(10,:))) = 0;
        var(13,isnan(var(10,:))) = 0;
        if i_sub == 26
            for i = 1:size(var,1)
                var(i,isnan(var(i,:))) = 0;
            end
        end
        % adapt script according to if you want to predict binary choices
        % (0=low E, 1=high E) or choices incorporating confidence
        % information (0=low E, high Conf; 0.25=low E, low conf; 0.75=high
        % E, low Conf; 1=high E, high Conf)
        switch binary_answers
            case false
                [posterior,out] = VBA_NLStateSpaceModel(all_choices, var, f_fname{i_model}, g_fname{i_model}, options.dim, options);
            case true
                [posterior,out] = VBA_NLStateSpaceModel(all_choices_rescaled, var, f_fname{i_model}, g_fname{i_model}, options.dim, options);
        end

        %put back NaN for easier processing afterwards
        all_choices(all_choices == -1) = NaN;
        
        
        % extract parameters
        sensitivitiesPhi{i_model,i_sub} = posterior.muPhi;
        sensitivitiesTheta{i_model,i_sub} = posterior.muTheta;
        
        % extract free energy matrix and AIC,BIC for model comparison
        free_energy_matrix(i_model, i_sub) = out.F;
        AIC(i_model, i_sub) = out.fit.AIC;
        BIC(i_model, i_sub) = out.fit.BIC;
        R2(i_model, i_sub) = out.fit.R2;
        
        % Extract model predictions for each choice.
        y_hat(:,i_model,i_sub) = out.suffStat.gx;
        switch i_model
            case {1,2} % no bias
                NV(:,i_model,i_sub) = -log((1-out.suffStat.gx)./out.suffStat.gx);
            case {3,4,5} % bias to remove
                NV(:,i_model,i_sub) = -posterior.muPhi(5) - log((1-out.suffStat.gx)./out.suffStat.gx);
        end
        
        %% transform p(choice = high E) extracted from VBA into a binary or a four-level variable like the input VBA received
        % (but in general, what you want should be the probability
        % extracted from the model)
        for i_choice = 1:n_trialsPerSession*4
            %             % extract predictions as only 2 options.
            %             if binary_answers == true
            %                 if y_hat(i_choice,i_model,i_sub) >=0 && y_hat(i_choice,i_model,i_sub) < 0.25
            %                     y_hat(i_choice,i_model,i_sub) = 0;
            %                 elseif y_hat(i_choice,i_model,i_sub) >=0.25 && y_hat(i_choice,i_model,i_sub) < 0.5
            %                     y_hat(i_choice,i_model,i_sub) = 0;
            %                 elseif y_hat(i_choice,i_model,i_sub) >=0.5 && y_hat(i_choice,i_model,i_sub) <= 0.75
            %                     y_hat(i_choice,i_model,i_sub) = 1;
            %                 elseif y_hat(i_choice,i_model,i_sub) >0.75 && y_hat(i_choice,i_model,i_sub) <= 1
            %                     y_hat(i_choice,i_model,i_sub) = 1;
            %                 end
            %
            %             else
            if do_smooth == true
                % smooth prediction choices in 4 category
                if y_hat(i_choice,i_model,i_sub) >=0 && y_hat(i_choice,i_model,i_sub) < 0.25
                    y_hat(i_choice,i_model,i_sub) = 0;
                elseif y_hat(i_choice,i_model,i_sub) >=0.25 && y_hat(i_choice,i_model,i_sub) < 0.5
                    y_hat(i_choice,i_model,i_sub) = 0.25;
                elseif y_hat(i_choice,i_model,i_sub) >=0.5 && y_hat(i_choice,i_model,i_sub) < 0.75
                    y_hat(i_choice,i_model,i_sub) = 0.75;
                elseif y_hat(i_choice,i_model,i_sub) >=0.75 && y_hat(i_choice,i_model,i_sub) <= 1
                    y_hat(i_choice,i_model,i_sub) = 1;
                end
            end
            %             end
        end
        
        % compute RMSE or MAE depending on the treatment applied to the data.
        if binary_answers == true
            % warning, make sure all_choices is the rescaled version
            RMSE(i_model,i_sub) = sqrt(nanmean(((all_choices-y_hat(:,i_model,i_sub)').^2)));
            MAE(i_model,i_sub) = nanmean(abs(all_choices-y_hat(:,i_model,i_sub)'));
        else
            
            RMSE(i_model,i_sub) = sqrt(nanmean(((all_choices(~isnan(all_choices))-y_hat(~isnan(all_choices),i_model,i_sub)').^2)));
            MAE(i_model,i_sub) = nanmean(abs(all_choices(~isnan(all_choices))-y_hat(~isnan(all_choices),i_model,i_sub)'));
            
            % Compute physical and mental MAE, perhaps it underfit one system.
            switch p_or_m(i_sub)
                % compute both physical and mental MAE, see where is has a hard time fitting. divide by 12
                % because normally its 216, and it's already grouped (divided by 2) and binned (divided by 9) = 12
                case "p"
                    MAE_physical(i_model,i_sub) = nansum(abs(all_choices(1:54)-y_hat(1:54,i_model,i_sub)')+abs(all_choices(109:162)-y_hat(109:162,i_model,i_sub)'))/108;
                    MAE_mental(i_model,i_sub) = nansum(abs(all_choices(55:108)-y_hat(55:108,i_model,i_sub)')+abs(all_choices(163:216)-y_hat(163:216,i_model,i_sub)'))/108;
                case "m"
                    MAE_mental(i_model,i_sub) = nansum(abs(all_choices(1:54)-y_hat(1:54,i_model,i_sub)')+abs(all_choices(109:162)-y_hat(109:162,i_model,i_sub)'))/108;
                    MAE_physical(i_model,i_sub) = nansum(abs(all_choices(55:108)-y_hat(55:108,i_model,i_sub)')+abs(all_choices(163:216)-y_hat(163:216,i_model,i_sub)'))/108;
            end
        end
        
        % extract all runs, and run 1 and 2 averaged together
        [nb_ND_physical(i_sub),nb_ND_mental(i_sub),bins_per_trial_choice_physical(:,i_sub),bins_per_trial_choice_mental(:,i_sub),...
            bins_per_trial_pred_physical(:,i_model,i_sub),bins_per_trial_pred_mental(:,i_model,i_sub),grouped_choice_physical_run(:,i_sub),...
            grouped_choice_mental_run(:,i_sub),grouped_pred_physical_run(:,i_model,i_sub),grouped_pred_mental_run(:,i_model,i_sub),...
            saturates_choice_once_per_sub,saturates_choice_twice_per_sub,saturates_pred_per_sub] = extract_individual_sessions(all_choices,sub_i_choices(:,i_sub),y_hat(:,i_model,i_sub),p_or_m(i_sub));
        
        %save if participant saturated
        idx_CID_saturates_choice_once(i_sub) = saturates_choice_once_per_sub;
        idx_CID_saturates_choice_twice(i_sub) = saturates_choice_twice_per_sub;
        idx_CID_saturates_pred(i_model,i_sub) = saturates_pred_per_sub;
        
        % as parameters are not reset in the struct option, clear it and put new data within
        clear options
        % compute participants choices and prediction from the model in heatmaps
        prediction_or_choice = true;
        [phys_choice_mat,ment_choice_mat,nb_trials_choice_mat_physical,nb_trials_choice_mat_mental] = create_choice_mat(y_hat(:,i_model,i_sub)',p_or_m(i_sub),deltaRP_idx,deltaRP,deltaE,prediction_or_choice);
        pred_mat_physical(:,:,i_sub,i_model) = phys_choice_mat;
        pred_mat_mental(:,:,i_sub,i_model) = ment_choice_mat;
        
    end
    

    
    %Extract choice matrix
    prediction_or_choice = false;
    [tmp_phys_choice_mat,tmp_ment_choice_mat,nb_trials_choice_mat_physical,nb_trials_choice_mat_mental] = create_choice_mat(all_choices,p_or_m(i_sub),deltaRP_idx,deltaRP,deltaE,prediction_or_choice);
    choice_mat_physical(:,:,i_sub) = tmp_phys_choice_mat;
    choice_mat_mental(:,:,i_sub) = tmp_ment_choice_mat;
    nb_trials_mat_physical(:,:,i_sub) = nb_trials_choice_mat_physical;
    nb_trials_mat_mental(:,:,i_sub) = nb_trials_choice_mat_mental;
    
    

end
% find CID to potentially exclude
CID_saturate_choice_once = [];
CID_saturate_choice_twice = [];
CID_saturate_pred = [];
for i = 1:length(RMSE)
    %find saturators
    if idx_CID_saturates_choice_once(i) == 1
        CID_saturate_choice_once(end+1) = CID_nb(i);
    end
        if idx_CID_saturates_choice_twice(i) == 1
        CID_saturate_choice_twice(end+1) = CID_nb(i);
    end
    for i_model = 1:size(idx_CID_saturates_pred,1)
        if idx_CID_saturates_pred(i_model,i) == 1
            CID_saturate_pred(i_model,i) = CID_nb(i);
        end
    end

end
% list of participant saturating only 1 run, not 1 or 2.
CID_saturate_choice_once = setdiff(CID_saturate_choice_once,CID_saturate_choice_twice);
unique(CID_saturate_pred);

% Compute confidence for each participant
conf = NaN(216,length(all_files));
for i_sub = 1:length(all_files)
    tmp = sub_i_choices(:,i_sub);
    low_conf = find(tmp == 0.25|tmp == 0.75);
    high_conf = find(tmp == 0|tmp == 1);
    find_NA = find(isnan(tmp));
    conf(low_conf,i_sub) = 0; conf(high_conf,i_sub) = 1; conf(find_NA,i_sub) = NaN;
end
%% compare models
model_i =5;
% [posterior_BMC, out_BMC] = VBA_groupBMC([free_energy_matrix(:,[1:25,27:69])]);
[posterior_BMC, out_BMC] = VBA_groupBMC([free_energy_matrix(:,[1:25,27:69])]);

%% prepare perf for plots

[m,p,tmp] = compute_perf_score(perf,incentive_idx,var_saved(2,1:54,:),var_saved(3,1:54,:),sub_i_choices,var_p_or_m,model_i,RT_mental,RT_physical);
%% plot results
if do_plot == true
% plot the results between models
%plot_model_selection(out_BMC,AIC,MAE,MAE_physical,MAE_mental,BIC,R2);

% plot the results between models
plot_model_selection(out_BMC,AIC,MAE,[],[],BIC,R2);

plot_figure_paper2(m,p,tmp,RT_mental)

% result for the specific model
plot_choice_and_prediction_matrix(choice_mat_physical,choice_mat_mental,pred_mat_physical,pred_mat_mental,model_i);


plot_choice_prediction_effort_incentive(choice_mat_physical,choice_mat_mental,pred_mat_physical,pred_mat_mental,nb_trials_mat_physical,nb_trials_mat_mental,model_i,CID_nb);


plot_per_trial_and_SV_choice_prediction(grouped_choice_physical_run,grouped_choice_mental_run,grouped_pred_physical_run(:,model_i,:),...
    grouped_pred_mental_run(:,model_i,:),bins_per_trial_choice_physical,bins_per_trial_choice_mental,bins_per_trial_pred_physical(:,model_i,:),...
    bins_per_trial_pred_mental(:,model_i,:),all_NMP(:,1))

plot_choice_conf_over_time(sub_i_choices,conf,p_or_m,perf)
end
%% Prepare perf dataset

if do_kEc_plot == true
for i_sub = 1:size(sub_i_choices,2)
    if var_p_or_m(1,i_sub) == 1
        m_HE(:,1,i_sub) = sub_i_choices(55:108,i_sub);
        m_HE(:,2,i_sub) = sub_i_choices(163:216,i_sub);
    else
        m_HE(:,1,i_sub) = sub_i_choices(1:54,i_sub);
        m_HE(:,2,i_sub) = sub_i_choices(109:162,i_sub);
    end
end

k_sub = 1;
j_sub = 1;
for i_sub = 1:length(perf)
    kEc(i_sub) = sensitivitiesPhi{i_model,i_sub}(7);
    if kEc(i_sub) <= 1.7852
        low_kEc_efficiency(:,:,j_sub) = perf(i_sub).m_ratio;
        low_kEc_energization(:,:,j_sub) = perf(i_sub).m_ratio*kEc(i_sub);
        low_kEc_performance(:,:,j_sub) = perf(i_sub).m_performance;
        low_kEc_HE(:,:,j_sub) = m_HE(:,:,i_sub);
        low_kEc_std_HE(:,j_sub) = nanstd(m_HE(:,:,i_sub),0,1);
        j_sub = j_sub+1;
    else
        high_kEc_efficiency(:,:,k_sub) = perf(i_sub).m_ratio;
        high_kEc_energization(:,:,k_sub) = perf(i_sub).m_ratio*kEc(i_sub);
        high_kEc_performance(:,:,k_sub) = perf(i_sub).m_performance;
        high_kEc_HE(:,:,k_sub) = m_HE(:,:,i_sub);
        high_kEc_std_HE(:,k_sub) = nanstd(m_HE(:,:,i_sub),0,1);
        k_sub = k_sub+1;
    end
end

% compute efficiency t-1 and t corr.
for i_run = 1:2
    for i_sub = 1:size(low_kEc_efficiency,3)
        low_kEc_eff_tp1(:,i_run,i_sub) = cat(1,0,low_kEc_efficiency(1:53,i_run,i_sub));
            tmp = corrcoef(low_kEc_efficiency(:,i_run,i_sub),low_kEc_eff_tp1(:,i_run,i_sub),'rows','pairwise');
        corr_lkEc_eff_efftp1(i_run,i_sub) = tmp(1,2);
                    tmp = corrcoef(low_kEc_efficiency(:,i_run,i_sub)*kEc(i_sub),low_kEc_eff_tp1(:,i_run,i_sub),'rows','pairwise');
        corr_lkEc_en_efftp1(i_run,i_sub) = tmp(1,2);
    end
        for i_sub = 1:size(high_kEc_efficiency,3)
        high_kEc_eff_tp1(:,i_run,i_sub) = cat(1,0,high_kEc_efficiency(1:53,i_run,i_sub));
            tmp = corrcoef(high_kEc_efficiency(:,i_run,i_sub),high_kEc_eff_tp1(:,i_run,i_sub),'rows','pairwise');
        corr_hkEc_eff_efftp1(i_run,i_sub) = tmp(1,2);
                    tmp = corrcoef(high_kEc_efficiency(:,i_run,i_sub)*kEc(i_sub),high_kEc_eff_tp1(:,i_run,i_sub),'rows','pairwise');
        corr_hkEc_en_efftp1(i_run,i_sub) = tmp(1,2);
    end
end
corr_lkEc_eff_efftp1 = nanmean(corr_lkEc_eff_efftp1,1);
corr_hkEc_eff_efftp1 = nanmean(corr_hkEc_eff_efftp1,1);

corr_lkEc_en_efftp1 = nanmean(corr_lkEc_en_efftp1,1);
corr_hkEc_en_efftp1 = nanmean(corr_hkEc_en_efftp1,1);
for i_sub = 1:length(perf)
   tmp1(:,:,i_sub) = perf(i_sub).m_ratio; 
end
for i_run = 1:2
    for i_trial = 1:54
    tmp = corrcoef(kEc,tmp1(i_trial,i_run,:),'rows','pairwise');
    corr_kEc_eff(i_trial,i_run) = tmp(1,2);
    end
end
for i_sub = 1:size(low_kEc_efficiency,3)
    bin_low_kEc_eff(1,:,i_sub) = nanmean(low_kEc_efficiency(1:9,:,i_sub));
    bin_low_kEc_eff(2,:,i_sub) = nanmean(low_kEc_efficiency(10:18,:,i_sub));
    bin_low_kEc_eff(3,:,i_sub) = nanmean(low_kEc_efficiency(19:27,:,i_sub));
    bin_low_kEc_eff(4,:,i_sub) = nanmean(low_kEc_efficiency(28:36,:,i_sub));
    bin_low_kEc_eff(5,:,i_sub) = nanmean(low_kEc_efficiency(37:45,:,i_sub));
    bin_low_kEc_eff(6,:,i_sub) = nanmean(low_kEc_efficiency(46:54,:,i_sub));
    
    bin_low_kEc_en(1,:,i_sub) = nanmean(low_kEc_energization(1:9,:,i_sub));
    bin_low_kEc_en(2,:,i_sub) = nanmean(low_kEc_energization(10:18,:,i_sub));
    bin_low_kEc_en(3,:,i_sub) = nanmean(low_kEc_energization(19:27,:,i_sub));
    bin_low_kEc_en(4,:,i_sub) = nanmean(low_kEc_energization(28:36,:,i_sub));
    bin_low_kEc_en(5,:,i_sub) = nanmean(low_kEc_energization(37:45,:,i_sub));
    bin_low_kEc_en(6,:,i_sub) = nanmean(low_kEc_energization(46:54,:,i_sub));
    
    bin_low_kEc_perf(1,:,i_sub) = nanmean(low_kEc_performance(1:9,:,i_sub));
    bin_low_kEc_perf(2,:,i_sub) = nanmean(low_kEc_performance(10:18,:,i_sub));
    bin_low_kEc_perf(3,:,i_sub) = nanmean(low_kEc_performance(19:27,:,i_sub));
    bin_low_kEc_perf(4,:,i_sub) = nanmean(low_kEc_performance(28:36,:,i_sub));
    bin_low_kEc_perf(5,:,i_sub) = nanmean(low_kEc_performance(37:45,:,i_sub));
    bin_low_kEc_perf(6,:,i_sub) = nanmean(low_kEc_performance(46:54,:,i_sub));
    
    bin_low_kEc_HE(1,:,i_sub) = nanmean(low_kEc_HE(1:9,:,i_sub));
    bin_low_kEc_HE(2,:,i_sub) = nanmean(low_kEc_HE(10:18,:,i_sub));
    bin_low_kEc_HE(3,:,i_sub) = nanmean(low_kEc_HE(19:27,:,i_sub));
    bin_low_kEc_HE(4,:,i_sub) = nanmean(low_kEc_HE(28:36,:,i_sub));
    bin_low_kEc_HE(5,:,i_sub) = nanmean(low_kEc_HE(37:45,:,i_sub));
    bin_low_kEc_HE(6,:,i_sub) = nanmean(low_kEc_HE(46:54,:,i_sub));
end
for i_sub = 1:size(high_kEc_efficiency,3)
    bin_high_kEc_eff(1,:,i_sub) = nanmean(high_kEc_efficiency(1:9,:,i_sub));
    bin_high_kEc_eff(2,:,i_sub) = nanmean(high_kEc_efficiency(10:18,:,i_sub));
    bin_high_kEc_eff(3,:,i_sub) = nanmean(high_kEc_efficiency(19:27,:,i_sub));
    bin_high_kEc_eff(4,:,i_sub) = nanmean(high_kEc_efficiency(28:36,:,i_sub));
    bin_high_kEc_eff(5,:,i_sub) = nanmean(high_kEc_efficiency(37:45,:,i_sub));
    bin_high_kEc_eff(6,:,i_sub) = nanmean(high_kEc_efficiency(46:54,:,i_sub));
    
    bin_high_kEc_en(1,:,i_sub) = nanmean(high_kEc_energization(1:9,:,i_sub));
    bin_high_kEc_en(2,:,i_sub) = nanmean(high_kEc_energization(10:18,:,i_sub));
    bin_high_kEc_en(3,:,i_sub) = nanmean(high_kEc_energization(19:27,:,i_sub));
    bin_high_kEc_en(4,:,i_sub) = nanmean(high_kEc_energization(28:36,:,i_sub));
    bin_high_kEc_en(5,:,i_sub) = nanmean(high_kEc_energization(37:45,:,i_sub));
    bin_high_kEc_en(6,:,i_sub) = nanmean(high_kEc_energization(46:54,:,i_sub));
    
    bin_high_kEc_perf(1,:,i_sub) = nanmean(high_kEc_performance(1:9,:,i_sub));
    bin_high_kEc_perf(2,:,i_sub) = nanmean(high_kEc_performance(10:18,:,i_sub));
    bin_high_kEc_perf(3,:,i_sub) = nanmean(high_kEc_performance(19:27,:,i_sub));
    bin_high_kEc_perf(4,:,i_sub) = nanmean(high_kEc_performance(28:36,:,i_sub));
    bin_high_kEc_perf(5,:,i_sub) = nanmean(high_kEc_performance(37:45,:,i_sub));
    bin_high_kEc_perf(6,:,i_sub) = nanmean(high_kEc_performance(46:54,:,i_sub));
    
    bin_high_kEc_HE(1,:,i_sub) = nanmean(high_kEc_HE(1:9,:,i_sub));
    bin_high_kEc_HE(2,:,i_sub) = nanmean(high_kEc_HE(10:18,:,i_sub));
    bin_high_kEc_HE(3,:,i_sub) = nanmean(high_kEc_HE(19:27,:,i_sub));
    bin_high_kEc_HE(4,:,i_sub) = nanmean(high_kEc_HE(28:36,:,i_sub));
    bin_high_kEc_HE(5,:,i_sub) = nanmean(high_kEc_HE(37:45,:,i_sub));
    bin_high_kEc_HE(6,:,i_sub) = nanmean(high_kEc_HE(46:54,:,i_sub));
end

median(kEc)

figure()
plot(mean(squeeze(nanmean(low_kEc_efficiency,2)),2))
hold on
plot(mean(squeeze(nanmean(high_kEc_efficiency,2)),2))
legend('low kEc','high kEc')
xlabel('trial_i')
ylabel('efficiency score')
title('mean over both run efficiency average over two runs, over trials')

figure()
plot(mean(squeeze(nanmean(high_kEc_performance,2)),2))
hold on
plot(mean(squeeze(nanmean(low_kEc_performance,2)),2))
legend('high kEc','low kEc')
xlabel('trial_i')
ylabel('performance score in %age')
title('mean over both run performance over trials')

figure()
plot(nanmean(squeeze(bin_low_kEc_eff(:,1,:)),2))
hold on
plot(nanmean(squeeze(bin_low_kEc_eff(:,2,:)),2))
plot(nanmean(squeeze(bin_high_kEc_eff(:,1,:)),2))
plot(nanmean(squeeze(bin_high_kEc_eff(:,2,:)),2))
legend('low kEc run1','low kEc run2','high kEc run1','high kEc run2')
xlabel('trial_i')
ylabel('efficiency score in %age')
title('mean over both run performance over trials')

figure()
plot(nanmean(squeeze(bin_low_kEc_eff(:,1,:)),2))
hold on
plot(nanmean(squeeze(bin_low_kEc_eff(:,2,:)),2))
plot(nanmean(squeeze(bin_high_kEc_eff(:,1,:)),2))
plot(nanmean(squeeze(bin_high_kEc_eff(:,2,:)),2))
legend('low kEc run1','low kEc run2','high kEc run1','high kEc run2')
xlabel('trial_i')
ylabel('efficiency score in %age')
title('mean over both run performance over trials')

figure()
errorbar(nanmean(nanmean(squeeze(nanmean(bin_low_kEc_perf(:,:,:),2)),2),2),std(squeeze(nanmean(bin_low_kEc_perf,2)),0,2))
hold on
errorbar(nanmean(nanmean(squeeze(nanmean(bin_high_kEc_perf(:,2,:),2)),2),2),std(squeeze(nanmean(bin_high_kEc_perf,2)),0,2))
legend('low kEc','high kEc')
xlabel('bin_i composed of 9 trials')
ylabel('perf score in %age')
title('mean over both run performance over trials')

figure()
errorbar(nanmean(nanmean(squeeze(nanmean(bin_low_kEc_eff(:,:,:),2)),2),2),std(squeeze(nanmean(bin_low_kEc_eff,2)),0,2))
hold on
errorbar(nanmean(nanmean(squeeze(nanmean(bin_high_kEc_eff(:,:,:),2)),2),2),std(squeeze(nanmean(bin_high_kEc_eff,2)),0,2))
legend('low kEc','high kEc')
xlabel('bin_i composed of 9 trials')
ylabel('efficiency score')
title('mean over both run efficiency over trials')

figure()
errorbar(nanmean(nanmean(squeeze(nanmean(bin_low_kEc_HE(:,:,:),2)),2),2),std(squeeze(nanmean(bin_low_kEc_HE,2)),0,2))
hold on
errorbar(nanmean(nanmean(squeeze(nanmean(bin_high_kEc_HE(:,:,:),2)),2),2),std(squeeze(nanmean(bin_high_kEc_HE,2)),0,2))
legend('low kEc','high kEc')
xlabel('bin_i composed of 9 trials')
ylabel('Number of HE selected')
title('mean over both run choices')

figure()
errorbar(nanmean(nanmean(squeeze(nanmean(bin_low_kEc_en(:,:,:),2)),2),2),std(squeeze(nanmean(bin_low_kEc_en,2)),0,2))
hold on
errorbar(nanmean(nanmean(squeeze(nanmean(bin_high_kEc_en(:,:,:),2)),2),2),std(squeeze(nanmean(bin_high_kEc_en,2)),0,2))
legend('low kEc','high kEc')
xlabel('bin_i composed of 9 trials')
ylabel('Energization score')
title('mean over both run choices')

figure()
boxchart(nanmean(low_kEc_std_HE,2))
hold on
boxchart(nanmean(high_kEc_std_HE,2))
legend('low kEc','high kEc')
xlabel('boxplot')
ylabel('fluctuation in HE (sd) score')
title('fluctuation (sd) in HE chosen depending on the groups kEc')

figure()
plot(corr_kEc_eff(:,1))
hold on
plot(corr_kEc_eff(:,2))
legend('run1','run2')
xlabel('tria_i')
ylabel('corr score')
title('correlation at each trial, between kEc and efficacy across participants, averaged over both run')

end

for i_sub =1:69
[binned_data1] = bin_data(perf(i_sub).AUC_overshoot(1:54,1),9)
[binned_data2] = bin_data(perf(i_sub).AUC_overshoot(1:54,2),9)
tmp(i_sub) = nansum(nanmean(perf(i_sub).AUC,2));
avg_AUC(:,i_sub) = (binned_data1+binned_data2)/2;
mdl = fitlm(1:9,avg_AUC(:,i_sub))
fit_time(i_sub) = table2array(mdl.Coefficients(2,1));
Phi = sensitivitiesPhi{5,i_sub};
kFp(i_sub) = log(1+exp(Phi(6)));
end
for i = 1:size(fit_time,1)
        fit_time(i,fit_time(i,:) >= 3*nanstd(fit_time(i,:))+nanmedian(fit_time(i,:))|fit_time(i,:) <= -3*nanstd(fit_time(i,:))+nanmedian(fit_time(i,:))) = NaN;
end
for i = 1:size(kFp,1)
        kFp(i,kFp(i,:) >= 3*nanstd(kFp(i,:))+nanmedian(kFp(i,:))|kFp(i,:) <= -3*nanstd(kFp(i,:))+nanmedian(kFp(i,:))) = NaN;
end
[h,p] = scatter_perso(fit_time',tmp','latency to start effort','kFp','Correlation')
clear kFp
clear fit_time

mdl = fitlm(1:6,avg_AUC(:,1))
figure()
plot(avg_AUC(:,1))
hold on
plot(mdl)

%% Save all computed datasets
% do you want to save the results ?
if do_save == true
    save_data(sensitivitiesPhi,sensitivitiesTheta,CID_nb,p,m,RTP,RTM,all_MVC,all_NMP,physical_IP,mental_IP)
    % for i_sub = 1:69
    %     ratio_D(i_sub,:,:) = perf(i_sub).ratio_simu.ratio_D;
    %     ratio_LE(i_sub,:,:) = perf(i_sub).ratio_simu.ratio_LE;
    %     ratio_ME(i_sub,:,:) = perf(i_sub).ratio_simu.ratio_ME;
    %     ratio_HE(i_sub,:,:) = perf(i_sub).ratio_simu.ratio_HE;
    % end
    % ratio_D = nanmean(squeeze(nanmean(ratio_D)),2);
    % ratio_LE = nanmean(squeeze(nanmean(ratio_LE)),2);
    % ratio_ME = nanmean(squeeze(nanmean(ratio_ME)),2);
    % ratio_HE = nanmean(squeeze(nanmean(ratio_HE)),2);
    % ratio = [ratio_D,ratio_LE,ratio_ME,ratio_HE];
    % save('M_ratio_for_Simu.mat','ratio')
else % transform parameters into correct space (for model 5)
    
    model_i = 5;
    mod = [];
    
    for i = 1:length(sensitivitiesPhi)
        % raw parameters
        %     mod.kR(i) =  sensitivitiesPhi{model_i,i}(1);
        %     mod.kP(i) =  sensitivitiesPhi{model_i,i}(2);
        %     mod.kEp(i) =  sensitivitiesPhi{model_i,i}(3);
        %     mod.kEm(i) =  sensitivitiesPhi{model_i,i}(4);
        %     mod.biais(i) = sensitivitiesPhi{model_i,i}(5);
        %     mod.kFp(i) =  sensitivitiesPhi{model_i,i}(6);
        %     mod.kLm(i) =  sensitivitiesPhi{model_i,i}(7);
        
        % parameters with positivity constraint (need to be retransformed)
        mod.kR(i) =  log(1+exp(sensitivitiesPhi{model_i,i}(1)));
        mod.kP(i) =  log(1+exp(sensitivitiesPhi{model_i,i}(2)));
        mod.kEp(i) =  log(1+exp(sensitivitiesPhi{model_i,i}(3)));
        mod.kEm(i) =  log(1+exp(sensitivitiesPhi{model_i,i}(4)));
        mod.biais(i) = sensitivitiesPhi{model_i,i}(5);
        mod.kFp(i) =  log(1+exp(sensitivitiesPhi{model_i,i}(6)));
        mod.kLm(i) =  log(1+exp(sensitivitiesPhi{model_i,i}(7)));
        
    end % loop over parameters
    
end
cd(results_folder);


end % function