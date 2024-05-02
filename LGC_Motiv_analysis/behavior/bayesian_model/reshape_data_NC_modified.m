function [var,all_choices,all_choices_rescaled,deltaRP_idx,deltaRP,deltaE,RT_physical,RT_mental,all_MVC,all_NMP,physical_IP,mental_IP,perf,...
    Incentive_idx] = reshape_data_NC_modified(i_sub,all_data,rescale_choice,binary_answers,p_or_m, include_currEff, include_prevEff)


% empty and prepare another var, in case it's size was not the same between subject (miss a session)
var = NaN(7,54,4);
% prepare also array of all the responses from participants, and it's increment
tmp_choices = [];
RT_physical = [];
RT_mental = [];
% add initial calibration to all calibrations
all_MVC = [all_data.(['sub_',num2str(i_sub)]).physicalCalib.MVC];
all_NMP = [all_data.(['sub_',num2str(i_sub)]).mentalCalib.NMP];
physical_IP = all_data.(['sub_',num2str(i_sub)]).IP.IP_variables.physicalDeltaIP;
mental_IP = all_data.(['sub_',num2str(i_sub)]).IP.IP_variables.mentalDeltaIP;
p_run_idx = 1;
m_run_idx = 1;
AUC = [];
for i_run = 1:4
    RP_trials = NaN(1,54);
    
    % Check run order and if it is physical, otherwise mental
    if ((i_run ==1 || i_run == 3) && strcmp('p',p_or_m(i_sub))==1) || ((i_run ==2 || i_run == 4) && strcmp('p',p_or_m(i_sub))==0)
        
        % extract all MVC
        all_MVC = cat(1,all_MVC,[all_data.(['sub_',num2str(i_sub)]).(['run_',num2str(i_run)]).start_maxPerf.MVC.MVC;all_data.(['sub_',num2str(i_sub)]).(['run_',num2str(i_run)]).end_maxPerf.MVC.MVC]);
        
        % Default LR, if -1 default = gauche.
        %Choice, if -1 = gauche
        E_chosen = all_data.(['sub_',num2str(i_sub)]).(['run_',num2str(i_run)]).physicalPerf.E_chosen;
        % always the same in all 4 runs
        Incentive_idx = abs(all_data.(['sub_',num2str(i_sub)]).(['run_',num2str(i_run)]).choice_opt.R.left-all_data.(['sub_',num2str(i_sub)]).(['run_',num2str(i_run)]).choice_opt.R.right);
        
        default_LR = all_data.(['sub_',num2str(i_sub)]).(['run_',num2str(i_run)]).physicalPerf.choiceOptions.default_LR;
        R_or_P = all_data.(['sub_',num2str(i_sub)]).(['run_',num2str(i_run)]).physicalPerf.choiceOptions.R_or_P;
        RT_physical =cat(1,RT_physical,all_data.(['sub_',num2str(i_sub)]).(['run_',num2str(i_run)]).physicalPerf.durations.dispChoiceOptions);
        % deltaR and delta E taking into account where the default option is. Like this, we
        % compute the proba of taking the big value, always.
        
        %multiplied by 100 to be in cents
        deltaRP = (all_data.(['sub_',num2str(i_sub)]).(['run_',num2str(i_run)]).physicalPerf.choiceOptions.monetary_amount.left - all_data.(['sub_',num2str(i_sub)]).(['run_',num2str(i_run)]).physicalPerf.choiceOptions.monetary_amount.right).*default_LR*100;
        deltaE = (all_data.(['sub_',num2str(i_sub)]).(['run_',num2str(i_run)]).physicalPerf.choiceOptions.E.left - all_data.(['sub_',num2str(i_sub)]).(['run_',num2str(i_run)]).physicalPerf.choiceOptions.E.right).*default_LR;
        Ep_or_Em_trials = ones(54,1)';
        tmp_choices = all_data.(['sub_',num2str(i_sub)]).(['run_',num2str(i_run)]).physicalPerf.choice.*default_LR*(-1);
        
        [perf.AUC(:,p_run_idx),perf.AUC_overshoot(:,p_run_idx),perf.force_peak(:,p_run_idx),perf.force_steepness(:,p_run_idx),...
            perf.p_latency(:,p_run_idx),~,t(:,p_run_idx),perf.p_performance(:,p_run_idx)] = extract_physical_perf(all_data.(['sub_',num2str(i_sub)]).(['run_',num2str(i_run)]).physicalPerf);
        
        if p_run_idx == 1
            t_post1 = t(p_run_idx).time_post_effort(:);
            t_pre1 = t(p_run_idx).time_pre_effort(:);
            p_AUC1(:) =perf.AUC(:,p_run_idx);
            p_force1(:) = perf.force_peak(:,p_run_idx)/50;
            p_run_idx = p_run_idx+1;
        elseif p_run_idx == 2
            t_post2 = t(p_run_idx).time_post_effort(:);
            t_pre2 = t(p_run_idx).time_pre_effort(:);
            p_AUC2(:) =perf.AUC(:,p_run_idx);
            p_force2(:) = perf.force_peak(:,p_run_idx)/50;
            p_run_idx = p_run_idx+1;
        end
        
    else
        % if subject is sub
        if i_sub == 26 && m_run_idx == 2
            % extract all NMP
            all_NMP = cat(1,all_NMP,[NaN;NaN]);
            E_chosen =  NaN(1,54);
            Incentive_idx = NaN(1,54);
            default_LR =  NaN(1,54);
            R_or_P =  NaN(1,54);
            RT_mental =  NaN(1,54);
            deltaRP =  NaN(1,54);
            deltaE =  NaN(1,54);
            Ep_or_Em_trials =  NaN(1,54);
            tmp_choices =  NaN(1,54);
            perf.NB_answer(:,m_run_idx) =   NaN(1,54);
            perf.NB_correct_answer(:,m_run_idx) =   NaN(1,54);
            perf.m_errors(:,m_run_idx) =   NaN(1,54);
            perf.m_performance(:,m_run_idx) =   NaN(1,54);
            perf.m_avg_rt(:,m_run_idx) =  NaN(1,54);
            perf.m_latency(:,m_run_idx) =   NaN(1,54);
            perf.m_ratio(:,m_run_idx) = NaN(1,54);
            m_AUC2(:) =  NaN(1,54);
            m_ratio2 = NaN(1,54);
            m_run_idx = m_run_idx +1;
            
        else
            % extract all NMP
            all_NMP = cat(1,all_NMP,[all_data.(['sub_',num2str(i_sub)]).(['run_',num2str(i_run)]).start_maxPerf.n_maxPerf;all_data.(['sub_',num2str(i_sub)]).(['run_',num2str(i_run)]).end_maxPerf.n_maxPerf]);
            
            % Default LR, if -1 default = gauche.
            %Choice, if -1 = gauche
            E_chosen = all_data.(['sub_',num2str(i_sub)]).(['run_',num2str(i_run)]).mentalE_perf.E_chosen;
            Incentive_idx = abs(all_data.(['sub_',num2str(i_sub)]).(['run_',num2str(i_run)]).choice_opt.R.left-all_data.(['sub_',num2str(i_sub)]).(['run_',num2str(i_run)]).choice_opt.R.right);

            default_LR = all_data.(['sub_',num2str(i_sub)]).(['run_',num2str(i_run)]).mentalE_perf.choiceOptions.default_LR;
            % prepare all the var informations
            R_or_P = all_data.(['sub_',num2str(i_sub)]).(['run_',num2str(i_run)]).mentalE_perf.choiceOptions.R_or_P;
            RT_mental = cat(1,RT_mental,all_data.(['sub_',num2str(i_sub)]).(['run_',num2str(i_run)]).mentalE_perf.durations.dispChoiceOptions);
            %multiplied by 100 to be in cents
            deltaRP = (all_data.(['sub_',num2str(i_sub)]).(['run_',num2str(i_run)]).mentalE_perf.choiceOptions.monetary_amount.left - all_data.(['sub_',num2str(i_sub)]).(['run_',num2str(i_run)]).mentalE_perf.choiceOptions.monetary_amount.right).*default_LR*100;
            deltaE = (all_data.(['sub_',num2str(i_sub)]).(['run_',num2str(i_run)]).mentalE_perf.choiceOptions.E.left - all_data.(['sub_',num2str(i_sub)]).(['run_',num2str(i_run)]).mentalE_perf.choiceOptions.E.right).*default_LR;
            Ep_or_Em_trials = ones(54,1)'-1;
            tmp_choices = (all_data.(['sub_',num2str(i_sub)]).(['run_',num2str(i_run)]).mentalE_perf.choice).*default_LR*(-1);
            
            [perf.NB_answer(:,m_run_idx),perf.NB_correct_answer(:,m_run_idx),perf.m_errors(:,m_run_idx),perf.m_performance(:,m_run_idx),perf.m_avg_rt(:,m_run_idx),perf.m_latency(:,m_run_idx)...
                ,perf.mental_rt(:,m_run_idx),total_time(:,m_run_idx),perf.m_ratio(:,m_run_idx)] = extract_mental_perf(all_data.(['sub_',num2str(i_sub)]).(['run_',num2str(i_run)]).mentalE_perf.perfSummary);

                       
            if m_run_idx == 1
                m_AUC1(:) = perf.NB_answer(:,m_run_idx)/all_NMP(1);
                m_ratio1(:) = perf.NB_correct_answer(:,m_run_idx)./total_time(:,m_run_idx);
                m_run_idx = m_run_idx +1;
            elseif m_run_idx == 2
                m_AUC2(:) = perf.NB_answer(:,m_run_idx)/all_NMP(1);
                m_ratio2(:) = perf.NB_correct_answer(:,m_run_idx)./total_time(:,m_run_idx);
                m_run_idx = m_run_idx +1;
            end
            
        end
    end
    
    % if we rescale choice between 0 and 1
    if rescale_choice == true
        for i_choice = 1:length(tmp_choices)
            switch tmp_choices(i_choice)
                
                case -2
                    tmp_choices(i_choice) = 0;
                case -1
                    tmp_choices(i_choice) = 0.25;
                case 0
                    tmp_choices(i_choice) = NaN;
                case 1
                    tmp_choices(i_choice) = 0.75;
                case 2
                    tmp_choices(i_choice) = 1;
                    
            end
        end
    end
    
    % Modify cell into classical array
    for iTrial = 1:54
        if i_sub == 26 && i_run == 4
        else
            switch R_or_P{iTrial}
                case 'R'
                    RP_trials(iTrial) = 1;
                case 'P'
                    RP_trials(iTrial) = 0;
            end
            
        end
    end
    
    % prepare the sum of effort a t-1
    sum_E_chosen=[];
    for i = 1:54
        if i == 1
            sum_E_chosen(1) = 0;
        else
            sum_E_chosen(i) = sum_E_chosen(i-1) +E_chosen(i-1);
        end
    end
    % prepare the matrix of our variables, for each sessions
    var(:,:,i_run) = [deltaRP; deltaE; RP_trials; Ep_or_Em_trials;(1:54)-1;[0 E_chosen(1:end-1)];sum_E_chosen/54];
    all_choices(1+(i_run-1)*54:(i_run)*54)=tmp_choices;
    

    
end


%for the simulation. compute deltaE
[perf.ratio_simu] = compute_ratio_over_time_for_simu(perf.m_ratio,E_chosen);

% sum elements for g_obs, shift of 1 as we use t-1 effort to compute sv(t)
for i = 1:54
    if i == 1
        p_AUC1_sum(1) = 0;
        p_AUC2_sum(1) = 0;
        time_post1_sum(1) = t_post1(1);
        time_post2_sum(1) = t_post2(1);
        time_pre1_sum(1) = t_pre1(1);
        time_pre2_sum(1) = t_pre2(1);
        m_AUC1_sum(1) = 0;
        m_AUC2_sum(1) = 0;
        m_ratio1_sum(1) = 0;
        m_ratio2_sum(2) = 0;
    else
        p_AUC1_sum(i) = p_AUC1_sum(i-1) + p_AUC1(i-1);
        p_AUC2_sum(i) = p_AUC2_sum(i-1) + p_AUC2(i-1);
        time_post1_sum(i) = time_post1_sum(i-1) + t_post1(i);
        time_post2_sum(i) = time_post2_sum(i-1) + t_post2(i);
        time_pre1_sum(i) = time_pre1_sum(i-1) + t_pre1(i);
        time_pre2_sum(i) = time_pre2_sum(i-1) + t_pre2(i);
        m_AUC1_sum(i) = m_AUC1_sum(i-1) + m_AUC1(i-1);
        m_AUC2_sum(i) = m_AUC2_sum(i-1) + m_AUC2(i-1);
        m_ratio1_sum(i) = m_ratio1_sum(i-1) + m_ratio1(i-1); 
        m_ratio2_sum(i) = m_ratio2_sum(i-1) + m_ratio2(i-1);
        
    end
end

if include_currEff == 1 && include_prevEff == 0 % current efficacy (for current trial)
    m_ratio1 = m_ratio1;
    m_ratio2 = m_ratio2;
    p_force1 = p_force1;
    p_force2 = p_force2;
    p_AUC1 = p_AUC1/100;
    p_AUC2 = p_AUC2/100;
elseif include_currEff == 0 && include_prevEff == 1 % previous efficacy (previous trial)
    m_ratio1 = [0 m_ratio1(1:end-1)];
    m_ratio2 = [0 m_ratio2(1:end-1)];
    p_force1 = [0 p_force1(1:end-1)];
    p_force2 = [0 p_force2(1:end-1)];
    p_AUC1 = [0 p_AUC1(1:end-1)]/100;
    p_AUC2 = [0 p_AUC2(1:end-1)]/100;
end
%% reshape it as one big session
var = reshape(permute(var,[2 3 1]),size(var,2)*4,[])';
% prepare measure of effort performed.
if strcmp('p',p_or_m(i_sub))==1
    var(8,1:54) = p_AUC1_sum/1000;
    var(8,55:108) = m_AUC1_sum/10;
    var(8,109:162) = p_AUC2_sum/1000;
    var(8,163:216) = m_AUC2_sum/10;
else
    var(8,1:54) = m_AUC1_sum/10;
    var(8,55:108) = p_AUC1_sum/1000;
    var(8,109:162) = m_AUC2_sum/10;
    var(8,163:216) = p_AUC2_sum/1000;
end
if strcmp('p',p_or_m(i_sub))==1
%     [data1,data2,~,data1_idx,data2_idx,~] = split_data(perf.AUC(:,1)',all_choices(1:54)','double',0.5);
%     mean(data1)
%     figure()
%     subplot(1,2,1)
%     hist(data1)
%     hold on
%     [N,X] = hist(data2)
%     xlabel('AUC value')
%     ylabel('nb_of AUC')
%     Bh = bar(X,N,'facecolor',[0.3 0.7 0.2]);
%     [data1,data2,~,data1_idx,data2_idx,~] = split_data(perf.AUC(:,2)',all_choices(1:54)','double',0.5);
%     mean(data1)
%     subplot(1,2,2)
%     xlabel('AUC value')
%     ylabel('nb_of AUC')
%     hist(data1)
%     hold on
%     [N,X] = hist(data2)
%     Bh = bar(X,N,'facecolor',[0.3 0.7 0.2]);

    var(9,1:54) = p_AUC1_sum/1000;
    var(9,55:108) = m_ratio1;
    var(9,109:162) = p_AUC2_sum/1000;
    var(9,163:216) = m_ratio2;
else
    var(9,1:54) = m_ratio1;
    var(9,55:108) = p_AUC1_sum/1000;
    var(9,109:162) = m_ratio2;
    var(9,163:216) = p_AUC2_sum/1000;
end
if strcmp('p',p_or_m(i_sub))==1
    var(10,1:54) = p_AUC1;
    var(10,55:108) = m_ratio1;
    var(10,109:162) = p_AUC2;
    var(10,163:216) = m_ratio2;
else
    var(10,1:54) = zeros(1,54);
    var(10,55:108) = p_AUC1;
    var(10,109:162) = zeros(1,54);
    var(10,163:216) = p_AUC2;
end
if strcmp('p',p_or_m(i_sub))==1
    var(11,1:54) = t_post1;
    var(11,55:108) = zeros(1,54);
    var(11,109:162) = t_post2;
    var(11,163:216) = zeros(1,54);
else
    var(11,1:54) = zeros(1,54);
    var(11,55:108) = t_post1;
    var(11,109:162) = zeros(1,54);
    var(11,163:216) = t_post2;
end
if strcmp('p',p_or_m(i_sub))==1
    var(12,1:54) = t_pre1;
    var(12,55:108) = zeros(1,54);
    var(12,109:162) = t_pre2;
    var(12,163:216) = zeros(1,54);
else
    var(12,1:54) = zeros(1,54);
    var(12,55:108) = t_pre1;
    var(12,109:162) = zeros(1,54);
    var(12,163:216) = t_pre2;
end
var(9,var(9,:) == Inf) = 0;
var(13,:) = all_choices;

var(1,:)=abs(var(1,:));
%look at the matrix of choices from participants
deltaRP_idx = unique(deltaRP);
% two cases occurs, sometimes reward was too high, or too small, meaning with only had 4 different
% values presented to participant, not 6. The rest of the code requires deltaRP to be of length 6
if  deltaRP_idx(1) == -0.2 && length(deltaRP_idx) == 4
    deltaRP_idx = [-0.0200 -0.0200 -0.0100 0.0100 0.0200 0.0200];
elseif length(deltaRP_idx) == 4
    %  vector for CID 064 has somekind of error in it. I have to manually recode deltRP and deltaRP_idx
    deltaRP_idx = [-2.000 -2.000 -1.000 1.000 2.000 2.000];
    deltaRP =[2.00000000000000,1.00000000000000,2.00000000000000,-1.00000000000000,1.00000000000000,-2.00000000000000,...
        -2.00000000000000,2.00000000000000,-2.00000000000000,-2.00000000000000,-2.00000000000000,-2.00000000000000,...
        2.00000000000000,2.00000000000000,-1.00000000000000,1.00000000000000,2.00000000000000,-1.00000000000000,2.00000000000000,...
        1.00000000000000,2.00000000000000,2.00000000000000,-2.00000000000000,-2.00000000000000,-2.00000000000000,-1.00000000000000,...
        1.00000000000000,-2.00000000000000,-2.00000000000000,2.00000000000000,2.00000000000000,-2.00000000000000,-1.00000000000000,...
        -1.00000000000000,2.00000000000000,1.00000000000000,-1.00000000000000,2.00000000000000,2.00000000000000,1.00000000000000,...
        -2.00000000000000,2.00000000000000,1.00000000000000,-1.00000000000000,-2.00000000000000,2.00000000000000,-2.00000000000000,...
        2.00000000000000,-2.00000000000000,-1.00000000000000,-2.00000000000000,2.00000000000000,1.00000000000000,-2.00000000000000];
end

% in case we want the answer as a binary.
all_choices_rescaled = [];
if binary_answers == true
    for i_choice = 1:length(all_choices)
        if all_choices(i_choice) >=0 && all_choices(i_choice) < 0.25
            all_choices_rescaled(i_choice) = 0;
        elseif all_choices(i_choice) >=0.25 && all_choices(i_choice) < 0.5
            all_choices_rescaled(i_choice) = 0;
        elseif all_choices(i_choice) >=0.5 && all_choices(i_choice) < 0.75
            all_choices_rescaled(i_choice) = 1;
        elseif all_choices(i_choice) >=0.75 && all_choices(i_choice) <= 1
            all_choices_rescaled(i_choice) = 1;
        elseif isnan(all_choices(i_choice))
            all_choices_rescaled(i_choice) = NaN;
        end
    end
end
end % function