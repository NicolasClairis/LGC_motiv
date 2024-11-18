function[r_corr, pval] = mental_efficiency_f_Ech(fig_disp)
% [r_corr, pval] = mental_efficiency_f_Ech(fig_disp)
% mental_efficiency_f_Ech tests whether the mental efficiency (defined as
% the number of correct answers in the mental task divided by the total
% time taken to finish the performance in the current trial (including the
% two first useless answers) varies with the level of effort chosen across
% trials.
%
% INPUTS
% fig_disp: display figure (1) or not (0)
%
% OUTPUTS
% r_corr: structure with correlation coefficient for each subject for mental efficiency in
% function of the effort chosen (r_corr.mentalEff_f_Ech) and of the choice
% made (r_corr.mentalEff_f_choice)
%
% pval: structure with p.value for the t.test looking at whether mental efficiency varies
% with the effort chosen (r_corr.mentalEff_f_Ech) and of the choice
% made (r_corr.mentalEff_f_choice)

%% default for inputs
if ~exist('fig_disp','var') || isempty(fig_disp) || ~ismember(fig_disp,[0,1])
   fig_disp = 1; 
end

%% subject selection
[study_nm, condition, ~, subject_id, NS] = sub_id;

%% working dir
rootPath = ['E:',filesep,study_nm,filesep];

%% initialize variables of interest
nTrialsPerRun = 54;
nMentalRuns = 2;
nTotalTrials = nTrialsPerRun*nMentalRuns;
[mental_eff, Ech, choice_hE,...
    mental_eff_Ech_fit, mental_eff_ch_fit] = deal(NaN(nTotalTrials,NS));
n_hE_levels = 4;
E_levels = [0,1,2,3];
[mental_eff_Ech_bin, mental_eff_Ech_fit_bin] = deal(NaN(n_hE_levels, NS));
n_ch_levels = 4; % 0/0.25/0.75/1
ch_levels = [0, 0.25, 0.75, 1];
[mental_eff_ch_bin, mental_eff_ch_fit_bin] = deal(NaN(n_ch_levels, NS));
[r_corr.mentalEff_f_Ech, r_corr.mentalEff_f_choice] = deal(NaN(1,NS));

%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    subBehaviorFolder = [fullfile(rootPath,['CID',sub_nm],'behavior'),filesep];
    
    [runs] = runs_definition(study_nm, sub_nm, condition);
    n_Em_runs = runs.nb_runs.Em;
    kRun = 0;
    for iRun = 1:n_Em_runs
        jRun = runs.Em.runsToKeep(iRun);
        run_nm = num2str(jRun);
        kRun = kRun + 1;
        run_trials_idx = (1:nTrialsPerRun) + nTrialsPerRun.*(kRun - 1);
        % extract data of interest
        % mental efficiency current run
        [~, ~, ~, RT_avg_tmp,...
            RTtotal_with2firstUseless_tmp,...
            RTtotal_pureNback_tmp,...
            efficacy_with2first_tmp,...
            efficacy_pureNback_tmp,...
            efficacy_bis_with2first_tmp,...
            efficacy_bis_pureNback_tmp, latency_tmp,...
            efficacy_ter_with2first_tmp,...
            efficacy_ter_pureNback_tmp,...
            RTtotal_with2firstUseless2_tmp,...
            RTtotal_pureNback2_tmp] = extract_mental_perf(subBehaviorFolder, sub_nm, run_nm);
        % E chosen current run
        [E_chosen_tmp] = extract_E_chosen(subBehaviorFolder, sub_nm, run_nm, 'mental');
        % choice current run
        [choice_highE_tmp] = extract_choice_hE_bis(subBehaviorFolder, sub_nm, run_nm, 'mental');
        
        % load data into global vectors
        mental_eff(run_trials_idx, iS) = efficacy_ter_with2first_tmp.allTrials;
        Ech(run_trials_idx, iS) = E_chosen_tmp;
        choice_hE(run_trials_idx, iS) = choice_highE_tmp;
    end % run loop
    
    %% perform correlation for each subject
    [r_corr.mentalEff_f_Ech(iS), ~, ~, mental_eff_Ech_fit(:,iS)] = glm_package(Ech(:, iS), mental_eff(:, iS), 'normal', 'on');
    [r_corr.mentalEff_f_choice(iS), ~, ~, mental_eff_ch_fit(:,iS)] = glm_package(choice_hE(:, iS), mental_eff(:, iS), 'normal', 'on');
    
    %% extract bins
    % Mental Eff = f(Ech)
    for iE = E_levels
        jE = iE + 1;
        curr_Ech_idx = Ech(:,iS) == iE;
        if sum(curr_Ech_idx) > 0
            mental_eff_Ech_bin(jE,iS) = mean(mental_eff(curr_Ech_idx, iS),'omitnan');
            mental_eff_Ech_fit_bin(jE,iS) = mean(mental_eff_Ech_fit(curr_Ech_idx, iS),'omitnan');
        end
    end % effort chosen level loop
    
    % Mental Eff = f(choice)
    for iCh_lvl = 1:n_ch_levels
        ch_lvl_tmp = ch_levels(iCh_lvl);
        curr_ch_idx = choice_hE(:,iS) == ch_lvl_tmp;
        if sum(curr_ch_idx) > 0
            mental_eff_ch_bin(iCh_lvl,iS) = mean(mental_eff(curr_ch_idx, iS),'omitnan');
            mental_eff_ch_fit_bin(iCh_lvl,iS) = mean(mental_eff_ch_fit(curr_ch_idx, iS),'omitnan');
        end
    end % effort chosen level loop
end % loop through subjects

%% test correlations
[~,pval.mentalEff_f_Ech] = ttest(r_corr.mentalEff_f_Ech);
[~,pval.mentalEff_f_choice] = ttest(r_corr.mentalEff_f_choice);

%% average across subjects
[mental_eff_Ech_bin_avg, mental_eff_Ech_bin_sem] = mean_sem_sd(mental_eff_Ech_bin, 2);
[mental_eff_Ech_fit_bin_avg, mental_eff_Ech_fit_bin_sem] = mean_sem_sd(mental_eff_Ech_fit_bin, 2);
[mental_eff_ch_bin_avg, mental_eff_ch_bin_sem] = mean_sem_sd(mental_eff_ch_bin, 2);
[mental_eff_ch_fit_bin_avg, mental_eff_ch_fit_bin_sem] = mean_sem_sd(mental_eff_ch_fit_bin, 2);

%% figure display
if fig_disp == 1
    fig;
    
    % mental efficiency = f(E chosen level)
    subplot(1,2,1); hold on;
    er_hdl = errorbar(E_levels, mental_eff_Ech_bin_avg, mental_eff_Ech_bin_sem);
    errorbar_hdl_upgrade(er_hdl);
    fit_hdl = plot(E_levels, mental_eff_Ech_fit_bin_avg);
    fit_hdl_upgrade(fit_hdl);
    xticks(E_levels);
    xticklabels({'0','1','2','3'});
    xlabel('E chosen');
    ylabel('Mental Efficiency');
    xlim([-0.01 3.1]);
    place_r_and_pval(mean(r_corr.mentalEff_f_Ech,2,'omitnan'), pval.mentalEff_f_Ech);
    
    % mental efficiency = f(choice)
    subplot(1,2,2); hold on;
    er_hdl = errorbar(ch_levels, mental_eff_ch_bin_avg, mental_eff_ch_bin_sem);
    errorbar_hdl_upgrade(er_hdl);
    fit_hdl = plot(ch_levels, mental_eff_ch_fit_bin_avg);
    fit_hdl_upgrade(fit_hdl);
    xticks(ch_levels);
    xticklabels({'0','0.25','0.75','1'});
    xlabel('Choice');
    ylabel('Mental Efficiency');
    xlim([0 1]);
    place_r_and_pval(mean(r_corr.mentalEff_f_choice,2,'omitnan'), pval.mentalEff_f_choice);
end % figure display

end % function