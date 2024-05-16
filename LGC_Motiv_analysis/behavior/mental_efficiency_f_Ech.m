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
% r_corr: correlation coefficient for each subject for mental efficiency in
% function of the effort chosen
%
% pval: p.value for the t.test looking at whether mental efficiency varies
% with the effort chosen

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
[mental_eff, Ech,...
    mental_eff_fit] = deal(NaN(nTotalTrials,NS));
n_hE_levels = 4;
E_levels = [0,1,2,3];
[mental_eff_bin, mental_eff_fit_bin] = deal(NaN(n_hE_levels, NS));
r_corr = deal(NaN(1,NS));

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
        [E_chosen_tmp] = extract_E_chosen(subBehaviorFolder, sub_nm, run_nm, 'mental');
        mental_eff(run_trials_idx, iS) = efficacy_ter_with2first_tmp.allTrials;
        Ech(run_trials_idx, iS) = E_chosen_tmp;
    end % run loop
    
    %% perform correlation for each subject
    [r_corr(iS), ~, ~, mental_eff_fit(:,iS)] = glm_package(Ech(:, iS), mental_eff(:, iS), 'normal', 'on');
    
    %% extract bins
    for iE = E_levels
        jE = iE + 1;
        curr_Ech_idx = Ech(:,iS) == iE;
        if sum(curr_Ech_idx) > 0
            mental_eff_bin(jE,iS) = mean(mental_eff(curr_Ech_idx, iS),'omitnan');
            mental_eff_fit_bin(jE,iS) = mean(mental_eff_fit(curr_Ech_idx, iS),'omitnan');
        end
    end % effort chosen level loop
end % loop through subjects

%% test correlation
[~,pval] = ttest(r_corr);

%% average across subjects
[mental_eff_bin_avg, mental_eff_bin_sem] = mean_sem_sd(mental_eff_bin, 2);
[mental_eff_fit_bin_avg, mental_eff_fit_bin_sem] = mean_sem_sd(mental_eff_fit_bin, 2);

%% figure display
if fig_disp == 1
    fig;
    er_hdl = errorbar(E_levels, mental_eff_bin_avg, mental_eff_bin_sem);
    errorbar_hdl_upgrade(er_hdl);
    fit_hdl = plot(E_levels, mental_eff_fit_bin_avg);
    fit_hdl_upgrade(fit_hdl);
    xticks(0:3);
    xticklabels({'0','1','2','3'});
    xlabel('E chosen');
    ylabel('Mental Efficiency');
    xlim([-0.01 3.1]);
end % figure display

end % function