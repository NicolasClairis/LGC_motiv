function[r_corr, pval] = prev_mental_efficiency_f_Ech(fig_disp)
% [r_corr, pval] = prev_mental_efficiency_f_Ech(fig_disp)
% prev_mental_efficiency_f_Ech tests whether the previous mental efficiency (defined as
% the number of correct answers in the mental task divided by the total
% time taken to finish the performance in the current trial (including the
% two first useless answers) varies with the level of effort chosen across
% trials.
%
% INPUTS
% fig_disp: display figure (1) or not (0)
%
% OUTPUTS
% r_corr: correlation coefficient for each subject for previous mental efficiency in
% function of the effort chosen
%
% pval: p.value for the t.test looking at whether previous mental efficiency varies
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
[prev_mental_eff, Ech,...
    prev_mental_eff_fit] = deal(NaN(nTotalTrials,NS));
n_hE_levels = 4;
E_levels = [0,1,2,3];
[prev_mental_eff_bin, prev_mental_eff_fit_bin] = deal(NaN(n_hE_levels, NS));
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
        [prevEfficacy_with2first_tmp,...
                    prevEfficacy_pureNback_tmp,...
                    prevEfficacy_bis_with2first_tmp,...
                    prevEfficacy_bis_pureNback_tmp,...
                    prevEfficacy_ter_with2first_tmp,...
                    prevEfficacy_ter_pureNback_tmp] = extract_mental_previous_efficacy(subBehaviorFolder, sub_nm, run_nm, 'mental');
        [E_chosen_tmp] = extract_E_chosen(subBehaviorFolder, sub_nm, run_nm, 'mental');
        prev_mental_eff(run_trials_idx, iS) = prevEfficacy_ter_with2first_tmp;
        Ech(run_trials_idx, iS) = E_chosen_tmp;
    end % run loop
    
    %% perform correlation for each subject
    [r_corr(iS), ~, ~, prev_mental_eff_fit(:,iS)] = glm_package(Ech(:, iS), prev_mental_eff(:, iS), 'normal', 'on');
    
    %% extract bins
    for iE = E_levels
        jE = iE + 1;
        curr_Ech_idx = Ech(:,iS) == iE;
        if sum(curr_Ech_idx) > 0
            prev_mental_eff_bin(jE,iS) = mean(prev_mental_eff(curr_Ech_idx, iS),'omitnan');
            prev_mental_eff_fit_bin(jE,iS) = mean(prev_mental_eff_fit(curr_Ech_idx, iS),'omitnan');
        end
    end % effort chosen level loop
end % loop through subjects

%% test correlation
[~,pval] = ttest(r_corr);

%% average across subjects
[prev_mental_eff_bin_avg, prev_mental_eff_bin_sem] = mean_sem_sd(prev_mental_eff_bin, 2);
[prev_mental_eff_fit_bin_avg, prev_mental_eff_fit_bin_sem] = mean_sem_sd(prev_mental_eff_fit_bin, 2);

%% figure display
if fig_disp == 1
    fig;
    er_hdl = errorbar(E_levels, prev_mental_eff_bin_avg, prev_mental_eff_bin_sem);
    errorbar_hdl_upgrade(er_hdl);
    fit_hdl = plot(E_levels, prev_mental_eff_fit_bin_avg);
    fit_hdl_upgrade(fit_hdl);
    xticks(0:3);
    xticklabels({'0','1','2','3'});
    xlabel('E chosen');
    ylabel('Previous Mental Efficiency');
    xlim([-0.01 3.1]);
end % figure display

end % function