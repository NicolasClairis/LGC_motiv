function[r_corr, pval] = Ech_f_prev_mental_efficiency(fig_disp)
% [r_corr, pval] = Ech_f_prev_mental_efficiency(fig_disp)
% Ech_f_prev_mental_efficiency tests whether the mental efficiency (defined as
% the number of correct answers in the mental task divided by the total
% time taken to finish the performance in the current trial (including the
% two first useless answers) in the previous trial predicts the level of effort chosen across
% trials.
%
% INPUTS
% fig_disp: display figure (1) or not (0)
%
% OUTPUTS
% r_corr: correlation coefficient for each subject for the effort chosen in
% function of the previous mental efficiency
%
% pval: p.value for the t.test looking at whether the effort chosen varies
% with the previous mental efficiency

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
[prev_mental_eff, Ech, Ech_fit] = deal(NaN(nTotalTrials,NS));
nBins = 6;
[prev_mental_eff_bin, Ech_bin, Ech_fit_bin] = deal(NaN(nBins, NS));
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
    [r_corr(iS), ~, ~, Ech_fit(:,iS)] = glm_package(prev_mental_eff(:, iS), Ech(:, iS), 'normal', 'on');
    
    %% extract bins
    [Ech_bin(:,iS), prev_mental_eff_bin(:,iS)] = do_bin2(Ech(:, iS), prev_mental_eff(:, iS), nBins, 0);
    [Ech_fit_bin(:,iS)] = do_bin2(Ech_fit(:, iS), prev_mental_eff(:, iS), nBins, 0);
end % loop through subjects

%% test correlation
[~,pval] = ttest(r_corr);

%% average across subjects
prev_mental_eff_bin_avg = mean_sem_sd(prev_mental_eff_bin, 2);
Ech_fit_bin_avg = mean_sem_sd(Ech_fit_bin, 2);
[Ech_bin_avg, Ech_bin_sem] = mean_sem_sd(Ech_bin, 2);

%% figure display
if fig_disp == 1
    fig;
    er_hdl = errorbar(prev_mental_eff_bin_avg,...
        Ech_bin_avg, Ech_bin_sem);
    errorbar_hdl_upgrade(er_hdl);
    fit_hdl = plot(prev_mental_eff_bin_avg, Ech_fit_bin_avg);
    fit_hdl_upgrade(fit_hdl);
    xlabel('Previous Mental Efficiency');
%     yticks(0:3);
%     yticklabels({'0','1','2','3'});
    ylabel('E chosen');
end % figure display

end % function