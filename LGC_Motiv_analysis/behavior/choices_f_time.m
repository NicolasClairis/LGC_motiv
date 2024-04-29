function[choice_f_trialN, pChoice_f_trialN,...
    trialN_f_trialN] = choices_f_time(study_nm, sub_nm,...
    condition, study_path, nBins, mdl_nm)
% [choice_f_trialN, pChoice_f_trialN,...
%     trialN_f_trialN] = choices_f_time(study_nm, sub_nm,...
%     condition, study_path, nBins, mdl_nm)
%choices_f_time will extract choice in function of time (ie trial number),
% and will also show how the model fits.
%
% INPUTS
% study_nm: study name
%
% sub_nm: subject name
%
% condition: condition
%
% study_path: study path
%
% nBins: number of bins
%
% mdl_nm: bayesian model to extract
%
% OUTPUTS
% choice_f_trialN: structure with choice in function of trial binned,
% split per run or averaging across runs
%
% pChoice_f_trialN: structure with p(choice) in function of trial binned,
% split per run or averaging across runs
%
% trialN_f_trialN: structure with trial number corresponding to each bin,
% split per run or averaging across runs

%% working directories
sub_CID_nm = ['CID',sub_nm];
subBehaviorFolder = [fullfile(study_path, sub_CID_nm, 'behavior'), filesep];
git_bayesian_data_folder = [fullfile('C:','Users','clairis','Desktop',...
    'GitHub','LGC_motiv','LGC_Motiv_results',...
    study_nm,'bayesian_modeling'), filesep];

%% define number of bins
if ~exist('nBins','var') || isempty(nBins)
    nBins = 6;
end

%% select relevant runs
[runs, n_runs] = runs_definition(study_nm, sub_nm, condition);

%% main parameters
task_for_bins = {'Ep','Em','EpEm'};
nTasks_for_bins = length(task_for_bins);
for iT = 1:nTasks_for_bins
    task_nm1 = task_for_bins{iT};
    switch task_nm1
        case {'Ep','Em'}
            [choice_f_trialN.(task_nm1).perRun.run1, pChoice_f_trialN.(task_nm1).perRun.run1,...
                trialN_f_trialN.(task_nm1).perRun.run1,...
                choice_f_trialN.(task_nm1).perRun.run2, pChoice_f_trialN.(task_nm1).perRun.run2,...
                trialN_f_trialN.(task_nm1).perRun.run2] = deal(NaN(nBins,1));
        case 'EpEm'
            [choice_f_trialN.(task_nm1).perRun.run1, pChoice_f_trialN.(task_nm1).perRun.run1,...
                trialN_f_trialN.(task_nm1).perRun.run1,...
                choice_f_trialN.(task_nm1).perRun.run2, pChoice_f_trialN.(task_nm1).perRun.run2,...
                trialN_f_trialN.(task_nm1).perRun.run2,...
                choice_f_trialN.(task_nm1).perRun.run3, pChoice_f_trialN.(task_nm1).perRun.run3,...
                trialN_f_trialN.(task_nm1).perRun.run3,...
                choice_f_trialN.(task_nm1).perRun.run4, pChoice_f_trialN.(task_nm1).perRun.run4,...
                trialN_f_trialN.(task_nm1).perRun.run4] = deal(NaN(nBins,1));
    end % task type
end
nTrials_perRun = 54;
trials = 1:nTrials_perRun;

for iR = 1:n_runs
    % extract corresponding run number
    jRun = runs.runsToKeep(iR);
    run_nm = num2str(jRun);
    EpEm_run_nm = ['run',run_nm];
    % run number for Ep/Em task (r1 or r2)
    switch jRun
        case {1,2}
            task_run = 1;
        case {3,4}
            task_run = 2;
    end
    task_nm2 = ['run',num2str(task_run)];
    % extract task type for the current run
    task_type_nm = runs.tasks{iR};
    switch task_type_nm
        case 'Ep'
            task_fullName = 'physical';
        case 'Em'
            task_fullName = 'mental';
    end
    
    %% extract choices
    [choices_tmp] = extract_choice_hE(subBehaviorFolder,...
        sub_nm, run_nm, task_fullName);

    %% load bayesian SV
    [~, ~, ~, pChoice_tmp] = extract_bayesian_mdl(git_bayesian_data_folder, subBehaviorFolder,...
        sub_nm, run_nm, task_fullName, mdl_nm);
    
    %% extract bins
    % extraction  per task
    [choice_f_trialN.(task_type_nm).perRun.(task_nm2)(:),...
        trialN_f_trialN.(task_type_nm).perRun.(task_nm2)(:)] = do_bin2(choices_tmp, trials, nBins, 0);
    [pChoice_f_trialN.(task_type_nm).perRun.(task_nm2)(:),...
        trialN_f_trialN.(task_type_nm).perRun.(task_nm2)(:)] = do_bin2(pChoice_tmp, trials, nBins, 0);
    % extraction global
    [choice_f_trialN.EpEm.perRun.(EpEm_run_nm)(:),...
        trialN_f_trialN.EpEm.perRun.(EpEm_run_nm)(:)] = do_bin2(choices_tmp, trials, nBins, 0);
    [pChoice_f_trialN.EpEm.perRun.(EpEm_run_nm)(:),...
        trialN_f_trialN.EpEm.perRun.(EpEm_run_nm)(:)] = do_bin2(pChoice_tmp, trials, nBins, 0);
end % run loop

%% average data across runs
for iT = 1:nTasks_for_bins
    task_nm3 = task_for_bins{iT};
    switch task_nm3
        case {'Ep','Em'}
            choice_f_trialN.(task_nm3).aRuns = mean([choice_f_trialN.(task_nm3).perRun.run1,...
                choice_f_trialN.(task_nm3).perRun.run2],2,'omitnan');
            pChoice_f_trialN.(task_nm3).aRuns = mean([pChoice_f_trialN.(task_nm3).perRun.run1,...
                pChoice_f_trialN.(task_nm3).perRun.run2],2,'omitnan');
            trialN_f_trialN.(task_nm3).aRuns = mean([trialN_f_trialN.(task_nm3).perRun.run1,...
                trialN_f_trialN.(task_nm3).perRun.run2],2,'omitnan');
        case 'EpEm'
            choice_f_trialN.EpEm.aRuns = mean([choice_f_trialN.EpEm.perRun.run1,...
                choice_f_trialN.EpEm.perRun.run2,...
                choice_f_trialN.EpEm.perRun.run3,...
                choice_f_trialN.EpEm.perRun.run4],2,'omitnan');
            pChoice_f_trialN.EpEm.aRuns = mean([pChoice_f_trialN.EpEm.perRun.run1,...
                pChoice_f_trialN.EpEm.perRun.run2,...
                pChoice_f_trialN.EpEm.perRun.run3,...
                pChoice_f_trialN.EpEm.perRun.run4],2,'omitnan');
            trialN_f_trialN.EpEm.aRuns = mean([trialN_f_trialN.EpEm.perRun.run1,...
                trialN_f_trialN.EpEm.perRun.run2,...
                trialN_f_trialN.EpEm.perRun.run3,...
                trialN_f_trialN.EpEm.perRun.run4],2,'omitnan');
    end % task type
end % loop over tasks
end % function