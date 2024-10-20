function[choice_f_pChoice, pChoice_f_pChoice,...
    choice_f_deltaNV, deltaNV_f_deltaNV,...
    pChoice_f_deltaNV] = choices_f_SV(study_nm, sub_nm,...
    condition, study_path, nBins, mdl_nm)
% [choice_f_pChoice, pChoice_f_pChoice,...
%     choice_f_deltaNV, deltaNV_f_deltaNV,...
%     pChoice_f_deltaNV] = choices_f_SV(study_nm, sub_nm,...
%     condition, study_path, nBins, mdl_nm)
%choices_f_SV will extract choice in function of SV, either with the delta
%between the subjective value of the two options or based on the
%probability of choosing (ie after the softmax transformation).
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
% choice_f_pChoice: structure with choice in function of p(choice) binned,
% split per run or averaging across runs
%
% pChoice_f_pChoice: structure with p(choice) in function of p(choice) binned,
% split per run or averaging across runs
%
% choice_f_deltaNV: structure with choice in function of deltaNV binned,
% split per run or averaging across runs
%
% deltaNV_f_deltaNV: structure with deltaNV in function of deltaNV binned,
% split per run or averaging across runs
% 
% pChoice_f_deltaNV: structure with p(choice) in function of deltaNV
% binned, split per run or averaging across runs

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
            [choice_f_pChoice.(task_nm1).perRun.run1, pChoice_f_pChoice.(task_nm1).perRun.run1,...
                choice_f_deltaNV.(task_nm1).perRun.run1, deltaNV_f_deltaNV.(task_nm1).perRun.run1,...
                pChoice_f_deltaNV.(task_nm1).perRun.run1,...
                choice_f_pChoice.(task_nm1).perRun.run2, pChoice_f_pChoice.(task_nm1).perRun.run2,...
                choice_f_deltaNV.(task_nm1).perRun.run2, deltaNV_f_deltaNV.(task_nm1).perRun.run2,...
                pChoice_f_deltaNV.(task_nm1).perRun.run2] = deal(NaN(nBins,1));
        case 'EpEm'
            [choice_f_pChoice.(task_nm1).perRun.run1, pChoice_f_pChoice.(task_nm1).perRun.run1,...
                choice_f_deltaNV.(task_nm1).perRun.run1, deltaNV_f_deltaNV.(task_nm1).perRun.run1,...
                pChoice_f_deltaNV.(task_nm1).perRun.run1,...
                choice_f_pChoice.(task_nm1).perRun.run2, pChoice_f_pChoice.(task_nm1).perRun.run2,...
                choice_f_deltaNV.(task_nm1).perRun.run2, deltaNV_f_deltaNV.(task_nm1).perRun.run2,...
                pChoice_f_deltaNV.(task_nm1).perRun.run2,...
                choice_f_pChoice.(task_nm1).perRun.run3, pChoice_f_pChoice.(task_nm1).perRun.run3,...
                choice_f_deltaNV.(task_nm1).perRun.run3, deltaNV_f_deltaNV.(task_nm1).perRun.run3,...
                pChoice_f_deltaNV.(task_nm1).perRun.run3,...
                choice_f_pChoice.(task_nm1).perRun.run4, pChoice_f_pChoice.(task_nm1).perRun.run4,...
                choice_f_deltaNV.(task_nm1).perRun.run4, deltaNV_f_deltaNV.(task_nm1).perRun.run4,...
                pChoice_f_deltaNV.(task_nm1).perRun.run4] = deal(NaN(nBins,1));
    end % task type
end

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
    [~, ~, ~, pChoice_tmp,...
        deltaNV_hE_min_lE_plus_bias] = extract_bayesian_mdl(git_bayesian_data_folder, subBehaviorFolder,...
        sub_nm, run_nm, task_fullName, mdl_nm);
    
    %% extract bins
    [choice_f_pChoice.(task_type_nm).perRun.(task_nm2)(:),...
        pChoice_f_pChoice.(task_type_nm).perRun.(task_nm2)(:)] = do_bin2(choices_tmp, pChoice_tmp, nBins, 0);
    [choice_f_pChoice.EpEm.perRun.(EpEm_run_nm)(:),...
        pChoice_f_pChoice.EpEm.perRun.(EpEm_run_nm)(:)] = do_bin2(choices_tmp, pChoice_tmp, nBins, 0);
    [choice_f_deltaNV.(task_type_nm).perRun.(task_nm2)(:),...
        deltaNV_f_deltaNV.(task_type_nm).perRun.(task_nm2)(:)] = do_bin2(choices_tmp, deltaNV_hE_min_lE_plus_bias, nBins, 0);
    [choice_f_deltaNV.EpEm.perRun.(EpEm_run_nm)(:),...
        deltaNV_f_deltaNV.EpEm.perRun.(EpEm_run_nm)(:)] = do_bin2(choices_tmp, deltaNV_hE_min_lE_plus_bias, nBins, 0);
    [pChoice_f_deltaNV.(task_type_nm).perRun.(task_nm2)(:)] = do_bin2(pChoice_tmp, deltaNV_hE_min_lE_plus_bias, nBins, 0);
    [pChoice_f_deltaNV.EpEm.perRun.(EpEm_run_nm)(:)] = do_bin2(pChoice_tmp, deltaNV_hE_min_lE_plus_bias, nBins, 0);
end % run loop

%% average data across runs
for iT = 1:nTasks_for_bins
    task_nm3 = task_for_bins{iT};
    switch task_nm3
        case {'Ep','Em'}
            choice_f_pChoice.(task_nm3).aRuns = mean([choice_f_pChoice.(task_nm3).perRun.run1,...
                choice_f_pChoice.(task_nm3).perRun.run2],2,'omitnan');
            pChoice_f_pChoice.(task_nm3).aRuns = mean([pChoice_f_pChoice.(task_nm3).perRun.run1,...
                pChoice_f_pChoice.(task_nm3).perRun.run2],2,'omitnan');
            
            choice_f_deltaNV.(task_nm3).aRuns = mean([choice_f_deltaNV.(task_nm3).perRun.run1,...
                choice_f_deltaNV.(task_nm3).perRun.run2],2,'omitnan');
            deltaNV_f_deltaNV.(task_nm3).aRuns = mean([deltaNV_f_deltaNV.(task_nm3).perRun.run1,...
                deltaNV_f_deltaNV.(task_nm3).perRun.run2],2,'omitnan');
            pChoice_f_deltaNV.(task_nm3).aRuns = mean([pChoice_f_deltaNV.(task_nm3).perRun.run1,...
                pChoice_f_deltaNV.(task_nm3).perRun.run2],2,'omitnan');
        case 'EpEm'
            choice_f_pChoice.EpEm.aRuns = mean([choice_f_pChoice.EpEm.perRun.run1,...
                choice_f_pChoice.EpEm.perRun.run2,...
                choice_f_pChoice.EpEm.perRun.run3,...
                choice_f_pChoice.EpEm.perRun.run4],2,'omitnan');
            pChoice_f_pChoice.EpEm.aRuns = mean([pChoice_f_pChoice.EpEm.perRun.run1,...
                pChoice_f_pChoice.EpEm.perRun.run2,...
                pChoice_f_pChoice.EpEm.perRun.run3,...
                pChoice_f_pChoice.EpEm.perRun.run4],2,'omitnan');
            choice_f_deltaNV.EpEm.aRuns = mean([choice_f_deltaNV.EpEm.perRun.run1,...
                choice_f_deltaNV.EpEm.perRun.run2,...
                choice_f_deltaNV.EpEm.perRun.run3,...
                choice_f_deltaNV.EpEm.perRun.run4],2,'omitnan');
            deltaNV_f_deltaNV.EpEm.aRuns = mean([deltaNV_f_deltaNV.EpEm.perRun.run1,...
                deltaNV_f_deltaNV.EpEm.perRun.run2,...
                deltaNV_f_deltaNV.EpEm.perRun.run3,...
                deltaNV_f_deltaNV.EpEm.perRun.run4],2,'omitnan');
            pChoice_f_deltaNV.EpEm.aRuns = mean([pChoice_f_deltaNV.EpEm.perRun.run1,...
                pChoice_f_deltaNV.EpEm.perRun.run2,...
                pChoice_f_deltaNV.EpEm.perRun.run3,...
                pChoice_f_deltaNV.EpEm.perRun.run4],2,'omitnan');
    end % task type
end % loop over tasks
end % function