% function[] = Conf_DT_tradeoff_basic_tests(conf_var, fig_disp)
%[] = Conf_DT_tradeoff_basic_tests(conf_var, fig_disp)
% Conf_DT_tradeoff_basic_tests will perform 4 basic tests required by Jean
% Daunizeau to verify if our data fits with the oMCD model of Confidence/RT
% tradeoff.
% 1) choice = f(dV) split depending on if low or high confidence trials
% 2) RT = f(dV) split depending on if low or high confidence trials
% 3) Confidence = f(dV) split between choices consistent with dV and those
% inconsistent with dV
% 4) Conf = f(RT) split into different dV bins
%
% Note: conf_var will determine if you perform the tests on the rated
% binary confidence (parallel/post-decision) or on the confidence predicted
% by the computational model (continuous between 0 and 1, pre-decisional)
%
% INPUTS
% fig_disp: display figure (1) or no (0)? (yes by default)
%
% conf_var = decide if you use the confidence ratings (0/1) which are done
% simultaneously with the choice (ie "post"-decision confidence) or the
% confidence inferred by the computational model (ie "pre"-decision
% confidence solely based on the dV between the two options)
% (0) use confidence ratings (~post-decision confidence)
% (1) use confidence inferred by the model (~pre-decision confidence)
%
% OUTPUTS
%


%% default parameters
% confidence variable (confidence ratings by default, ie post-decision
% confidence)
if ~exist('conf_var','var') || isempty(conf_var) 
    conf_var = 0;
end
% figure display (yes by default)
if ~exist('fig_disp','var') || isempty(fig_disp) 
    fig_disp = 1;
end

%% subject selection
[study_nm, condition, ~, subject_id, NS] = sub_id;

%% working dir
computerRoot = LGCM_root_paths;
switch computerRoot
    case [fullfile('E:'),filesep] % lab
        gitFolder = [fullfile('C:','Users','clairis','Desktop','GitHub',...
            'LGC_motiv','LGC_Motiv_results','study1','bayesian_modeling'),filesep];
    case [filesep,filesep,fullfile('sv-nas1.rcp.epfl.ch','Sandi-lab','human_data_private','raw_data_subject'),filesep] % home
        gitFolder = [fullfile('C:','Users','Nicolas Clairis','Documents','GitHub',...
            'LGC_motiv','LGC_Motiv_results','study1','bayesian_modeling'),filesep];
end

%% general parameters
mdl_nm = 'mdl_6';
nTrialsPerRun = 54;
n_totalRuns = 4;
n_runsPerTask = 2;
n_bins = 6;
n_dV_bins = 3;
task_names = {'EpEm','Ep','Em'}; nTasks = length(task_names);
for iT = 1:nTasks
    task_nm = task_names{iT};
    switch task_nm
        case 'EpEm'
            n_runs0 = n_totalRuns;
        otherwise
            n_runs0 = n_runsPerTask;
    end
    [choice_f_dV.low_Conf.(task_nm).perSub.perRun, choice_f_dV.high_Conf.(task_nm).perSub.perRun,...
        pChoice_f_dV.low_Conf.(task_nm).perSub.perRun, pChoice_f_dV.high_Conf.(task_nm).perSub.perRun,...
        DT_f_dV.low_Conf.(task_nm).perSub.perRun, DT_f_dV.high_Conf.(task_nm).perSub.perRun,...
        dV_f_dV.low_Conf.(task_nm).perSub.perRun, dV_f_dV.high_Conf.(task_nm).perSub.perRun,...
        Conf_f_dV.consistent.(task_nm).perSub.perRun, Conf_f_dV.inconsistent.(task_nm).perSub.perRun,...
        dV_f_dV.consistent.(task_nm).perSub.perRun, dV_f_dV.inconsistent.(task_nm).perSub.perRun] = deal(NaN(n_bins, NS, n_runs0));
    [choice_f_dV.low_Conf.(task_nm).perSub.aRuns, choice_f_dV.high_Conf.(task_nm).perSub.aRuns,...
        pChoice_f_dV.low_Conf.(task_nm).perSub.aRuns, pChoice_f_dV.high_Conf.(task_nm).perSub.aRuns,...
        DT_f_dV.low_Conf.(task_nm).perSub.aRuns, DT_f_dV.high_Conf.(task_nm).perSub.aRuns,...
        dV_f_dV.low_Conf.(task_nm).perSub.aRuns, dV_f_dV.high_Conf.(task_nm).perSub.aRuns,...
        Conf_f_dV.consistent.(task_nm).perSub.aRuns, Conf_f_dV.inconsistent.(task_nm).perSub.aRuns,...
        dV_f_dV.consistent.(task_nm).perSub.aRuns, dV_f_dV.inconsistent.(task_nm).perSub.aRuns] = deal(NaN(n_bins, NS));
    
    dV_names = cell(1,n_dV_bins);
    for i_dV_bins = 1:n_dV_bins
        dV_nm = ['dV',num2str(i_dV_bins)];
        dV_names{i_dV_bins} = dV_nm;
        [Conf_f_DT.(task_nm).perSub.perRun.(dV_nm),...
            DT_f_DT.(task_nm).perSub.perRun.(dV_nm)] = deal((NaN(n_bins, NS, n_runs0)));
        [Conf_f_DT.(task_nm).perSub.aRuns.(dV_nm),...
            DT_f_DT.(task_nm).perSub.aRuns.(dV_nm)] = deal((NaN(n_bins, NS)));
    end % dV loop
end % task loop

%% loop over subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    subBehaviorFolder = [fullfile(computerRoot,study_nm,['CID',sub_nm],'behavior'),filesep];
    [runs, n_runs] = runs_definition(study_nm, sub_nm, condition);
    
    for iR = 1:n_runs
        jR = runs.runsToKeep(iR);
        run_nm = num2str(jR);
        run_task_nm = runs.tasks{iR};
        [task_fullName] = task_fullName_extraction(run_task_nm);
        switch jR
            case {1,2}
                kRun = 1;
            case {3,4}
                kRun = 2;
        end
        
        % extract actual choice
        [choice_hE_tmp] = extract_choice_hE(subBehaviorFolder,...
            sub_nm, run_nm, task_fullName);

        % extract predicted confidence + subjective value
        [~, deltaNV_hE_min_lE_tmp, confidence_highE_tmp, pChoice_tmp,...
            deltaNV_hE_min_lE_plus_bias_tmp] = extract_bayesian_mdl(gitFolder, subBehaviorFolder,...
            sub_nm, run_nm, task_fullName, mdl_nm);
        dV_tmp = deltaNV_hE_min_lE_plus_bias_tmp;
        
        % extract RT
        [DT_tmp] = extract_RT(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        
        % extract confidence variable
        switch conf_var
            case 0 % confidence ratings
                [conf_rtg_tmp] = extract_confidence_rating(subBehaviorFolder, sub_nm, run_nm, task_fullName);
                conf_tmp = conf_rtg_tmp;
            case 1 % confidence inferred by the model
                conf_tmp = confidence_highE_tmp;
        end
        
        %% extract the data following the rule of the test
        %% 1) choice = f(dV) split depending on if low or high confidence trials
        % identify low vs high confidence trials
        switch conf_var
            case 0 % binary ratings
                hConf_trials = conf_tmp == 1;
                lConf_trials = conf_tmp == 0;
            case 1 % continuous prediction from the model
                hConf_trials = conf_tmp > 0.5;
                lConf_trials = conf_tmp <= 0.5;
        end
        % pool across physical and mental tasks
        % extract data for low confidence trials
        [choice_f_dV.low_Conf.EpEm.perSub.perRun(:,iS,jR),...
            dV_f_dV.low_Conf.EpEm.perSub.perRun(:,iS,jR)] = do_bin2(choice_hE_tmp(lConf_trials), dV_tmp(lConf_trials), n_bins, 0);
        [pChoice_f_dV.low_Conf.EpEm.perSub.perRun(:,iS,jR)] = do_bin2(pChoice_tmp(lConf_trials), dV_tmp(lConf_trials), n_bins, 0);
        % extract data for high confidence trials
        [choice_f_dV.high_Conf.EpEm.perSub.perRun(:,iS,jR),...
            dV_f_dV.high_Conf.EpEm.perSub.perRun(:,iS,jR)] = do_bin2(choice_hE_tmp(hConf_trials), dV_tmp(hConf_trials), n_bins, 0);
        [pChoice_f_dV.high_Conf.EpEm.perSub.perRun(:,iS,jR)] = do_bin2(pChoice_tmp(hConf_trials), dV_tmp(hConf_trials), n_bins, 0);
        
        % split per task
        % extract data for low confidence trials
        [choice_f_dV.low_Conf.(run_task_nm).perSub.perRun(:,iS,kRun),...
            dV_f_dV.low_Conf.(run_task_nm).perSub.perRun(:,iS,kRun)] = do_bin2(choice_hE_tmp(lConf_trials), dV_tmp(lConf_trials), n_bins, 0);
        [pChoice_f_dV.low_Conf.(run_task_nm).perSub.perRun(:,iS,kRun)] = do_bin2(pChoice_tmp(lConf_trials), dV_tmp(lConf_trials), n_bins, 0);
        % extract data for high confidence trials
        [choice_f_dV.high_Conf.(run_task_nm).perSub.perRun(:,iS,kRun),...
            dV_f_dV.high_Conf.(run_task_nm).perSub.perRun(:,iS,kRun)] = do_bin2(choice_hE_tmp(hConf_trials), dV_tmp(hConf_trials), n_bins, 0);
        [pChoice_f_dV.high_Conf.(run_task_nm).perSub.perRun(:,iS,kRun)] = do_bin2(pChoice_tmp(hConf_trials), dV_tmp(hConf_trials), n_bins, 0);
        
        %% 2) RT = f(dV) split depending on if low or high confidence trials
        % pool across physical and mental tasks
        % extract data for low confidence trials
        [DT_f_dV.low_Conf.EpEm.perSub.perRun(:,iS,jR)] = do_bin2(DT_tmp(lConf_trials), dV_tmp(lConf_trials), n_bins, 0);
        % extract data for high confidence trials
        [DT_f_dV.high_Conf.EpEm.perSub.perRun(:,iS,jR)] = do_bin2(DT_tmp(hConf_trials), dV_tmp(hConf_trials), n_bins, 0);
        
        % split per task
        % extract data for low confidence trials
        [DT_f_dV.low_Conf.(run_task_nm).perSub.perRun(:,iS,kRun)] = do_bin2(DT_tmp(lConf_trials), dV_tmp(lConf_trials), n_bins, 0);
        % extract data for high confidence trials
        [DT_f_dV.high_Conf.(run_task_nm).perSub.perRun(:,iS,kRun)] = do_bin2(DT_tmp(hConf_trials), dV_tmp(hConf_trials), n_bins, 0);
        
        %% 3) Confidence = f(dV) split between choices consistent with dV and those
        % inconsistent with dV
        % identify trials consistent vs inconsistent with dV
        consistent_trials = ((choice_hE_tmp' == 1) & (dV_tmp >= 0)) | ((choice_hE_tmp' == 0) & (dV_tmp < 0)); % consistent trials = choice high Effort while dV >= 0 or choice low Effort while dV < 0
        inconsistent_trials = ((choice_hE_tmp' == 0) & (dV_tmp >= 0)) | ((choice_hE_tmp' == 1) & (dV_tmp < 0)); % inconsistent trials = choice high Effort while dV < 0 or choice low Effort while dV >= 0
        % pool across physical and mental tasks
        [Conf_f_dV.consistent.EpEm.perSub.perRun(:,iS,jR),...
            dV_f_dV.consistent.EpEm.perSub.perRun(:,iS,jR)] = do_bin2(conf_tmp(consistent_trials), dV_tmp(consistent_trials), n_bins, 0);
        [Conf_f_dV.inconsistent.EpEm.perSub.perRun(:,iS,jR),...
            dV_f_dV.inconsistent.EpEm.perSub.perRun(:,iS,jR)] = do_bin2(conf_tmp(inconsistent_trials), dV_tmp(inconsistent_trials), n_bins, 0);
        
        % split per task
        [Conf_f_dV.consistent.(run_task_nm).perSub.perRun(:,iS,kRun),...
            dV_f_dV.consistent.(run_task_nm).perSub.perRun(:,iS,kRun)] = do_bin2(conf_tmp(consistent_trials), dV_tmp(consistent_trials), n_bins, 0);
        [Conf_f_dV.inconsistent.(run_task_nm).perSub.perRun(:,iS,kRun),...
            dV_f_dV.inconsistent.(run_task_nm).perSub.perRun(:,iS,kRun)] = do_bin2(conf_tmp(inconsistent_trials), dV_tmp(inconsistent_trials), n_bins, 0);
        
        %% 4) Conf = f(RT) split into different dV bins
        % identify trial index for each bin
        [~, dV_bin_dV_tmp, dV_bin_trials_tmp] = do_bin2(1:nTrialsPerRun, dV_tmp, n_dV_bins,0);
        for i_dV_bins = 1:n_dV_bins
            dV_nm = ['dV',num2str(i_dV_bins)];
            dV_trials = dV_bin_trials_tmp == i_dV_bins;
            
            % pool across physical and mental tasks
            [Conf_f_DT.EpEm.perSub.perRun.(dV_nm)(:,iS,jR),...
                DT_f_DT.EpEm.perSub.perRun.(dV_nm)(:,iS,jR)] = do_bin2(conf_tmp(dV_trials), DT_tmp(dV_trials), n_bins, 0);
            
            % split per task
            [Conf_f_DT.(run_task_nm).perSub.perRun.(dV_nm)(:,iS,kRun),...
                DT_f_DT.(run_task_nm).perSub.perRun.(dV_nm)(:,iS,kRun)] = do_bin2(conf_tmp(dV_trials), DT_tmp(dV_trials), n_bins, 0);
        end % dV bin loop
    end % run loop
    
    %% average data per subject (per task and across tasks)
    for iT = 1:nTasks
        task_nm = task_names{iT};
        % 1) choice = f(dV)
        % low conf
        choice_f_dV.low_Conf.(task_nm).perSub.aRuns(:,iS) = mean(choice_f_dV.low_Conf.(task_nm).perSub.perRun(:,iS,:), 3,'omitnan');
        dV_f_dV.low_Conf.(task_nm).perSub.aRuns(:,iS) = mean(dV_f_dV.low_Conf.(task_nm).perSub.perRun(:,iS,:), 3,'omitnan');
        pChoice_f_dV.low_Conf.(task_nm).perSub.aRuns(:,iS) = mean(pChoice_f_dV.low_Conf.(task_nm).perSub.perRun(:,iS,:), 3,'omitnan');
        % high conf
        choice_f_dV.high_Conf.(task_nm).perSub.aRuns(:,iS) = mean(choice_f_dV.high_Conf.(task_nm).perSub.perRun(:,iS,:), 3,'omitnan');
        dV_f_dV.high_Conf.(task_nm).perSub.aRuns(:,iS) = mean(dV_f_dV.high_Conf.(task_nm).perSub.perRun(:,iS,:), 3,'omitnan');
        pChoice_f_dV.high_Conf.(task_nm).perSub.aRuns(:,iS) = mean(pChoice_f_dV.high_Conf.(task_nm).perSub.perRun(:,iS,:), 3,'omitnan');
        
        % 2) RT = f(dV)
        % low conf
        DT_f_dV.low_Conf.(task_nm).perSub.aRuns(:,iS) = mean(DT_f_dV.low_Conf.(task_nm).perSub.perRun(:,iS,:),3,'omitnan');
        % high conf
        DT_f_dV.high_Conf.(task_nm).perSub.aRuns(:,iS) = mean(DT_f_dV.high_Conf.(task_nm).perSub.perRun(:,iS,:),3,'omitnan');
        
        % 3) Confidence = f(dV)
        % consistent trials
        Conf_f_dV.consistent.(task_nm).perSub.aRuns(:,iS) = mean(Conf_f_dV.consistent.(task_nm).perSub.perRun(:,iS,:),3,'omitnan');
        dV_f_dV.consistent.(task_nm).perSub.aRuns(:,iS) = mean(dV_f_dV.consistent.(task_nm).perSub.perRun(:,iS,:),3,'omitnan');
        % inconsistent trials
        Conf_f_dV.inconsistent.(task_nm).perSub.aRuns(:,iS) = mean(Conf_f_dV.inconsistent.(task_nm).perSub.perRun(:,iS,:),3,'omitnan');
        dV_f_dV.inconsistent.(task_nm).perSub.aRuns(:,iS) = mean(dV_f_dV.inconsistent.(task_nm).perSub.perRun(:,iS,:),3,'omitnan');
        
        % 4) Conf = f(RT)
        for i_dV_bins = 1:n_dV_bins
            dV_nm = ['dV',num2str(i_dV_bins)];
            Conf_f_DT.(task_nm).perSub.aRuns.(dV_nm)(:,iS) = mean(Conf_f_DT.(task_nm).perSub.perRun.(dV_nm)(:,iS,:),3,'omitnan');
            DT_f_DT.(task_nm).perSub.aRuns.(dV_nm)(:,iS) = mean(DT_f_DT.(task_nm).perSub.perRun.(dV_nm)(:,iS,:),3,'omitnan');
        end % dV bin
    end % task loop
    
end % subject loop

%% average + SEM/SD across subjects
for iT = 1:nTasks
    task_nm = task_names{iT};
    
    % 1) choice = f(dV) split high/low Conf
    % low conf
    [choice_f_dV.low_Conf.(task_nm).aSubs.mean,...
        choice_f_dV.low_Conf.(task_nm).aSubs.sem,...
        choice_f_dV.low_Conf.(task_nm).aSubs.sd] = mean_sem_sd(choice_f_dV.low_Conf.(task_nm).perSub.aRuns, 2);
    [dV_f_dV.low_Conf.(task_nm).aSubs.mean,...
        dV_f_dV.low_Conf.(task_nm).aSubs.sem,...
        dV_f_dV.low_Conf.(task_nm).aSubs.sd] = mean_sem_sd(dV_f_dV.low_Conf.(task_nm).perSub.aRuns, 2);
    [pChoice_f_dV.low_Conf.(task_nm).aSubs.mean,...
        pChoice_f_dV.low_Conf.(task_nm).aSubs.sem,...
        pChoice_f_dV.low_Conf.(task_nm).aSubs.sd] = mean_sem_sd(pChoice_f_dV.low_Conf.(task_nm).perSub.aRuns, 2);
    % high conf
    [choice_f_dV.high_Conf.(task_nm).aSubs.mean,...
        choice_f_dV.high_Conf.(task_nm).aSubs.sem,...
        choice_f_dV.high_Conf.(task_nm).aSubs.sd] = mean_sem_sd(choice_f_dV.high_Conf.(task_nm).perSub.aRuns, 2);
    [dV_f_dV.high_Conf.(task_nm).aSubs.mean,...
        dV_f_dV.high_Conf.(task_nm).aSubs.sem,...
        dV_f_dV.high_Conf.(task_nm).aSubs.sd] = mean_sem_sd(dV_f_dV.high_Conf.(task_nm).perSub.aRuns, 2);
    [pChoice_f_dV.high_Conf.(task_nm).aSubs.mean,...
        pChoice_f_dV.high_Conf.(task_nm).aSubs.sem,...
        pChoice_f_dV.high_Conf.(task_nm).aSubs.sd] = mean_sem_sd(pChoice_f_dV.high_Conf.(task_nm).perSub.aRuns, 2);
    
    % 2) RT = f(dV) split high/low Conf
    % low conf
    [DT_f_dV.low_Conf.(task_nm).aSubs.mean,...
        DT_f_dV.low_Conf.(task_nm).aSubs.sem,...
        DT_f_dV.low_Conf.(task_nm).aSubs.sd] = mean_sem_sd(DT_f_dV.low_Conf.(task_nm).perSub.aRuns, 2);
    % high conf
    [DT_f_dV.high_Conf.(task_nm).aSubs.mean,...
        DT_f_dV.high_Conf.(task_nm).aSubs.sem,...
        DT_f_dV.high_Conf.(task_nm).aSubs.sd] = mean_sem_sd(DT_f_dV.high_Conf.(task_nm).perSub.aRuns, 2);
    
    % 3) Confidence = f(dV) split consistent/inconsistent
    [Conf_f_dV.consistent.(task_nm).aSubs.mean,...
        Conf_f_dV.consistent.(task_nm).aSubs.sem,...
        Conf_f_dV.consistent.(task_nm).aSubs.sd] = mean_sem_sd(Conf_f_dV.consistent.(task_nm).perSub.aRuns,2);
    [dV_f_dV.consistent.(task_nm).aSubs.mean,...
        dV_f_dV.consistent.(task_nm).aSubs.sem,...
        dV_f_dV.consistent.(task_nm).aSubs.sd] = mean_sem_sd(dV_f_dV.consistent.(task_nm).perSub.aRuns,2);
    [Conf_f_dV.inconsistent.(task_nm).aSubs.mean,...
        Conf_f_dV.inconsistent.(task_nm).aSubs.sem,...
        Conf_f_dV.inconsistent.(task_nm).aSubs.sd] = mean_sem_sd(Conf_f_dV.inconsistent.(task_nm).perSub.aRuns,2);
    [dV_f_dV.inconsistent.(task_nm).aSubs.mean,...
        dV_f_dV.inconsistent.(task_nm).aSubs.sem,...
        dV_f_dV.inconsistent.(task_nm).aSubs.sd] = mean_sem_sd(dV_f_dV.inconsistent.(task_nm).perSub.aRuns,2);
    
    % 4) Confidence = f(RT) split into different dV bins
    for i_dV_bins = 1:n_dV_bins
        dV_nm = ['dV',num2str(i_dV_bins)];
        [Conf_f_DT.(task_nm).aSubs.mean.(dV_nm),...
            Conf_f_DT.(task_nm).aSubs.sem.(dV_nm),...
            Conf_f_DT.(task_nm).aSubs.sd.(dV_nm)] = mean_sem_sd(Conf_f_DT.(task_nm).perSub.aRuns.(dV_nm),2);
        [DT_f_DT.(task_nm).aSubs.mean.(dV_nm),...
            DT_f_DT.(task_nm).aSubs.sem.(dV_nm),...
            DT_f_DT.(task_nm).aSubs.sd.(dV_nm)] = mean_sem_sd(DT_f_DT.(task_nm).perSub.aRuns.(dV_nm),2);
    end
end % task loop
%% figures
if fig_disp == 1
    [pSize, lWidth, col, mSize, fontName] = general_fig_prm;
    
    fig1_choice_f_dV = fig;
    fig2_RT_f_dV = fig;
    fig3_Conf_f_dV = fig;
    fig4_Conf_f_RT = fig;
    
    for iT = 1:nTasks
        task_nm = task_names{iT};
        switch task_nm
            case 'EpEm'
                jT = 1.5;
            otherwise
                jT = 1 + iT;
        end
        
        %% 1) choice = f(dV) split high/low Conf
        figure(fig1_choice_f_dV);
        subplot(2,2,jT); hold on;
        ch_lC_hdl = errorbar(dV_f_dV.low_Conf.(task_nm).aSubs.mean,...
            choice_f_dV.low_Conf.(task_nm).aSubs.mean,...
            choice_f_dV.low_Conf.(task_nm).aSubs.sem,...
            choice_f_dV.low_Conf.(task_nm).aSubs.sem,...
            dV_f_dV.low_Conf.(task_nm).aSubs.sem,...
            dV_f_dV.low_Conf.(task_nm).aSubs.sem);
        ch_hC_hdl = errorbar(dV_f_dV.high_Conf.(task_nm).aSubs.mean,...
            choice_f_dV.high_Conf.(task_nm).aSubs.mean,...
            choice_f_dV.high_Conf.(task_nm).aSubs.sem,...
            choice_f_dV.high_Conf.(task_nm).aSubs.sem,...
            dV_f_dV.high_Conf.(task_nm).aSubs.sem,...
            dV_f_dV.high_Conf.(task_nm).aSubs.sem);
        ch_lC_hdl.LineStyle = 'none';
        ch_hC_hdl.LineStyle = 'none';
        ch_lC_hdl.Color = 'r';
        ch_hC_hdl.Color = 'b';
        ch_lC_hdl.LineWidth = lWidth;
        ch_hC_hdl.LineWidth = lWidth;
        % overlay p(choice)
        pCh_lC_hdl = plot(dV_f_dV.low_Conf.(task_nm).aSubs.mean,...
            pChoice_f_dV.low_Conf.(task_nm).aSubs.mean);
        pCh_hC_hdl = plot(dV_f_dV.high_Conf.(task_nm).aSubs.mean,...
            pChoice_f_dV.high_Conf.(task_nm).aSubs.mean);
        pCh_lC_hdl.LineStyle = '--';
        pCh_hC_hdl.LineStyle = '--';
        pCh_lC_hdl.Color = 'r';
        pCh_hC_hdl.Color = 'b';
        pCh_lC_hdl.LineWidth = lWidth;
        pCh_hC_hdl.LineWidth = lWidth;
        xlabel('dV');
        ylabel('Choice (%)');
        % add reference lines
        line([0 0],[0 1],'LineStyle','-','LineWidth',1,'Color','k');
        xlim_vals = xlim();
        line(xlim_vals,[0.5 0.5],'LineStyle','-','LineWidth',1,'Color','k');
        line(xlim_vals,[1 1],'LineStyle','-','LineWidth',1,'Color','k');
        legend([ch_hC_hdl, ch_lC_hdl],...
            {'high conf.','low conf.'});
        legend('boxoff');
        
        %% 2) RT = f(dV) split high/low Conf
        figure(fig2_RT_f_dV);
        subplot(2,2,jT); hold on;
        DT_lC_hdl = errorbar(dV_f_dV.low_Conf.(task_nm).aSubs.mean,...
            DT_f_dV.low_Conf.(task_nm).aSubs.mean,...
            DT_f_dV.low_Conf.(task_nm).aSubs.sem,...
            DT_f_dV.low_Conf.(task_nm).aSubs.sem,...
            dV_f_dV.low_Conf.(task_nm).aSubs.sem,...
            dV_f_dV.low_Conf.(task_nm).aSubs.sem);
        DT_hC_hdl = errorbar(dV_f_dV.high_Conf.(task_nm).aSubs.mean,...
            DT_f_dV.high_Conf.(task_nm).aSubs.mean,...
            DT_f_dV.high_Conf.(task_nm).aSubs.sem,...
            DT_f_dV.high_Conf.(task_nm).aSubs.sem,...
            dV_f_dV.high_Conf.(task_nm).aSubs.sem,...
            dV_f_dV.high_Conf.(task_nm).aSubs.sem);
        DT_lC_hdl.LineStyle = 'none';
        DT_hC_hdl.LineStyle = 'none';
        DT_lC_hdl.Color = 'r';
        DT_hC_hdl.Color = 'b';
        DT_lC_hdl.LineWidth = lWidth;
        DT_hC_hdl.LineWidth = lWidth;
        xlabel('dV');
        ylabel('DT (s)');
        % add reference lines
        ylim_vals = ylim();
        line([0 0],ylim_vals,'LineStyle','-','LineWidth',1,'Color','k');
        legend([DT_hC_hdl, DT_lC_hdl],...
            {'high conf.','low conf.'});
        legend('boxoff');
        
        %% 3) Confidence = f(dV) split consistent/inconsistent
        figure(fig3_Conf_f_dV);
        subplot(2,2,jT); hold on;
        Conf_consistent_hdl = errorbar(dV_f_dV.consistent.(task_nm).aSubs.mean,...
            Conf_f_dV.consistent.(task_nm).aSubs.mean,...
            Conf_f_dV.consistent.(task_nm).aSubs.sem,...
            Conf_f_dV.consistent.(task_nm).aSubs.sem,...
            dV_f_dV.consistent.(task_nm).aSubs.sem,...
            dV_f_dV.consistent.(task_nm).aSubs.sem);
        Conf_inconsistent_hdl = errorbar(dV_f_dV.inconsistent.(task_nm).aSubs.mean,...
            Conf_f_dV.inconsistent.(task_nm).aSubs.mean,...
            Conf_f_dV.inconsistent.(task_nm).aSubs.sem,...
            Conf_f_dV.inconsistent.(task_nm).aSubs.sem,...
            dV_f_dV.inconsistent.(task_nm).aSubs.sem,...
            dV_f_dV.inconsistent.(task_nm).aSubs.sem);
        Conf_consistent_hdl.LineStyle = 'none';
        Conf_inconsistent_hdl.LineStyle = 'none';
        Conf_consistent_hdl.Color = 'b';
        Conf_inconsistent_hdl.Color = 'r';
        Conf_consistent_hdl.LineWidth = lWidth;
        Conf_inconsistent_hdl.LineWidth = lWidth;
        xlabel('dV');
        ylabel('Confidence');
        % add reference lines
        ylim_vals = ylim();
        line([0 0],ylim_vals,'LineStyle','-','LineWidth',1,'Color','k');
        legend([Conf_consistent_hdl, Conf_inconsistent_hdl],...
            {'value-consistent','value-inconsistent'});
        legend('boxoff');
        
        %% 4) Confidence = f(RT) split into different dV bins
        figure(fig4_Conf_f_RT);
        subplot(2,2,jT); hold on;
        for i_dV_bins = 1:n_dV_bins
            dV_nm = ['dV',num2str(i_dV_bins)];
            dV_bin_hdl.(dV_nm) = errorbar(DT_f_DT.(task_nm).aSubs.mean.(dV_nm),...
                Conf_f_DT.(task_nm).aSubs.mean.(dV_nm),...
                Conf_f_DT.(task_nm).aSubs.sem.(dV_nm),...
                Conf_f_DT.(task_nm).aSubs.sem.(dV_nm),...
                DT_f_DT.(task_nm).aSubs.sem.(dV_nm),...
                DT_f_DT.(task_nm).aSubs.sem.(dV_nm));
            dV_bin_hdl.(dV_nm).LineStyle = 'none';
            switch i_dV_bins
                case 1
                    dV_bin_hdl.(dV_nm).Color = [255 255 178]./255;
                case 2
                    dV_bin_hdl.(dV_nm).Color = [254 204 92]./255;
                case 3
                    dV_bin_hdl.(dV_nm).Color = [253 141 60]./255;
                case 4
                    dV_bin_hdl.(dV_nm).Color = [240 59 32]./255;
                case 5
                    dV_bin_hdl.(dV_nm).Color = [189 0 38]./255;
                otherwise
                    dV_bin_hdl.(dV_nm).Color = [189 0 38]./255;
            end
            dV_bin_hdl.(dV_nm).LineWidth = lWidth;
        end % dV bin
        xlabel('RT (s)');
        ylabel('Confidence');
        switch n_dV_bins
            case 2
                legend([dV_bin_hdl.dV1, dV_bin_hdl.dV2],...
                    {'dV low','dV high'});
            case 3
                legend([dV_bin_hdl.dV1, dV_bin_hdl.dV2, dV_bin_hdl.dV3],...
                    {'dV low','dV middle','dV high'});
            otherwise
                error(['Please adapt case where n_dV_bins = ',num2str(n_dV_bins),' for legend']);
        end
        legend('boxoff');
    end % task loop
end % figure display

% end % function