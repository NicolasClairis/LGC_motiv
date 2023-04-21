function[b_forceLatency_perEch_f_time, pval, b_forceLatency_allE_f_time] = forceLatency_f_time(study_nm, subject_id, condition, figDisp)
%[b_forceLatency_perEch_f_time, pval, b_forceLatency_allE_f_time] = forceLatency_f_time(study_nm, subject_id, condition, figDisp)
% forceLatency_f_time will look, for each effort level, how the peak of force
% decreases over time
%
% INPUTS
% study_nm: study name 'study1'/'study2'
%
% subject_id: list of subjects
%
% condition: condition for subjects and runs to include
%
% figDisp: display figure (1) or not (0)?
%
% OUTPUTS
% b_forceLatency_perEch_f_time: nb effort level chosen*NS matrix with slope for
% each effort level chosen per subject
%
% pval: structure with p.value for the slope for each level of effort
% chosen
%
% b_forceLatency_allE_f_time: 1*NS matrix with slope across all efforts per subject


%% subject selection
if ~exist('study_nm','var') || isempty(study_nm) ||...
        ~ismember(study_nm,{'study1','study2'})
    study_nm = 'study1';
end
if ~exist('condition','var') || isempty(condition)
    condition = subject_condition;
end
if ~exist('subject_id','var') || isempty(subject_id)
    [subject_id, NS] = LGCM_subject_selection(study_nm, condition);
else
    NS = length(subject_id);
end
if ~exist('figDisp','var') || isempty(figDisp)
    figDisp = 1;
end

%% working directories
computerRoot = LGCM_root_paths;
dataRoot = [computerRoot, filesep, study_nm, filesep];

%% extract force peak
Ech_levels = 0:3;
n_Ech_lvl = length(Ech_levels);
b_forceLatency_perEch_f_time = NaN(n_Ech_lvl, NS);
b_forceLatency_allE_f_time = NaN(1,NS);
task_fullName = 'physical';
nTrialsPerRun = 54;
trialN = 1:nTrialsPerRun;
n_trialBins = 6;
for iEch = Ech_levels
    Ech_nm = ['Ech',num2str(iEch)];
    [forcePeak_f_time.(Ech_nm),...
        trialN_f_time.(Ech_nm),...
        forcePeak_fit_f_time.(Ech_nm)] = deal(NaN(n_trialBins, NS));
end
[forcePeak_f_time.allE,...
        trialN_f_time.allE,...
        forcePeak_fit_f_time.allE] = deal(NaN(n_trialBins, NS));
for iS = 1:NS
    sub_nm = subject_id{iS};
    subBehaviorFolder = [dataRoot, filesep, 'CID',sub_nm, filesep, 'behavior',filesep];
    runs = runs_definition(study_nm, sub_nm, condition);
    b_forceLatency_Sub = NaN(2,n_Ech_lvl,2);
    b_forceLatency_allE_Sub = NaN(2,2);
    [forceLatency_tmp.run1, trialN_tmp.run1, latency_fit_tmp.run1, ...
        forceLatency_tmp.run2, trialN_tmp.run2, latency_fit_tmp.run2] = deal(NaN(1,n_trialBins,n_Ech_lvl));
    [forceLatency_allE_tmp.run1, trialN_allE_tmp.run1, forceLatency_fit_allE_tmp.run1, ...
        forceLatency_allE_tmp.run2, trialN_allE_tmp.run2, forceLatency_fit_allE_tmp.run2] = deal(NaN(1,n_trialBins));
    for iRun = 1:runs.nb_runs.Ep
        jRun = runs.Ep.runsToKeep(iRun);
        run_nm = num2str(jRun);
        run_nm_bis = ['run',num2str(iRun)];
        latency_tmp = extract_grip_force(subBehaviorFolder, sub_nm, run_nm);
        [E_chosen_tmp] = extract_E_chosen(subBehaviorFolder, sub_nm, run_nm, task_fullName);
        % fit decrease of peak force per effort chosen across time
        for iEch = Ech_levels
            jEch = iEch + 1;
            Ech_idx = (E_chosen_tmp == iEch).*(~isnan(latency_tmp.allTrials)) == 1;
            betas_tmp = glmfit(trialN(Ech_idx), latency_tmp.allTrials(Ech_idx),'normal');
            b_forceLatency_Sub(:,jEch,iRun) = betas_tmp;
            forceLatency_fit_tmp = glmval(betas_tmp, trialN(Ech_idx), 'identity');
            % extract bins
            [forceLatency_tmp.(run_nm_bis)(1,:,jEch),...
                trialN_tmp.(run_nm_bis)(1,:,jEch)] = do_bin2(latency_tmp.allTrials(Ech_idx), trialN(Ech_idx), n_trialBins,0);
            latency_fit_tmp.(run_nm_bis)(1,:,jEch) = do_bin2(forceLatency_fit_tmp, trialN(Ech_idx), n_trialBins,0);
        end % effort level chosen loop
        % do it as well for pool across all efforts
        okTrials = ~isnan(latency_tmp.allTrials);
        betas_tmp = glmfit(trialN(okTrials), latency_tmp.allTrials(okTrials),'normal');
        b_forceLatency_allE_Sub(:,iRun) = betas_tmp;
        forceLatency_fit_tmp = glmval(betas_tmp, trialN, 'identity');
        % extract bins
        [forceLatency_allE_tmp.(run_nm_bis)(1,:),...
            trialN_allE_tmp.(run_nm_bis)(1,:)] = do_bin2(latency_tmp.allTrials, trialN, n_trialBins,0);
        forceLatency_fit_allE_tmp.(run_nm_bis)(1,:) = do_bin2(forceLatency_fit_tmp, trialN, n_trialBins,0);
    end % run loop

    %% average [decrease of peak force per effort chosen across time] slope
    % across runs
    b_forceLatency_perEch_f_time(:,iS) = mean(b_forceLatency_Sub(2,:,:),3,'omitnan');
    for iEch = Ech_levels
        Ech_nm = ['Ech',num2str(iEch)];
        jEch = iEch + 1;
        forcePeak_f_time.(Ech_nm)(:,iS) = mean([forceLatency_tmp.run1(1,:,jEch); forceLatency_tmp.run2(1,:,jEch)],1,'omitnan');
        trialN_f_time.(Ech_nm)(:,iS) = mean([trialN_tmp.run1(1,:,jEch); trialN_tmp.run2(1,:,jEch)],1,'omitnan');
        forcePeak_fit_f_time.(Ech_nm)(:,iS) = mean([latency_fit_tmp.run1(1,:,jEch); latency_fit_tmp.run2(1,:,jEch)],1,'omitnan');
    end
    b_forceLatency_allE_f_time(1,iS) = mean(b_forceLatency_allE_Sub(2,:),2,'omitnan');
    forcePeak_f_time.allE(:,iS) = mean([forceLatency_allE_tmp.run1(1,:); forceLatency_allE_tmp.run2(1,:)],1,'omitnan');
    trialN_f_time.allE(:,iS) = mean([trialN_allE_tmp.run1(1,:); trialN_allE_tmp.run2(1,:)],1,'omitnan');
    forcePeak_fit_f_time.allE(:,iS) = mean([forceLatency_fit_allE_tmp.run1(1,:); forceLatency_fit_allE_tmp.run2(1,:)],1,'omitnan');
end % subject loop

%% average across subjects
for iEch = Ech_levels
    Ech_nm = ['Ech',num2str(iEch)];
    [m_forcePeak_f_time.(Ech_nm), sem_forcePeak_f_time.(Ech_nm)] = mean_sem_sd(forcePeak_f_time.(Ech_nm),2);
    [m_trialN_f_time.(Ech_nm), sem_trialN_f_time.(Ech_nm)] = mean_sem_sd(trialN_f_time.(Ech_nm),2);
    [m_forcePeak_fit_f_time.(Ech_nm), sem_forcePeak_fit_f_time.(Ech_nm)] = mean_sem_sd(forcePeak_fit_f_time.(Ech_nm),2);

    % check how significant the decrease is
    jEch = iEch + 1;
    [~,pval.(Ech_nm)] = ttest(b_forceLatency_perEch_f_time(jEch, :));
end
% all efforts 
[m_forcePeak_f_time.allE, sem_forcePeak_f_time.allE] = mean_sem_sd(forcePeak_f_time.allE,2);
[m_trialN_f_time.allE, sem_trialN_f_time.allE] = mean_sem_sd(trialN_f_time.allE,2);
[m_forcePeak_fit_f_time.allE, sem_forcePeak_fit_f_time.allE] = mean_sem_sd(forcePeak_fit_f_time.allE,2);
[~,pval.allE] = ttest(b_forceLatency_allE_f_time);

%% figure: performance peak decrease over time
if figDisp == 1
    pSize = 30;
    lWidth = 3;
    grey = [143 143 143]./255;

    % split per effort level
    fig;
    for iEch = Ech_levels
        Ech_nm = ['Ech',num2str(iEch)];
        jEch = iEch + 1;
        subplot(2,2,jEch);
        hold on;
        er_hdl = errorbar(m_trialN_f_time.(Ech_nm),...
            m_forcePeak_f_time.(Ech_nm),...
            sem_forcePeak_f_time.(Ech_nm));
        er_hdl.LineWidth = lWidth;
        er_hdl.Color = 'k';
        er_hdl.LineStyle = 'none';
        fit_hdl = plot(m_trialN_f_time.(Ech_nm), m_forcePeak_fit_f_time.(Ech_nm));
        fit_hdl.LineWidth = lWidth;
        fit_hdl.Color = grey;
        fit_hdl.LineStyle = '--';
        xlabel('Trial number');
        ylabel({'Latency (s)';Ech_nm});
        legend_size(pSize);
    end % effort chosen

    % general effect of time
    fig;
    hold on;
    er_hdl = errorbar(m_trialN_f_time.allE,...
        m_forcePeak_f_time.allE,...
        sem_forcePeak_f_time.allE);
    er_hdl.LineWidth = lWidth;
    er_hdl.Color = 'k';
    er_hdl.LineStyle = 'none';
    fit_hdl = plot(m_trialN_f_time.allE, m_forcePeak_fit_f_time.allE);
    fit_hdl.LineWidth = lWidth;
    fit_hdl.Color = grey;
    fit_hdl.LineStyle = '--';
    xlabel('Trial number');
    ylabel('Latency (s)');
    legend_size(pSize);
end

end % function