function[] = choices_f_time_group()
% [] = choices_f_time_group()
% choices_f_time_group: display choices in function of time 
% extracted by choices_f_time.m

%% subject selection
[study_nm, condition, ~, subject_id, NS] = sub_id;
%% working directory
rootPath = LGCM_root_paths;
study_path = [rootPath, study_nm, filesep];

%% main  parameters
mdl_nm = 'mdl_3'; % bayesian model to extract
nBins = 6;
[choices_f_trialN.Ep,...
    choices_f_trialN.Em,...
    choices_f_trialN.EpEm,...
    pChoice_f_trialN.Ep,...
    pChoice_f_trialN.Em,...
    pChoice_f_trialN.EpEm,...
    trialN_f_trialN.Ep,...
    trialN_f_trialN.Em,...
    trialN_f_trialN.EpEm] = deal(NaN(nBins, NS));
nTrials_perRun = 54;

%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    [choice_f_trialN_tmp, pChoice_f_trialN_tmp,...
        trialN_f_trialN_tmp] = choices_f_time(study_nm, sub_nm,...
        condition, study_path, nBins, mdl_nm);
    
    % physical effort
    trialN_f_trialN.Ep(:,iS) = trialN_f_trialN_tmp.Ep.aRuns;
    choices_f_trialN.Ep(:,iS) = choice_f_trialN_tmp.Ep.aRuns;
    pChoice_f_trialN.Ep(:,iS) = pChoice_f_trialN_tmp.Ep.aRuns;
    % mental effort
    trialN_f_trialN.Em(:,iS) = trialN_f_trialN_tmp.Em.aRuns;
    choices_f_trialN.Em(:,iS) = choice_f_trialN_tmp.Em.aRuns;
    pChoice_f_trialN.Em(:,iS) = pChoice_f_trialN_tmp.Em.aRuns;
    % physical + mental effort
    trialN_f_trialN.EpEm(:,iS) = trialN_f_trialN_tmp.EpEm.aRuns;
    choices_f_trialN.EpEm(:,iS) = choice_f_trialN_tmp.EpEm.aRuns;
    pChoice_f_trialN.EpEm(:,iS) = pChoice_f_trialN_tmp.EpEm.aRuns;
end % subject loop

%% average across subjects
% delta(NV)
[m_trialN_f_trialN.Ep, ~] = mean_sem_sd(trialN_f_trialN.Ep,2);
[m_trialN_f_trialN.Em, ~] = mean_sem_sd(trialN_f_trialN.Em,2);
[m_trialN_f_trialN.EpEm, ~] = mean_sem_sd(trialN_f_trialN.EpEm,2);
% choices
[m_choices_f_trialN.Ep, sem_choices_f_trialN.Ep] = mean_sem_sd(choices_f_trialN.Ep,2);
[m_choices_f_trialN.Em, sem_choices_f_trialN.Em] = mean_sem_sd(choices_f_trialN.Em,2);
[m_choices_f_trialN.EpEm, sem_choices_f_trialN.EpEm] = mean_sem_sd(choices_f_trialN.EpEm,2);
% p(choices)
[m_pChoice_f_trialN.Ep, ~] = mean_sem_sd(pChoice_f_trialN.Ep,2);
[m_pChoice_f_trialN.Em, ~] = mean_sem_sd(pChoice_f_trialN.Em,2);
[m_pChoice_f_trialN.EpEm, ~] = mean_sem_sd(pChoice_f_trialN.EpEm,2);

%% display figures
[pSize, ~, col, ~] = general_fig_prm;
lWidth = 1;
%% Ep/Em (each in one handle)
fig;

% Ep
subplot(1,2,1); hold on;
choice_hdl = errorbar(m_trialN_f_trialN.Ep,...
    m_choices_f_trialN.Ep,...
    sem_choices_f_trialN.Ep);
pChoice_hdl = plot(m_trialN_f_trialN.Ep, m_pChoice_f_trialN.Ep);
% figure properties
errorbar_hdl_upgrade(choice_hdl, col.black);
fit_hdl_upgrade(pChoice_hdl, col.grey)
legend_size(pSize);
xlabel('Trial');
ylabel('Choices (%)');

% Em
subplot(1,2,2); hold on;
choice_hdl = errorbar(m_trialN_f_trialN.Em,...
    m_choices_f_trialN.Em,...
    sem_choices_f_trialN.Em);
pChoice_hdl = plot(m_trialN_f_trialN.Em, m_pChoice_f_trialN.Em);
% figure properties
errorbar_hdl_upgrade(choice_hdl, col.black);
fit_hdl_upgrade(pChoice_hdl, col.grey)
legend_size(pSize);
xlabel('Trial');
ylabel('Choices (%)');


%% Ep+Em (total choices)
fig;
choice_hdl = errorbar(m_trialN_f_trialN.EpEm,...
    m_choices_f_trialN.EpEm,...
    sem_choices_f_trialN.EpEm);
pChoice_hdl = plot(m_trialN_f_trialN.EpEm, m_pChoice_f_trialN.EpEm);
% figure properties
errorbar_hdl_upgrade(choice_hdl, col.black);
fit_hdl_upgrade(pChoice_hdl, col.grey)
legend_size(pSize);
xlabel('Trial');
ylabel('Choices (%)');

%% Ep/Em in the same graph
Ep_col = [215 25 28]./255;
Em_col = [166 217 106]./255;
% Ep/Em
fig;

% show indifference point
xlim_vals = [0 nTrials_perRun];
ylim_vals = [40 80];

% Ep
Ep_choice_hdl = errorbar(m_trialN_f_trialN.Ep,...
    m_choices_f_trialN.Ep.*100,...
    sem_choices_f_trialN.Ep.*100);
Ep_pChoice_hdl = plot(m_trialN_f_trialN.Ep,...
    m_pChoice_f_trialN.Ep.*100);
% figure properties
errorbar_hdl_upgrade(Ep_choice_hdl, Ep_col);
fit_hdl_upgrade(Ep_pChoice_hdl, Ep_col);

% Em
Em_choice_hdl = errorbar(m_trialN_f_trialN.Em,...
    m_choices_f_trialN.Em.*100,...
    sem_choices_f_trialN.Em.*100);
Em_pChoice_hdl = plot(m_trialN_f_trialN.Em,...
    m_pChoice_f_trialN.Em.*100);
% figure properties
errorbar_hdl_upgrade(Em_choice_hdl, Em_col);
fit_hdl_upgrade(Em_pChoice_hdl, Em_col);
xlabel('Trial');
ylabel('Choices (%)');
legend([Ep_choice_hdl, Em_choice_hdl], 'physical','mental',...
    'Location','SouthWest');
legend('boxoff');
xlim(xlim_vals);
ylim(ylim_vals);
legend_size(60);

end % function