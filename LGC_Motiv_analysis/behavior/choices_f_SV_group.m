function[] = choices_f_SV_group()
% [] = choices_f_SV_group()
% choices_f_SV_group: display choices in function of the subjective value 
% extracted by choices_f_SV.m

%% subject selection
[study_nm, condition, ~, subject_id, NS] = sub_id;
%% working directory
rootPath = LGCM_root_paths;
study_path = [rootPath, study_nm, filesep];

%% main  parameters
mdl_nm = 'mdl_3'; % bayesian model to extract
nBins = 6;
[choices_f_deltaNV.Ep,...
    choices_f_deltaNV.Em,...
    choices_f_deltaNV.EpEm,...
    pChoice_f_deltaNV.Ep,...
    pChoice_f_deltaNV.Em,...
    pChoice_f_deltaNV.EpEm] = deal(NaN(nBins, NS));

%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    [~, ~,...
        choices_f_deltaNV_tmp, deltaNV_f_deltaNV_tmp,...
        pChoice_f_deltaNV_tmp] = choices_f_SV(study_nm, sub_nm,...
        condition, study_path, nBins, mdl_nm);
    
    % physical effort
    deltaNV_f_deltaNV.Ep(:,iS) = deltaNV_f_deltaNV_tmp.Ep.aRuns;
    choices_f_deltaNV.Ep(:,iS) = choices_f_deltaNV_tmp.Ep.aRuns;
    pChoice_f_deltaNV.Ep(:,iS) = pChoice_f_deltaNV_tmp.Ep.aRuns;
    % mental effort
    deltaNV_f_deltaNV.Em(:,iS) = deltaNV_f_deltaNV_tmp.Em.aRuns;
    choices_f_deltaNV.Em(:,iS) = choices_f_deltaNV_tmp.Em.aRuns;
    pChoice_f_deltaNV.Em(:,iS) = pChoice_f_deltaNV_tmp.Em.aRuns;
    % physical + mental effort
    deltaNV_f_deltaNV.EpEm(:,iS) = deltaNV_f_deltaNV_tmp.EpEm.aRuns;
    choices_f_deltaNV.EpEm(:,iS) = choices_f_deltaNV_tmp.EpEm.aRuns;
    pChoice_f_deltaNV.EpEm(:,iS) = pChoice_f_deltaNV_tmp.EpEm.aRuns;
end % subject loop

%% average across subjects
% delta(NV)
[m_deltaNV_f_deltaNV.Ep, ~] = mean_sem_sd(deltaNV_f_deltaNV.Ep,2);
[m_deltaNV_f_deltaNV.Em, ~] = mean_sem_sd(deltaNV_f_deltaNV.Em,2);
[m_deltaNV_f_deltaNV.EpEm, ~] = mean_sem_sd(deltaNV_f_deltaNV.EpEm,2);
% choices
[m_choices_f_deltaNV.Ep, sem_choices_f_deltaNV.Ep] = mean_sem_sd(choices_f_deltaNV.Ep,2);
[m_choices_f_deltaNV.Em, sem_choices_f_deltaNV.Em] = mean_sem_sd(choices_f_deltaNV.Em,2);
[m_choices_f_deltaNV.EpEm, sem_choices_f_deltaNV.EpEm] = mean_sem_sd(choices_f_deltaNV.EpEm,2);
% p(choices)
[m_pChoice_f_deltaNV.Ep, ~] = mean_sem_sd(pChoice_f_deltaNV.Ep,2);
[m_pChoice_f_deltaNV.Em, ~] = mean_sem_sd(pChoice_f_deltaNV.Em,2);
[m_pChoice_f_deltaNV.EpEm, ~] = mean_sem_sd(pChoice_f_deltaNV.EpEm,2);

%% display figures
[pSize, ~, col, ~] = general_fig_prm;
lWidth = 1;
%% Ep/Em (each in one handle)
fig;

% Ep
subplot(1,2,1); hold on;
choice_hdl = errorbar(m_deltaNV_f_deltaNV.Ep,...
    m_choices_f_deltaNV.Ep,...
    sem_choices_f_deltaNV.Ep);
pChoice_hdl = plot(m_deltaNV_f_deltaNV.Ep, m_pChoice_f_deltaNV.Ep);
% figure properties
errorbar_hdl_upgrade(choice_hdl, col.black);
fit_hdl_upgrade(pChoice_hdl, col.grey)
legend_size(pSize);
xlabel('NV');
ylabel('Choices (%)');

% Em
subplot(1,2,2); hold on;
choice_hdl = errorbar(m_deltaNV_f_deltaNV.Em,...
    m_choices_f_deltaNV.Em,...
    sem_choices_f_deltaNV.Em);
pChoice_hdl = plot(m_deltaNV_f_deltaNV.Em, m_pChoice_f_deltaNV.Em);
% figure properties
errorbar_hdl_upgrade(choice_hdl, col.black);
fit_hdl_upgrade(pChoice_hdl, col.grey)
legend_size(pSize);
xlabel('NV');
ylabel('Choices (%)');


%% Ep+Em (total choices)
fig;
choice_hdl = errorbar(m_deltaNV_f_deltaNV.EpEm,...
    m_choices_f_deltaNV.EpEm,...
    sem_choices_f_deltaNV.EpEm);
pChoice_hdl = plot(m_deltaNV_f_deltaNV.EpEm, m_pChoice_f_deltaNV.EpEm);
% figure properties
errorbar_hdl_upgrade(choice_hdl, col.black);
fit_hdl_upgrade(pChoice_hdl, col.grey)
legend_size(pSize);
xlabel('NV');
ylabel('Choices (%)');

%% Ep/Em in the same graph
Ep_col = [215 25 28]./255;
Em_col = [166 217 106]./255;
% Ep/Em
fig;

% show indifference point
xlim_vals = [-5 10];
ylim_vals = [0 100];
line([0 0],ylim_vals, 'LineWidth',lWidth, 'Color',col.black);
line(xlim_vals,[0.5 0.5], 'LineWidth',lWidth, 'Color',col.black);

% Ep
Ep_choice_hdl = errorbar(m_deltaNV_f_deltaNV.Ep,...
    m_choices_f_deltaNV.Ep.*100,...
    sem_choices_f_deltaNV.Ep.*100);
Ep_pChoice_hdl = plot(m_deltaNV_f_deltaNV.Ep,...
    m_pChoice_f_deltaNV.Ep.*100);
% figure properties
errorbar_hdl_upgrade(Ep_choice_hdl, Ep_col);
fit_hdl_upgrade(Ep_pChoice_hdl, Ep_col);

% Em
Em_choice_hdl = errorbar(m_deltaNV_f_deltaNV.Em,...
    m_choices_f_deltaNV.Em.*100,...
    sem_choices_f_deltaNV.Em.*100);
Em_pChoice_hdl = plot(m_deltaNV_f_deltaNV.Em,...
    m_pChoice_f_deltaNV.Em.*100);
% figure properties
errorbar_hdl_upgrade(Em_choice_hdl, Em_col);
fit_hdl_upgrade(Em_pChoice_hdl, Em_col);
xlabel('ΔSV + bias');
ylabel('Choices (%)');
legend([Ep_choice_hdl, Em_choice_hdl], 'physical','mental',...
    'Location','NorthWest');
legend('boxoff');
xlim(xlim_vals);
ylim(ylim_vals);
legend_size(60);

end % function