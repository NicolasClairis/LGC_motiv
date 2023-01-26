function[study_nm, cond, subject_id, sub_males, sub_females,...
    MVC_all, PCSA_all,...
    MVC_m, PCSA_m,...
    MVC_f, PCSA_f,...
    stats_all, stats_m, stats_f,...
    b_all, b_m, b_f] = grip_MVC_vs_PCSA(dispFig)
%[study_nm, cond, subject_id, sub_males, sub_females,...
%     MVC_all, PCSA_all,...
%     MVC_m, PCSA_m,...
%     MVC_f, PCSA_f,...
%     stats_all, stats_m, stats_f] = grip_MVC_vs_PCSA(dispFig)
% grip_MVC_vs_PCSA.m will test how much the theoretical maximal voluntary
% force inferred by measuring the forearm (PCSA) and the maximul voluntary
% contraction (MVC) force performed by the subjects during calibration of
% the physical grip task correlate with each other.
%
% INPUT
% dispFig: display figure (1) or not (0)
%
% OUTPUT
% study_nm: study name
%
% cond: condition used for the extraction
%
% subject_id: list of all subjects extracted
%
% sub_males: list of male subjects extracted
%
% sub_females: list of female subjects extracted
%
% MVC_all: maximum voluntary contraction force for all subjects
%
% PCSA_all: maximum theoretical force for all subjects
%
% MVC_m: maximum voluntary contraction force for male subjects only
%
% PCSA_m: maximum theoretical force for male subjects only
%
% MVC_f: maximum voluntary contraction force for female subjects only
%
% PCSA_f: maximum theoretical force for female subjects only
%
% stats_all: statistics for all subjects
%
% stats_m, stats_f: statistics for males (m) and females (f) only
%
% b_all, b_m, b_f: betas for the linear correlation for all, males (m) and
% females (f)

%% working directory
list_pcs = {'Lab','Home'};
which_pc_idx = listdlg('PromptString',{'Lab or home pc?'},...
    'SelectionMode','single','ListString',list_pcs);
curr_pc = list_pcs{which_pc_idx};
switch curr_pc
    case 'Lab'
        gitPath = [fullfile('C:','Users','clairis','Desktop','GitHub',...
            'LGC_motiv','LGC_Motiv_results'),filesep];
        dataRoot = ['E:',filesep];
    case 'Home'
        gitPath = [fullfile('C:','Users','Loco','Documents','GitHub',...
            'LGC_motiv','LGC_Motiv_results'),filesep];
        dataRoot = [fullfile('L:','human_data_private','raw_data_subject'),filesep];
end
%% load subjects
study_nm = 'study1';
cond = 'behavior';
[subject_id, NS] = LGCM_subject_selection(study_nm, cond);
% check also gender to know who is who
[sub_males, NS_m] = LGCM_subject_selection(study_nm, cond,'males');
[sub_females, NS_f] = LGCM_subject_selection(study_nm, cond,'females');

%% load infos (including info about maximal theoretical force)
excelReadInfosFile = readtable([gitPath,study_nm,filesep,'summary_participants_infos.xlsx'],...
    'Sheet','Sheet1');

%% loop through subjects
[MVC_all, PCSA_all] = deal(NaN(1,NS));
[MVC_m, PCSA_m] = deal(NaN(1,NS_m));
[MVC_f, PCSA_f] = deal(NaN(1,NS_f));
jM = 0;
jF = 0;
for iS = 1:NS
    sub_nm = subject_id{iS};
    % extract theoretical Fmax
    sub_info_idx = strcmp(['CID',sub_nm],excelReadInfosFile.CID);
    pli_a = excelReadInfosFile.PliCutan_Ant_rieur_mm_(sub_info_idx);
    pli_p = excelReadInfosFile.PliCutan_Post_rieur_mm_(sub_info_idx);
    circ = excelReadInfosFile.Circonf_renceDeL_avant_bras_mm_(sub_info_idx);
    length = excelReadInfosFile.LongeurDeL_avant_bras_mm_(sub_info_idx);
    PCSA_all(iS) = Emax_morpho(pli_a, pli_p, circ, length);
    % extract calibrated Fmax
    subDataPath = [dataRoot,study_nm,filesep,'CID',sub_nm,filesep,'behavior',filesep];
    MVC_volts = getfield(load([subDataPath,'CID',sub_nm,'_physicalCalib.mat'],'MVC'),'MVC'); % MVC in volts
    MVC_all(iS) = grip_biopac_volts_to_newtons_conversion(MVC_volts); % convert in Newtons
    
    % same but split by gender
    if ismember(sub_nm, sub_males)
        jM = jM + 1;
        PCSA_m(jM) = PCSA_all(iS);
        MVC_m(jM) = MVC_all(iS);
    elseif ismember(sub_nm, sub_females)
        jF = jF + 1;
        PCSA_f(jF) = PCSA_all(iS);
        MVC_f(jF) = MVC_all(iS);
    end
end % subject loop

%% perform correlation
% all subjects
[b_all,~,stats_all] = glmfit(PCSA_all, MVC_all, 'normal');
PCSA_all_ascOrder = sort(PCSA_all);
MVC_all_fit = glmval(b_all, PCSA_all_ascOrder, 'identity');
% males
[b_m,~,stats_m] = glmfit(PCSA_m, MVC_m, 'normal');
PCSA_m_ascOrder = sort(PCSA_m);
MVC_m_fit = glmval(b_m, PCSA_m_ascOrder, 'identity');

% females
[b_f,~,stats_f] = glmfit(PCSA_f, MVC_f, 'normal');
PCSA_f_ascOrder = sort(PCSA_f);
MVC_f_fit = glmval(b_all, PCSA_f_ascOrder, 'identity');

%% display result
if dispFig == 1
    pSize = 40;
    lWidth = 3;
    black = [0 0 0];
    grey = [143 143 143]./255;
    blue = [0 0 1];
    green = [0 1 0];
    lightBlue = [0 255 255]./255;
    lightGreen = [106 204 132]./255;
    
    % all subjects
    fig;
    scat_hdl = scatter(PCSA_all, MVC_all);
    scat_hdl.LineWidth = lWidth;
    scat_hdl.MarkerEdgeColor = black;
    fit_hdl = plot(PCSA_all_ascOrder, MVC_all_fit);
    fit_hdl.LineWidth = lWidth;
    fit_hdl.Color = grey;
    fit_hdl.LineStyle = '--';
    xlabel('Theoretical maximal force (N)');
    ylabel('Calibrated force (N)');
    legend_size(pSize);
    
    % split by gender
    fig;
    scat_m_hdl = scatter(PCSA_m, MVC_m);
    scat_m_hdl.LineWidth = lWidth;
    scat_m_hdl.MarkerEdgeColor = blue;
    scat_f_hdl = scatter(PCSA_f, MVC_f);
    scat_f_hdl.LineWidth = lWidth;
    scat_f_hdl.MarkerEdgeColor = green;
    fit_m_hdl = plot(PCSA_m_ascOrder, MVC_m_fit);
    fit_m_hdl.LineWidth = lWidth;
    fit_m_hdl.Color = lightBlue;
    fit_m_hdl.LineStyle = '--';
    fit_f_hdl = plot(PCSA_f_ascOrder, MVC_f_fit);
    fit_f_hdl.LineWidth = lWidth;
    fit_f_hdl.Color = lightGreen;
    fit_f_hdl.LineStyle = '--';
    legend([scat_m_hdl, scat_f_hdl],{'males','females'});
    legend('boxoff');
    legend('Location','NorthWest');
    xlabel('Theoretical maximal force (N)');
    ylabel('Calibrated force (N)');
    legend_size(pSize);
end % display figure
end % function