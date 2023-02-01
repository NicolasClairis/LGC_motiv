function[betas, pval] = grip_MVC_PCSA_vs_bhvPrm()
% [betas, pval] = grip_MVC_PCSA_vs_bhvPrm()
% grip_MVC_vs_bhvPrm compare maximum voluntary contraction force (MVC),
% theoretical force (thFmax) and difference between theoretical and actual 
% force(MVC-thFmax) and behavioral parameters.
%
% INPUTS
% 
% OUTPUTS
% betas: betas for the linear fit of each test
%
% pval: p.value for the linear fit of each test
%

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

%% load behavioral parameters
[parameters] = prm_extraction(study_nm,subject_id);
parameter_names = fieldnames(parameters);
behavPrm_CID = parameters.CID;
prm_names = parameter_names;
prm_names(strcmp(prm_names,'CID')) = []; % remove indication of subject ID
n_bhvPrm = length(prm_names);

%% load infos (including info about maximal theoretical force)
excelReadInfosFile = readtable([gitPath,study_nm,filesep,'summary_participants_infos.xlsx'],...
    'Sheet','Sheet1');

%% loop through subjects
[MVC_all, thFmax_all, MVC_min_thFmax_all] = deal(NaN(1,NS));
for iBhvPrm = 1:n_bhvPrm
    prm.(prm_names{iBhvPrm}) = NaN(1,NS);
end
for iS = 1:NS
    sub_nm = subject_id{iS};
    % extract theoretical Fmax
    sub_info_idx = strcmp(['CID',sub_nm],excelReadInfosFile.CID);
    pli_a = excelReadInfosFile.PliCutan_Ant_rieur_mm_(sub_info_idx);
    pli_p = excelReadInfosFile.PliCutan_Post_rieur_mm_(sub_info_idx);
    circ = excelReadInfosFile.Circonf_renceDeL_avant_bras_mm_(sub_info_idx);
    armLength = excelReadInfosFile.LongeurDeL_avant_bras_mm_(sub_info_idx);
    thFmax_all(iS) = Emax_morpho(pli_a, pli_p, circ, armLength);
    
    % extract calibrated Fmax
    subDataPath = [dataRoot,study_nm,filesep,'CID',sub_nm,filesep,'behavior',filesep];
    MVC_volts = getfield(load([subDataPath,'CID',sub_nm,'_physicalCalib.mat'],'MVC'),'MVC'); % MVC in volts
    MVC_all(iS) = grip_biopac_volts_to_newtons_conversion(MVC_volts); % convert in Newtons
    
    % extract difference
    MVC_min_thFmax_all(iS) = MVC_all(iS) - thFmax_all(iS);
    
    % extract behavioral parameter
    bhv_sub_idx = strcmp(behavPrm_CID,sub_nm);
    for iBhvPrm = 1:n_bhvPrm
        prm_nm = prm_names{iBhvPrm};
        prm.(prm_nm)(iS) = parameters.(prm_nm)(bhv_sub_idx);
    end
end % subject loop

%% figure parameters
pSize = 40;
lWidth = 3;
black = [0 0 0];
grey = [143 143 143]./255;

%% perform correlation
for iPrm = 1:n_bhvPrm
    prm_nm = prm_names{iPrm};
    
    % filter correct subjects
    goodSubs_MVC = ~isnan(prm.(prm_nm)).*~isnan(MVC_all) == 1;
    goodSubs_thFmax = ~isnan(prm.(prm_nm)).*~isnan(thFmax_all) == 1;
    goodSubs_MVC_thFmax = goodSubs_MVC.*goodSubs_thFmax == 1;
    
    % extract people in ascending order
    prm_MVC_all_ascOrder = sort(prm.(prm_nm)(goodSubs_MVC));
    prm_thFmax_all_ascOrder = sort(prm.(prm_nm)(goodSubs_thFmax));
    prm_MVC_min_thFmax_all_ascOrder = sort(prm.(prm_nm)(goodSubs_MVC_thFmax));
    
    % MVC = f(parameter)
    [betas.MVC_all.(prm_nm),~,stats_MVC_all.(prm_nm)] = glmfit(prm.(prm_nm)(goodSubs_MVC), MVC_all(goodSubs_MVC), 'normal');
    pval.MVC.(prm_nm) = stats_MVC_all.(prm_nm).p;
    MVC_all_fit.(prm_nm) = glmval(betas.MVC_all.(prm_nm), prm_MVC_all_ascOrder, 'identity');
    
    % theoretical Fmax = f(parameter)
    [betas.thFmax_all.(prm_nm),~,stats_thFmax_all.(prm_nm)] = glmfit(prm.(prm_nm)(goodSubs_thFmax), thFmax_all(goodSubs_thFmax), 'normal');
    pval.thFmax.(prm_nm) = stats_thFmax_all.(prm_nm).p;
    thFmax_all_fit.(prm_nm) = glmval(betas.thFmax_all.(prm_nm), prm_thFmax_all_ascOrder, 'identity');
    
    % MVC - thFmax = f(parameter)
    [betas.MVC_min_thFmax_all.(prm_nm),~,stats_MVC_min_thFmax_all.(prm_nm)] = glmfit(prm.(prm_nm)(goodSubs_MVC_thFmax), MVC_min_thFmax_all(goodSubs_MVC_thFmax), 'normal');
    pval.MVC_min_thFmax.(prm_nm) = stats_MVC_min_thFmax_all.(prm_nm).p;
    MVC_min_thFmax_all_fit.(prm_nm) = glmval(betas.MVC_min_thFmax_all.(prm_nm), prm_MVC_min_thFmax_all_ascOrder, 'identity');
    
    %% display result
    
    % all subjects
    fig;
    
    % MVC = f(parameter)
    subplot(1,3,1);
    scat_hdl = scatter(prm.(prm_nm)(goodSubs_MVC), MVC_all(goodSubs_MVC));
    hold on;
    scat_hdl.LineWidth = lWidth;
    scat_hdl.MarkerEdgeColor = black;
    fit_hdl = plot(prm_MVC_all_ascOrder, MVC_all_fit.(prm_nm));
    fit_hdl.LineWidth = lWidth;
    fit_hdl.Color = grey;
    fit_hdl.LineStyle = '--';
    xlabel(prm_nm);
    ylabel('Calibrated force (N)');
    legend_size(pSize);
    
    % theoretical Fmax = f(parameter)
    subplot(1,3,2);
    scat_hdl = scatter(prm.(prm_nm)(goodSubs_thFmax), thFmax_all(goodSubs_thFmax));
    hold on;
    scat_hdl.LineWidth = lWidth;
    scat_hdl.MarkerEdgeColor = black;
    fit_hdl = plot(prm_thFmax_all_ascOrder, thFmax_all_fit.(prm_nm));
    fit_hdl.LineWidth = lWidth;
    fit_hdl.Color = grey;
    fit_hdl.LineStyle = '--';
    xlabel(prm_nm);
    ylabel('Theoretical maximal force (N)');
    legend_size(pSize);
    
    % MVC - thFmax = f(parameter)
    subplot(1,3,3);
    scat_hdl = scatter(prm.(prm_nm)(goodSubs_MVC_thFmax), MVC_min_thFmax_all(goodSubs_MVC_thFmax));
    hold on;
    scat_hdl.LineWidth = lWidth;
    scat_hdl.MarkerEdgeColor = black;
    fit_hdl = plot(prm_MVC_min_thFmax_all_ascOrder, MVC_min_thFmax_all_fit.(prm_nm));
    fit_hdl.LineWidth = lWidth;
    fit_hdl.Color = grey;
    fit_hdl.LineStyle = '--';
    xlabel(prm_nm);
    ylabel('Calibrated force - Theoretical maximal force (N)');
    legend_size(pSize);
end % parameter loop

end % function