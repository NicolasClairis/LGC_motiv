function[betas, pval] = mental_nMax_vs_bhvPrm()
% [betas, pval] = mental_nMax_vs_bhvPrm()
% mental_nMax_vs_bhvPrm compare n.max of calibration and behavioral 
% parameters.
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
        dataRoot = ['E:',filesep];
    case 'Home'
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

%% loop through subjects
[nMax_all] = deal(NaN(1,NS));
for iBhvPrm = 1:n_bhvPrm
    prm.(prm_names{iBhvPrm}) = NaN(1,NS);
end
for iS = 1:NS
    sub_nm = subject_id{iS};
    
    % extract calibrated Fmax
    subDataPath = [dataRoot,study_nm,filesep,'CID',sub_nm,filesep,'behavior',filesep];
    nMax_all(iS) = getfield(load([subDataPath,'CID',sub_nm,'_mentalCalib.mat'],'NMP'),'NMP'); % maximum number of answers
    
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
    goodSubs = ~isnan(prm.(prm_nm)).*~isnan(nMax_all) == 1;
    
    % extract people in ascending order
    prm_nMax_all_ascOrder = sort(prm.(prm_nm)(goodSubs));
    
    % MVC = f(parameter)
    [betas.nMax_all.(prm_nm),~,stats_nMax_all.(prm_nm)] = glmfit(prm.(prm_nm)(goodSubs), nMax_all(goodSubs), 'normal');
    pval.nMax.(prm_nm) = stats_nMax_all.(prm_nm).p;
    nMax_all_fit.(prm_nm) = glmval(betas.nMax_all.(prm_nm), prm_nMax_all_ascOrder, 'identity');
    
    %% display result
    
    % all subjects
    fig;
    
    % MVC = f(parameter)
    scat_hdl = scatter(prm.(prm_nm)(goodSubs), nMax_all(goodSubs));
    hold on;
    scat_hdl.LineWidth = lWidth;
    scat_hdl.MarkerEdgeColor = black;
    fit_hdl = plot(prm_nMax_all_ascOrder, nMax_all_fit.(prm_nm));
    fit_hdl.LineWidth = lWidth;
    fit_hdl.Color = grey;
    fit_hdl.LineStyle = '--';
    xlabel(prm_nm);
    ylabel('Calibration number of good answers');
    legend_size(pSize);
end % parameter loop

end % function