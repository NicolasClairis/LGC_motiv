%% script performing the mediation between brain metabolites and 
% behavioural parameters through BOLD regression estimates of specific
% regions.
% You need to define all participants of the mediation (brain area to
% select for metabolites, which metabolite, which GLM for fMRI, which
% regression estimate of the fMRI GLM, which behavioural model and which
% behavioural parameter). Then the script will perform the mediation for
% you.

%% define subjects to include
study_nm = 'study1';
condition = subject_condition();
[subject_id, NS] = LGCM_subject_selection(study_nm, condition, 'all');

%% define metabolite and ROI you want to focus on
% ROI
if strcmp(study_nm,'study1')
    ROIs = {'dmPFC','aIns'};
    nROIs = length(ROIs);
    ROI_idx = spm_input('Metabolites in which brain area?',1,'m',...
        ROIs,1:nROIs,0);
    ROI_nm = ROIs{ROI_idx};
    % select metabolite of interest
    metabolites = {'Mac','Ala','Asp','PCho','Cr','PCr','GABA',...
        'Gln','Glu','GSH','Gly','Ins','Lac','NAA','Scyllo','Tau',...
        'Asc','Glc','NAAG','GPC','PE','Ser',...
        'NAA_NAAG','Glu_Gln','GPC_PCho','Cr_PCr','Gly_Ins','Gln_div_Glu'};
    n_met = length(metabolites);
    metabolite_idx = spm_input('Which metabolite to focus on?',1,'m',...
        metabolites,1:n_met,0);
    metabolite_nm = metabolites{metabolite_idx};
else
    error('not ready yet for study2');
end

%% extract corresponding metabolites across individuals
[metabolites] = metabolite_load(subject_id);
% focus on metabolite and brain area selected
metabolite_allSubs = metabolites.(ROI_nm).(metabolite_nm);

%% define GLM to work on
GLM_str = inputdlg('Which GLM?');
GLM = str2double(GLM_str);

%% define fMRI ROI to use
[con_vec_all,...
    ~, ~, ~,...
    con_names,...
    ROI_coords, ttest_ROI] = ROI_extraction_group('study1', GLM,...
    subject_id, condition, 0);
n_cons = size(con_vec_all, 1);
n_ROIs = size(con_vec_all,3);
if n_ROIs > 1
    error(['more than 1 ROI selected, mediation script cannot work that way',...
        'please focus on one and do it again.']);
end

%% define regression estimate to look for in the fMRI GLM
con_idx = listdlg('promptstring','Which contrast?',...
    'ListString',con_names);
con_nm = con_names{con_idx};
con_data = NaN(1,NS);
con_data(:) = con_vec_all(con_idx, :, 1);

%% extract behavioural parameters
[prm] = prm_extraction(subject_id);
parameters = fieldnames(prm);
behavPrm_CID = prm.CID;
parameter_names = parameters;
parameter_names(strcmp(parameter_names,'CID'))=[]; % remove indication of subject ID

%% define behavioural parameters to check
behavPrm_idx = listdlg('promptstring','Which behavioural parameter?',...
    'ListString',parameter_names);
prm_nm = parameter_names{behavPrm_idx};
behavPrm = prm.(prm_nm);

%% perform each path of the mediation
goodSubs = ~isnan(metabolite_allSubs);

% test correlation between metabolites and fMRI BOLD
[betas_1,~,stats_1] = glmfit(metabolite_allSubs(goodSubs), con_data(goodSubs),'normal');

% test correlation between fMRI BOLD, behavioural parameter and direct path
[betas_2,~,stats_2] = glmfit([con_data(goodSubs)',...
    metabolite_allSubs(goodSubs)'],...
    behavPrm(goodSubs),'normal');

% test also direct path
[betas_3,~,stats_3] = glmfit(metabolite_allSubs(goodSubs),...
    behavPrm(goodSubs),'normal');

% remove variance related to first path
con_data_residual = con_data(goodSubs) - betas_1(2).*metabolite_allSubs(goodSubs);
[betas_4,~,stats_4] = glmfit([con_data_residual',...
    metabolite_allSubs(goodSubs)'],...
    behavPrm(goodSubs),'normal');

% display relevant p.values
disp(['metabolites -> fMRI ',...
    'fMRI= ',num2str(round(betas_1(1),3)),...
    ' + (',num2str(round(betas_1(2),3)),')*metabolites; ',...
    'p=',num2str(stats_1.p(2))]);
disp(['fMRI -> behaviour ',...
    prm_nm,'= ',num2str(round(betas_2(1),3)),...
    ' + (',num2str(round(betas_2(2),3)),')*fMRI; ',...
    ' + (',num2str(round(betas_2(3),3)),')*metabolites; ',...
    'pfMRI=',num2str(stats_2.p(2)),';',...
    'pMetabolites=',num2str(stats_2.p(3))]);
disp(['metabolites -> behaviour ',...
    'behaviour= ',num2str(round(betas_3(1),3)),...
    ' + (',num2str(round(betas_3(2),3)),')*metabolites; ',...
    'p=',num2str(stats_3.p(2))]);