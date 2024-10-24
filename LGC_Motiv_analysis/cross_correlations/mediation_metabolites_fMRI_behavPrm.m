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
% + extract corresponding metabolites across individuals
[metabolite_allSubs, MRS_ROI_nm, metabolite_nm] = metabolite_extraction(study_nm, subject_id);

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
% con_nm = con_names{con_idx};
con_nm = inputdlg('Contrast short name?');
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

%% perform the mediation
goodSubs = ~isnan(metabolite_allSubs);

X_nm = [MRS_ROI_nm,'-',metabolite_nm];
M_nm = ['fMRI-',ROI_coords.ROI_nm.ROI_1_shortName,'-',con_nm{1}];
Y_nm = prm_nm;
[a,b,c,c_prime, pval] = mediation(metabolite_allSubs(goodSubs),...
    con_data(goodSubs),...
    behavPrm(goodSubs),...
    X_nm, M_nm, Y_nm);