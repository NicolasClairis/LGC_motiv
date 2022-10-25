%% script performing the mediation between brain metabolites and 
% behavioural parameters through BOLD regression estimates of specific
% regions.
% You need to define all participants of the mediation (which GLM for fMRI, which
% regression estimate of the fMRI GLM, which behavioural model and which
% behavioural parameter). Then the script will perform the mediation for
% you. This script will test all the metabolites of all brain areas, while
% mediation_metabolites_fMRI_behavPrm.m is targeted to one single
% metabolite.

%% define subjects to include
study_nm = 'study1';
condition = subject_condition();
[subject_id, NS] = LGCM_subject_selection(study_nm, condition, 'all');

%% load metabolites for all individuals and all brain areas
[metabolites] = metabolite_load(subject_id);
switch study_nm
    case 'study1'
        ROIs = {'dmPFC','aIns'};
    case 'study2'
        error('not ready yet');
end
nROIs = length(ROIs);
for iROI = 1:nROIs
    metabolite_names.(ROIs{iROI}) = fieldnames(metabolites.(ROIs{iROI}));
    n_metabolites.(ROIs{iROI}) = length(metabolite_names.(ROIs{iROI}));
end % roi loop

%% define fMRI GLM to work on
GLM_str = inputdlg('Which fMRI GLM?');
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
for iROI = 1:nROIs
    ROI_nm = ROIs{iROI};
    for iMb = 1:n_metabolites
        metabolite_nm = metabolites{iMb};
        switch metabolite_nm
            case 'Glu_Gln'
                metabolite_nm_bis = 'Glx';
            otherwise
                metabolite_nm_bis = metabolite_nm;
        end
        metabolite_allSubs = metabolites.(ROI_nm).(metabolite_nm);
        goodSubs = ~isnan(metabolite_allSubs);
        
        X_nm = [MRS_ROI_nm,'-',metabolite_nm_bis];
        M_nm = ['fMRI-',ROI_coords.ROI_nm.ROI_1_shortName,'-',con_nm{1}];
        Y_nm = prm_nm;
        [a.(ROI_nm).(metabolite_nm),...
            b.(ROI_nm).(metabolite_nm),...
            c.(ROI_nm).(metabolite_nm),...
            c_prime.(ROI_nm).(metabolite_nm),...
            pval.(ROI_nm).(metabolite_nm)] = mediation(metabolite_allSubs(goodSubs),...
            con_data(goodSubs),...
            behavPrm(goodSubs),...
            X_nm, M_nm, Y_nm);
    end % metabolites loop
end % ROI loop