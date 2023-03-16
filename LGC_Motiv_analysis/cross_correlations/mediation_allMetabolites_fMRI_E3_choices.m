%% script performing the mediation between brain metabolites and 
% the proportion of choices for the E3 option through BOLD regression 
% estimates of specific regions.
% You need to define all participants of the mediation (which GLM for fMRI, 
% which regression estimate of the fMRI GLM). Then the script will perform 
% the mediation for you. This script will test all the metabolites of all 
% brain areas, while mediation_metabolites_fMRI_behavPrm.m is targeted to 
% one single metabolite.

%% define subjects to include
study_nm = 'study1';
condition = subject_condition();
[subject_id, NS] = LGCM_subject_selection(study_nm, condition, 'all');

%% load metabolites for all individuals and all brain areas
[metabolites] = metabolite_load(subject_id);
switch study_nm
    case 'study1'
        MRS_ROIs = {'dmPFC','aIns'};
    case 'study2'
        error('not ready yet');
end
nROIs = length(MRS_ROIs);
for iROI = 1:nROIs
    metabolite_names.(MRS_ROIs{iROI}) = fieldnames(metabolites.(MRS_ROIs{iROI}));
    n_metabolites.(MRS_ROIs{iROI}) = length(metabolite_names.(MRS_ROIs{iROI}));
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

%% extract proportion of high E choices for each E level
tasks = {'Ep','Em'};
task_idx = listdlg('promptstring','Which behavioral task?',...
    'ListString',tasks);
bhv_task_nm = tasks{task_idx};
[choice_hE] = extract_proportion_choice_hE_perSub(study_nm, subject_id,...
    condition, bhv_task_nm);
choice_hE_E3 = choice_hE.(bhv_task_nm).allSubs(3,:);
bhv_nm = [bhv_task_nm,'_E3'];
%% perform the mediation
pval.signif = struct;
dispMed = 0; % do not display mediation (too many plots)
for iROI = 1:nROIs
    MRS_ROI_nm = MRS_ROIs{iROI};
    for iMb = 1:n_metabolites.(MRS_ROI_nm)
        metabolite_nm = metabolite_names.(MRS_ROI_nm){iMb};
        metabolite_allSubs = metabolites.(MRS_ROI_nm).(metabolite_nm);
        goodSubs = ~isnan(metabolite_allSubs);
        
        X_nm = [MRS_ROI_nm,'-',metabolite_nm];
        M_nm = ['fMRI-',ROI_coords.ROI_nm.ROI_1_shortName,'-',con_nm{1}];
        Y_nm = ['choice hE - E3 - ',bhv_task_nm,' task'];
        [mediation_path.a.(MRS_ROI_nm).(metabolite_nm).(bhv_nm),...
            mediation_path.b.(MRS_ROI_nm).(metabolite_nm).(bhv_nm),...
            mediation_path.c.(MRS_ROI_nm).(metabolite_nm).(bhv_nm),...
            mediation_path.c_prime.(MRS_ROI_nm).(metabolite_nm).(bhv_nm),...
            pval.(MRS_ROI_nm).(metabolite_nm).(bhv_nm)] = mediation(metabolite_allSubs(goodSubs),...
            con_data(goodSubs),...
            choice_hE_E3(goodSubs),...
            X_nm, M_nm, Y_nm, dispMed);
        
        % store when significant
        if pval.(MRS_ROI_nm).(metabolite_nm).(bhv_nm).a < 0.05 &&...
                pval.(MRS_ROI_nm).(metabolite_nm).(bhv_nm).b < 0.05
            pval.signif.([MRS_ROI_nm,'_',metabolite_nm]).(bhv_nm) = max(pval.(MRS_ROI_nm).(metabolite_nm).(bhv_nm).a,...
                pval.(MRS_ROI_nm).(metabolite_nm).(bhv_nm).b);
        end
        
        
        %% perform the same but removing "outliers" (><mean*3SD)
        [~, ~, metabolite_clean] = rmv_outliers_3sd(metabolite_allSubs);
        [~, ~, con_data_clean] = rmv_outliers_3sd(con_data);
        [~, ~, behavPrm_clean] = rmv_outliers_3sd(choice_hE_E3);
        goodSubs_bis = ~isnan(metabolite_clean).*~isnan(con_data_clean).*~isnan(behavPrm_clean) == 1;
        
        [mediation_path.no_outliers.a.(MRS_ROI_nm).(metabolite_nm).(bhv_nm),...
            mediation_path.no_outliers.b.(MRS_ROI_nm).(metabolite_nm).(bhv_nm),...
            mediation_path.no_outliers.c.(MRS_ROI_nm).(metabolite_nm).(bhv_nm),...
            mediation_path.no_outliers.c_prime.(MRS_ROI_nm).(metabolite_nm).(bhv_nm),...
            pval.no_outliers.(MRS_ROI_nm).(metabolite_nm).(bhv_nm)] = mediation(metabolite_allSubs(goodSubs_bis),...
            con_data(goodSubs_bis),...
            choice_hE_E3(goodSubs_bis),...
            X_nm, M_nm, Y_nm, dispMed);
        
        % store when significant
        if pval.no_outliers.(MRS_ROI_nm).(metabolite_nm).(bhv_nm).a < 0.05 &&...
                pval.no_outliers.(MRS_ROI_nm).(metabolite_nm).(bhv_nm).b < 0.05
            pval.no_outliers.signif.([MRS_ROI_nm,'_',metabolite_nm]).(bhv_nm) = max(pval.no_outliers.(MRS_ROI_nm).(metabolite_nm).(bhv_nm).a,...
                pval.no_outliers.(MRS_ROI_nm).(metabolite_nm).(bhv_nm).b);
        end
    end % metabolites loop
end % ROI loop