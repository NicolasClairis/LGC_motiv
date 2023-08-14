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

%% extract proportion of choices across individuals and tasks
fig_disp = 0;
[choice_hE] = choice_hE_proportion(study_nm, condition, subject_id, fig_disp);
parameter_names = {'Ep','Em','EpEm'};
nPrm = length(parameter_names);

%% launch this before to avoid case where nothing is significant for the
% current BOLD contrast but the information from the previous test was kept
clear('mediation_path','pval','N_goodSubs');
%% perform the mediation
pval.signif = struct;
dispMed = 0; % do not display mediation (too many plots)
for iROI = 1:nROIs
    MRS_ROI_nm = MRS_ROIs{iROI};
    for iMb = 1:n_metabolites.(MRS_ROI_nm)
        metabolite_nm = metabolite_names.(MRS_ROI_nm){iMb};
        metabolite_nm_bis = strrep(metabolite_nm,'_div_','/');
        metabolite_allSubs = metabolites.(MRS_ROI_nm).(metabolite_nm);
        goodSubs = ~isnan(metabolite_allSubs);
        
        for iPrm = 1:nPrm
            prm_nm = parameter_names{iPrm};
            behavPrm = choice_hE.(prm_nm);
            
            X_nm = [MRS_ROI_nm,'-',metabolite_nm_bis];
            M_nm = ['fMRI-',ROI_coords.ROI_nm.ROI_1_shortName,'-',con_nm{1}];
            Y_nm = prm_nm;
            [mediation_path.(MRS_ROI_nm).(metabolite_nm).(prm_nm).a,...
                mediation_path.(MRS_ROI_nm).(metabolite_nm).(prm_nm).b,...
                mediation_path.(MRS_ROI_nm).(metabolite_nm).(prm_nm).c,...
                mediation_path.(MRS_ROI_nm).(metabolite_nm).(prm_nm).c_prime,...
                pval.(MRS_ROI_nm).(metabolite_nm).(prm_nm),...
                stats.(MRS_ROI_nm).(metabolite_nm).(prm_nm)] = mediation(metabolite_allSubs(goodSubs),...
                con_data(goodSubs),...
                behavPrm(goodSubs),...
                X_nm, M_nm, Y_nm, dispMed);
            % store how many subjects were kept
            N_goodSubs.(MRS_ROI_nm).(metabolite_nm) = sum(goodSubs);
            
            % store when significant
            if pval.(MRS_ROI_nm).(metabolite_nm).(prm_nm).a < 0.05 &&...
                    pval.(MRS_ROI_nm).(metabolite_nm).(prm_nm).b < 0.05
                pval.signif.([MRS_ROI_nm,'_',metabolite_nm]).(prm_nm) = max(pval.(MRS_ROI_nm).(metabolite_nm).(prm_nm).a,...
                    pval.(MRS_ROI_nm).(metabolite_nm).(prm_nm).b);
            end
            
            
            %% perform the same but removing "outliers" (><mean*3SD)
            [~, ~, metabolite_clean] = rmv_outliers_3sd(metabolite_allSubs);
            [~, ~, con_data_clean] = rmv_outliers_3sd(con_data);
            [~, ~, behavPrm_clean] = rmv_outliers_3sd(behavPrm);
            goodSubs_bis = ~isnan(metabolite_clean).*~isnan(con_data_clean).*~isnan(behavPrm_clean) == 1;
            
            [mediation_path.no_outliers.(MRS_ROI_nm).(metabolite_nm).(prm_nm).a,...
                mediation_path.no_outliers.(MRS_ROI_nm).(metabolite_nm).(prm_nm).b,...
                mediation_path.no_outliers.(MRS_ROI_nm).(metabolite_nm).(prm_nm).c,...
                mediation_path.no_outliers.(MRS_ROI_nm).(metabolite_nm).(prm_nm).c_prime,...
                pval.no_outliers.(MRS_ROI_nm).(metabolite_nm).(prm_nm),...
                stats.no_outliers.(MRS_ROI_nm).(metabolite_nm).(prm_nm)] = mediation(metabolite_allSubs(goodSubs_bis),...
                con_data(goodSubs_bis),...
                behavPrm(goodSubs_bis),...
                X_nm, M_nm, Y_nm, dispMed);
            % store how many subjects were kept
            N_goodSubs.no_outliers.(MRS_ROI_nm).(metabolite_nm).(prm_nm) = sum(goodSubs_bis);
            
            % store when significant
            if pval.no_outliers.(MRS_ROI_nm).(metabolite_nm).(prm_nm).a < 0.05 &&...
                    pval.no_outliers.(MRS_ROI_nm).(metabolite_nm).(prm_nm).b < 0.05
                pval.no_outliers.signif.([MRS_ROI_nm,'_',metabolite_nm]).(prm_nm) = max(pval.no_outliers.(MRS_ROI_nm).(metabolite_nm).(prm_nm).a,...
                    pval.no_outliers.(MRS_ROI_nm).(metabolite_nm).(prm_nm).b);
            end
        end % parameter loop
    end % metabolites loop
end % ROI loop

%% line to launch to run again on a different contrast
% launch this before to avoid case where nothing is significant for the
% current BOLD contrast but the information from the previous test was kept
% clear('mediation_path','pval','N_goodSubs');

%% lines to launch to display metabolite of interest
MRS_ROI_nm='dmPFC';
metabolite_nm='Glu_div_GSH';
metabolite_nm_bis = strrep(metabolite_nm,'_div_','/');
prm_nm='Ep';
metabolite_allSubs = metabolites.(MRS_ROI_nm).(metabolite_nm);
        goodSubs = ~isnan(metabolite_allSubs);
behavPrm = choice_hE.(prm_nm);
X_nm = [MRS_ROI_nm,'-',metabolite_nm_bis];
M_nm='dmPFC';
Y_nm = prm_nm;
dispMed = 1;
[~, ~, metabolite_clean] = rmv_outliers_3sd(metabolite_allSubs);
[~, ~, con_data_clean] = rmv_outliers_3sd(con_data);
[~, ~, behavPrm_clean] = rmv_outliers_3sd(behavPrm);
goodSubs_bis = ~isnan(metabolite_clean).*~isnan(con_data_clean).*~isnan(behavPrm_clean) == 1;

[mediation_path.no_outliers.a.(MRS_ROI_nm).(metabolite_nm).(prm_nm),...
    mediation_path.no_outliers.b.(MRS_ROI_nm).(metabolite_nm).(prm_nm),...
    mediation_path.no_outliers.c.(MRS_ROI_nm).(metabolite_nm).(prm_nm),...
    mediation_path.no_outliers.c_prime.(MRS_ROI_nm).(metabolite_nm).(prm_nm),...
    pval.no_outliers.(MRS_ROI_nm).(metabolite_nm).(prm_nm)] = mediation(metabolite_allSubs(goodSubs_bis),...
    con_data(goodSubs_bis),...
    behavPrm(goodSubs_bis),...
    X_nm, M_nm, Y_nm, dispMed);