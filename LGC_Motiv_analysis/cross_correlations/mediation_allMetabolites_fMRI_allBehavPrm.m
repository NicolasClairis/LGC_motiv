%% script performing the mediation between brain metabolites and
% behavioural parameters through BOLD regression estimates of specific
% regions.
% You need to define all participants of the mediation (which GLM for fMRI, which
% regression estimate of the fMRI GLM, which behavioural model and which
% behavioural parameter). Then the script will perform the mediation for
% you. This script will test all the metabolites of all brain areas, while
% mediation_metabolites_fMRI_behavPrm.m is targeted to one single
% metabolite.

%% did you already launch the ROI extraction (1) or not (0)?
ROI_already_launched = 1;
if ~exist('con_vec_all','var') && ROI_already_launched ~= 0
    error(['ROI_already_launched = ',num2str(ROI_already_launched),...
        ' while data not in workspace. Please fix it.']);
end

if ROI_already_launched == 0
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
    
end % ROI already launched

% %% define regression estimate to look for in the fMRI GLM
% con_idx = listdlg('promptstring','Which contrast?',...
%     'ListString',con_names);
% % con_nm = con_names{con_idx};
% con_nm_str = inputdlg('Contrast short name?');
% con_nm = con_nm_str{1};
% con_data = NaN(1,NS);
% con_data(:) = con_vec_all(con_idx, :, 1);

%% extract behavioural parameters
[prm] = prm_extraction(study_nm, subject_id);
parameters = fieldnames(prm);
behavPrm_CID = prm.CID;
parameter_names = parameters;
parameter_names(strcmp(parameter_names,'CID'))=[]; % remove indication of subject ID
nPrm = length(parameter_names);

%% launch this before to avoid case where nothing is significant for the
% current BOLD contrast but the information from the previous test was kept
clear('mediation_path','pval','N_goodSubs','stats');
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
            behavPrm = prm.(prm_nm);
            
            X_nm = [MRS_ROI_nm,'-',metabolite_nm_bis];
            M_nm = ['fMRI-',ROI_coords.ROI_nm.ROI_1_shortName,'-',con_nm];
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
            
            % store when mediation is significant
            if pval.(MRS_ROI_nm).(metabolite_nm).(prm_nm).a < 0.05 &&...
                    pval.(MRS_ROI_nm).(metabolite_nm).(prm_nm).b < 0.05
                pval.signif.([MRS_ROI_nm,'_',metabolite_nm]).(prm_nm) = max(pval.(MRS_ROI_nm).(metabolite_nm).(prm_nm).a,...
                    pval.(MRS_ROI_nm).(metabolite_nm).(prm_nm).b);
            end
            
            % store when direct path is significant
            if pval.(MRS_ROI_nm).(metabolite_nm).(prm_nm).c < 0.05
                pval.direct_path.signif.([MRS_ROI_nm,'_',metabolite_nm]).(prm_nm) = pval.(MRS_ROI_nm).(metabolite_nm).(prm_nm).c;
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
            
            % store when mediation is significant
            if pval.no_outliers.(MRS_ROI_nm).(metabolite_nm).(prm_nm).a < 0.05 &&...
                    pval.no_outliers.(MRS_ROI_nm).(metabolite_nm).(prm_nm).b < 0.05
                pval.no_outliers.signif.([MRS_ROI_nm,'_',metabolite_nm]).(prm_nm) = max(pval.no_outliers.(MRS_ROI_nm).(metabolite_nm).(prm_nm).a,...
                    pval.no_outliers.(MRS_ROI_nm).(metabolite_nm).(prm_nm).b);
            end
            
            % store when direct path is significant
            if pval.no_outliers.(MRS_ROI_nm).(metabolite_nm).(prm_nm).c < 0.05
                pval.direct_path.no_outliers.signif.([MRS_ROI_nm,'_',metabolite_nm]).(prm_nm) = pval.no_outliers.(MRS_ROI_nm).(metabolite_nm).(prm_nm).c;
            end
            
            %% same but with boxcox transformation of behavioral parameters
            if ~strcmp(prm_nm,'kBiasM')
                behavPrm_boxcox = (boxcox(prm.(prm_nm)'))';
                [mediation_path.boxcox.(MRS_ROI_nm).(metabolite_nm).(prm_nm).a,...
                    mediation_path.boxcox.(MRS_ROI_nm).(metabolite_nm).(prm_nm).b,...
                    mediation_path.boxcox.(MRS_ROI_nm).(metabolite_nm).(prm_nm).c,...
                    mediation_path.boxcox.(MRS_ROI_nm).(metabolite_nm).(prm_nm).c_prime,...
                    pval.boxcox.(MRS_ROI_nm).(metabolite_nm).(prm_nm),...
                    stats.boxcox.(MRS_ROI_nm).(metabolite_nm).(prm_nm)] = mediation(metabolite_allSubs(goodSubs),...
                    con_data(goodSubs),...
                    behavPrm_boxcox(goodSubs),...
                    X_nm, M_nm, Y_nm, dispMed);
                % store how many subjects were kept
                N_goodSubs.boxcox.(MRS_ROI_nm).(metabolite_nm) = sum(goodSubs);
                
                % store when mediation is significant
                if pval.boxcox.(MRS_ROI_nm).(metabolite_nm).(prm_nm).a < 0.05 &&...
                        pval.boxcox.(MRS_ROI_nm).(metabolite_nm).(prm_nm).b < 0.05
                    pval.boxcox.signif.([MRS_ROI_nm,'_',metabolite_nm]).(prm_nm) = max(pval.boxcox.(MRS_ROI_nm).(metabolite_nm).(prm_nm).a,...
                        pval.boxcox.(MRS_ROI_nm).(metabolite_nm).(prm_nm).b);
                end
                
                % store when direct path is significant
                if pval.boxcox.(MRS_ROI_nm).(metabolite_nm).(prm_nm).c < 0.05
                    pval.direct_path.boxcox.signif.([MRS_ROI_nm,'_',metabolite_nm]).(prm_nm) = pval.boxcox.(MRS_ROI_nm).(metabolite_nm).(prm_nm).c;
                end
                
                %% perform the same but removing "outliers" (><mean*3SD)
                [~, ~, behavPrm_boxcox_clean] = rmv_outliers_3sd(behavPrm_boxcox);
                goodSubs_ter = ~isnan(metabolite_clean).*~isnan(con_data_clean).*~isnan(behavPrm_boxcox_clean) == 1;
                
                [mediation_path.boxcox.no_outliers.(MRS_ROI_nm).(metabolite_nm).(prm_nm).a,...
                    mediation_path.boxcox.no_outliers.(MRS_ROI_nm).(metabolite_nm).(prm_nm).b,...
                    mediation_path.boxcox.no_outliers.(MRS_ROI_nm).(metabolite_nm).(prm_nm).c,...
                    mediation_path.boxcox.no_outliers.(MRS_ROI_nm).(metabolite_nm).(prm_nm).c_prime,...
                    pval.boxcox.no_outliers.(MRS_ROI_nm).(metabolite_nm).(prm_nm),...
                    stats.boxcox.no_outliers.(MRS_ROI_nm).(metabolite_nm).(prm_nm)] = mediation(metabolite_allSubs(goodSubs_ter),...
                    con_data(goodSubs_ter),...
                    behavPrm_boxcox(goodSubs_ter),...
                    X_nm, M_nm, Y_nm, dispMed);
                % store how many subjects were kept
                N_goodSubs.boxcox.no_outliers.(MRS_ROI_nm).(metabolite_nm).(prm_nm) = sum(goodSubs_ter);
                
                % store when mediation is significant
                if pval.boxcox.no_outliers.(MRS_ROI_nm).(metabolite_nm).(prm_nm).a < 0.05 &&...
                        pval.boxcox.no_outliers.(MRS_ROI_nm).(metabolite_nm).(prm_nm).b < 0.05
                    pval.boxcox.no_outliers.signif.([MRS_ROI_nm,'_',metabolite_nm]).(prm_nm) = max(pval.boxcox.no_outliers.(MRS_ROI_nm).(metabolite_nm).(prm_nm).a,...
                        pval.boxcox.no_outliers.(MRS_ROI_nm).(metabolite_nm).(prm_nm).b);
                end
                
                % store when direct path is significant
                if pval.boxcox.no_outliers.(MRS_ROI_nm).(metabolite_nm).(prm_nm).c < 0.05
                    pval.direct_path.boxcox.no_outliers.signif.([MRS_ROI_nm,'_',metabolite_nm]).(prm_nm) = pval.boxcox.no_outliers.(MRS_ROI_nm).(metabolite_nm).(prm_nm).c;
                end
            end % filter bias parameter because bias is both negative and positive + already follows normal distribution
        end % parameter loop
    end % metabolites loop
end % ROI loop

%% line to launch to run again on a different contrast
% launch this before to avoid case where nothing is significant for the
% current BOLD contrast but the information from the previous test was kept
% clear('mediation_path','pval','N_goodSubs');

% %% lines to launch to display metabolite of interest
% MRS_ROI_nm = 'dmPFC';
% metabolite_nm = 'Glu_div_GSH';
% [metabolite_nm_bis] = metab_div_rnm(metabolite_nm);
% metabolite_allSubs = metabolites.(MRS_ROI_nm).(metabolite_nm);
% dispMed = 1;
% X_nm = [MRS_ROI_nm,'-',metabolite_nm_bis];
% M_nm = 'dmPFC=f(Ech)';
% 
% % kEp
% prm_nm='kEp';
% behavPrm = prm.(prm_nm);
% Y_nm = prm_nm;
% [~, ~, metabolite_clean] = rmv_outliers_3sd(metabolite_allSubs);
% [~, ~, con_data_clean] = rmv_outliers_3sd(con_data);
% [~, ~, behavPrm_clean] = rmv_outliers_3sd(behavPrm);
% goodSubs_bis = ~isnan(metabolite_clean).*~isnan(con_data_clean).*~isnan(behavPrm_clean) == 1;
% 
% mediation(metabolite_allSubs(goodSubs_bis),...
%     con_data(goodSubs_bis),...
%     behavPrm(goodSubs_bis),...
%     X_nm, M_nm, Y_nm, dispMed);
% 
% % kEm
% prm_nm='kEm';
% behavPrm = prm.(prm_nm);
% Y_nm = prm_nm;
% [~, ~, metabolite_clean] = rmv_outliers_3sd(metabolite_allSubs);
% [~, ~, con_data_clean] = rmv_outliers_3sd(con_data);
% [~, ~, behavPrm_clean] = rmv_outliers_3sd(behavPrm);
% goodSubs_bis = ~isnan(metabolite_clean).*~isnan(con_data_clean).*~isnan(behavPrm_clean) == 1;
% 
% mediation(metabolite_allSubs(goodSubs_bis),...
%     con_data(goodSubs_bis),...
%     behavPrm(goodSubs_bis),...
%     X_nm, M_nm, Y_nm, dispMed);
% 
% % boxcox kEp
% prm_nm = 'kEp';
% behavPrm_boxcox = (boxcox(prm.(prm_nm)'))';
% Y_nm = prm_nm;
% [~, ~, metabolite_clean] = rmv_outliers_3sd(metabolite_allSubs);
% [~, ~, con_data_clean] = rmv_outliers_3sd(con_data);
% [~, ~, behavPrm_clean] = rmv_outliers_3sd(behavPrm_boxcox);
% goodSubs_bis = ~isnan(metabolite_clean).*~isnan(con_data_clean).*~isnan(behavPrm_clean) == 1;
% 
% mediation(metabolite_allSubs(goodSubs_bis),...
%     con_data(goodSubs_bis),...
%     behavPrm_boxcox(goodSubs_bis),...
%     X_nm, M_nm, Y_nm, dispMed);
% 
% % boxcox kEm
% prm_nm='kEm';
% behavPrm_boxcox = (boxcox(prm.(prm_nm)'))';
% Y_nm = prm_nm;
% [~, ~, metabolite_clean] = rmv_outliers_3sd(metabolite_allSubs);
% [~, ~, con_data_clean] = rmv_outliers_3sd(con_data);
% [~, ~, behavPrm_clean] = rmv_outliers_3sd(behavPrm_boxcox);
% goodSubs_bis = ~isnan(metabolite_clean).*~isnan(con_data_clean).*~isnan(behavPrm_clean) == 1;
% 
% mediation(metabolite_allSubs(goodSubs_bis),...
%     con_data(goodSubs_bis),...
%     behavPrm_boxcox(goodSubs_bis),...
%     X_nm, M_nm, Y_nm, dispMed);

% disp(['dmPFC Glu/GSH => fMRI ',con_nm,'=> kEp (no outliers): p = ',...
%     num2str(max(pval.no_outliers.dmPFC.Glu_div_GSH.kEp.a,...
%     pval.no_outliers.dmPFC.Glu_div_GSH.kEp.b))]);
% disp(['dmPFC Glu/GSH => fMRI ',con_nm,'=> kEm (no outliers): p = ',...
%     num2str(max(pval.no_outliers.dmPFC.Glu_div_GSH.kEm.a,...
%     pval.no_outliers.dmPFC.Glu_div_GSH.kEm.b))]);

% kEp
disp(['kEp dmPFC-Glu/GSH: ',...
    'r(a) = ',num2str(round(stats.no_outliers.dmPFC.Glu_div_GSH.kEp.r.a,3)),...
    '; p(a) =  ',num2str(round(pval.no_outliers.dmPFC.Glu_div_GSH.kEp.a,3))]);
disp(['kEp dmPFC-Glu/GSH: ',...
    'r(b) = ',num2str(round(stats.no_outliers.dmPFC.Glu_div_GSH.kEp.r.b,3)),...
    '; p(b) =  ',num2str(round(pval.no_outliers.dmPFC.Glu_div_GSH.kEp.b,3))]);
disp(['kEp dmPFC-Glu/GSH: ',...
    'r(c) = ',num2str(round(stats.no_outliers.dmPFC.Glu_div_GSH.kEp.r.c,3)),...
    '; p(c) =  ',num2str(round(pval.no_outliers.dmPFC.Glu_div_GSH.kEp.c,3))]);
% kEm
disp(['kEm dmPFC-Glu/GSH: ',...
    'r(a) = ',num2str(round(stats.no_outliers.dmPFC.Glu_div_GSH.kEm.r.a,3)),...
    '; p(a) =  ',num2str(round(pval.no_outliers.dmPFC.Glu_div_GSH.kEm.a,3))]);
disp(['kEm dmPFC-Glu/GSH: ',...
    'r(b) = ',num2str(round(stats.no_outliers.dmPFC.Glu_div_GSH.kEm.r.b,3)),...
    '; p(b) =  ',num2str(round(pval.no_outliers.dmPFC.Glu_div_GSH.kEm.b,3))]);
disp(['kEm dmPFC-Glu/GSH: ',...
    'r(c) = ',num2str(round(stats.no_outliers.dmPFC.Glu_div_GSH.kEm.r.c,3)),...
    '; p(c) =  ',num2str(round(pval.no_outliers.dmPFC.Glu_div_GSH.kEm.c,3))]);