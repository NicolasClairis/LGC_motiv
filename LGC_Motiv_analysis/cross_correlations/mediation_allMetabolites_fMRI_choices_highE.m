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
% ROI_already_launched = 1;
ROI_already_launched_quest = questdlg('ROI already extracted?',...
    'ROI extracted?','Yes','No','No');
ROI_already_launched = strcmp(ROI_already_launched_quest,'Yes');
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
end
%% define regression estimate to look for in the fMRI GLM
con_idx = listdlg('promptstring','Which contrast?',...
    'ListString',con_names);
% con_nm = con_names{con_idx};
con_nm_str = inputdlg('Contrast short name?');
con_nm = con_nm_str{1};
con_data = NaN(1,NS);
con_data(:) = con_vec_all(con_idx, :, 1);

%% extract proportion of choices across individuals and tasks
fig_disp = 0;
[choice_hE] = choice_hE_proportion(study_nm, condition, subject_id, fig_disp);
task_names = {'Ep','Em','EpEm'};
nTasks = length(task_names);

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
        metabolite_nm_bis = metab_div_rnm(metabolite_nm);
        
        metabolite_allSubs = metabolites.(MRS_ROI_nm).(metabolite_nm);
        goodSubs = ~isnan(metabolite_allSubs);
        
        for iPrm = 1:nTasks
            prm_nm = task_names{iPrm};
            behavPrm = choice_hE.(prm_nm);
            
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
            behavPrm_boxcox = (boxcox(choice_hE.(prm_nm)'))';
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
        end % parameter loop
    end % metabolites loop
end % ROI loop

%% line to launch to run again on a different contrast
% launch this before to avoid case where nothing is significant for the
% current BOLD contrast but the information from the previous test was kept
% clear('mediation_path','pval','N_goodSubs');

% %% lines to launch to display metabolite of interest without outliers (but without boxcox transformation)
% MRS_ROI_nm='dmPFC';
% metabolite_nm='Glu_div_GSH';
% [metabolite_nm_bis] = metab_div_rnm(metabolite_nm);
% metabolite_allSubs = metabolites.(MRS_ROI_nm).(metabolite_nm);
% dispMed = 1;
% X_nm = [MRS_ROI_nm,'-',metabolite_nm_bis];
% M_nm='dmPFC=f(Ech)';
% 
% % Ep
% prm_nm='Ep';
% behavPrm = choice_hE.(prm_nm);
% Y_nm = ['choices ',prm_nm,' (%)'];
% 
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
% % Em
% prm_nm='Em';
% behavPrm = choice_hE.(prm_nm);
% Y_nm = ['choices ',prm_nm,' (%)'];
% 
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
% % pool Ep+Em
% prm_nm='EpEm';
% behavPrm = choice_hE.(prm_nm);
% Y_nm = ['choices ',prm_nm,' (%)'];
% 
% [~, ~, metabolite_clean] = rmv_outliers_3sd(metabolite_allSubs);
% [~, ~, con_data_clean] = rmv_outliers_3sd(con_data);
% [~, ~, behavPrm_clean] = rmv_outliers_3sd(behavPrm);
% goodSubs_bis = ~isnan(metabolite_clean).*~isnan(con_data_clean).*~isnan(behavPrm_clean) == 1;
% 
% mediation(metabolite_allSubs(goodSubs_bis),...
%     con_data(goodSubs_bis),...
%     behavPrm(goodSubs_bis),...
%     X_nm, M_nm, Y_nm, dispMed);

% disp(['dmPFC Glu/GSH => fMRI ',con_nm,'=> Ep+Em (no outliers): p = ',...
%     num2str(max(pval.no_outliers.dmPFC.Glu_div_GSH.EpEm.a,...
%     pval.no_outliers.dmPFC.Glu_div_GSH.EpEm.b))]);
% disp(['dmPFC Glu/GSH => fMRI ',con_nm,'=> Ep (no outliers): p = ',...
%     num2str(max(pval.no_outliers.dmPFC.Glu_div_GSH.Ep.a,...
%     pval.no_outliers.dmPFC.Glu_div_GSH.Ep.b))]);
% disp(['dmPFC Glu/GSH => fMRI ',con_nm,'=> Em (no outliers): p = ',...
%     num2str(max(pval.no_outliers.dmPFC.Glu_div_GSH.Em.a,...
%     pval.no_outliers.dmPFC.Glu_div_GSH.Em.b))]);

metab_nm = 'GSH'; % 'Glu_div_GSH'
metab_nm_bis = strrep(metab_nm,'_div_','/');
% Ep+Em
disp(['EpEm dmPFC-',metab_nm_bis,': ',...
    'a = ',num2str(round(stats.no_outliers.dmPFC.(metab_nm).EpEm.r.a,3)),...
    '; p(a) =  ',num2str(round(pval.no_outliers.dmPFC.(metab_nm).EpEm.a,3))]);
disp(['EpEm dmPFC-',metab_nm_bis,': ',...
    'b = ',num2str(round(stats.no_outliers.dmPFC.(metab_nm).EpEm.r.b,3)),...
    '; p(b) =  ',num2str(round(pval.no_outliers.dmPFC.(metab_nm).EpEm.b,3))]);
disp(['EpEm dmPFC-',metab_nm_bis,': ',...
    'c = ',num2str(round(stats.no_outliers.dmPFC.(metab_nm).EpEm.r.c,3)),...
    '; p(c) =  ',num2str(round(pval.no_outliers.dmPFC.(metab_nm).EpEm.c,3))]);

% Ep
disp(['Ep dmPFC-',metab_nm_bis,': ',...
    'a = ',num2str(round(stats.no_outliers.dmPFC.(metab_nm).Ep.r.a,3)),...
    '; p(a) =  ',num2str(round(pval.no_outliers.dmPFC.(metab_nm).Ep.a,3))]);
disp(['Ep dmPFC-',metab_nm_bis,': ',...
    'b = ',num2str(round(stats.no_outliers.dmPFC.(metab_nm).Ep.r.b,3)),...
    '; p(b) =  ',num2str(round(pval.no_outliers.dmPFC.(metab_nm).Ep.b,3))]);
disp(['Ep dmPFC-',metab_nm_bis,': ',...
    'c = ',num2str(round(stats.no_outliers.dmPFC.(metab_nm).Ep.r.c,3)),...
    '; p(c) =  ',num2str(round(pval.no_outliers.dmPFC.(metab_nm).Ep.c,3))]);

% Em
disp(['Em dmPFC-',metab_nm_bis,': ',...
    'r(a) = ',num2str(round(stats.no_outliers.dmPFC.(metab_nm).Em.r.a,3)),...
    '; p(a) =  ',num2str(round(pval.no_outliers.dmPFC.(metab_nm).Em.a,3))]);
disp(['Em dmPFC-',metab_nm_bis,': ',...
    'r(b) = ',num2str(round(stats.no_outliers.dmPFC.(metab_nm).Em.r.b,3)),...
    '; p(b) =  ',num2str(round(pval.no_outliers.dmPFC.(metab_nm).Em.b,3))]);
disp(['Em dmPFC-',metab_nm_bis,': ',...
    'r(c) = ',num2str(round(stats.no_outliers.dmPFC.(metab_nm).Em.r.c,3)),...
    '; p(c) =  ',num2str(round(pval.no_outliers.dmPFC.(metab_nm).Em.c,3))]);


%% save results automatically
savePath = fullfile('P:','boulot','postdoc_CarmenSandi','results','ROI');
file_nm = [savePath,filesep,'GLM',num2str(GLM),'_',...
    ROI_coords.ROI_nm.ROI_1_shortName,'_ROI_',num2str(NS),'subs.mat'];
if ~exist(file_nm,'file')
    save(file_nm);
end