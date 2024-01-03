% script to look at the mediation going from blood lactate to fMRI to
% motivated behavior (1) PHE/MHE and 2) behavioral parameters)

%% did you already launch the ROI extraction (1) or not (0)?
% ROI_already_launched = 1;
ROI_already_launched_quest = questdlg('ROI already extracted in the workspace?',...
    'ROI extracted in the workspace?','Yes','No','No');
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
    
    %% load blood lactate
    plasma_Lac_struct = load_plasma_Lac(subject_id);
    plasma_Lac = plasma_Lac_struct.Lac;
    
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

%% selection of participants
goodSubs = ~isnan(plasma_Lac);

%% mediation variables
X_nm = 'plasma_Lac';
M_nm = ['fMRI-',ROI_coords.ROI_nm.ROI_1_shortName,'-',con_nm];

%% choices HPE/HME/HE
for iTasks = 1:nTasks
    task_nm = task_names{iTasks};
    choicePrm = choice_hE.(task_nm);
    
    Y_nm = task_nm;
    [mediation_path.plasma_Lac.choices.(task_nm).a,...
        mediation_path.plasma_Lac.choices.(task_nm).b,...
        mediation_path.plasma_Lac.choices.(task_nm).c,...
        mediation_path.plasma_Lac.choices.(task_nm).c_prime,...
        pval.plasma_Lac.choices.(task_nm),...
        stats.plasma_Lac.choices.(task_nm)] = mediation(plasma_Lac(goodSubs),...
        con_data(goodSubs),...
        choicePrm(goodSubs),...
        X_nm, M_nm, Y_nm, dispMed);
    % store how many subjects were kept
    N_goodSubs.plasma_Lac.choices = sum(goodSubs);
    
    % store when mediation is significant
    if pval.plasma_Lac.choices.(task_nm).a < 0.05 &&...
            pval.plasma_Lac.choices.(task_nm).b < 0.05
        pval.signif.plasma_Lac.choices.(task_nm) = max(pval.plasma_Lac.choices.(task_nm).a,...
            pval.plasma_Lac.choices.(task_nm).b);
    end
    
    % store when direct path is significant
    if pval.plasma_Lac.choices.(task_nm).c < 0.05
        pval.direct_path.signif.plasma_Lac.choices.(task_nm) = pval.plasma_Lac.choices.(task_nm).c;
    end
    
    %% perform the same but removing "outliers" (><mean*3SD)
    [~, ~, plasma_Lac_clean] = rmv_outliers_3sd(plasma_Lac);
    [~, ~, con_data_clean] = rmv_outliers_3sd(con_data);
    [~, ~, behavPrm_clean] = rmv_outliers_3sd(choicePrm);
    goodSubs_bis_choices = ~isnan(plasma_Lac_clean).*~isnan(con_data_clean).*~isnan(behavPrm_clean) == 1;
    
    [mediation_path.no_outliers.plasma_Lac.choices.(task_nm).a,...
        mediation_path.no_outliers.plasma_Lac.choices.(task_nm).b,...
        mediation_path.no_outliers.plasma_Lac.choices.(task_nm).c,...
        mediation_path.no_outliers.plasma_Lac.choices.(task_nm).c_prime,...
        pval.no_outliers.plasma_Lac.choices.(task_nm),...
        stats.no_outliers.plasma_Lac.choices.(task_nm)] = mediation(plasma_Lac(goodSubs_bis_choices),...
        con_data(goodSubs_bis_choices),...
        choicePrm(goodSubs_bis_choices),...
        X_nm, M_nm, Y_nm, dispMed);
    % store how many subjects were kept
    N_goodSubs.no_outliers.plasma_Lac.choices.(task_nm) = sum(goodSubs_bis_choices);
    
    % store when mediation is significant
    if pval.no_outliers.plasma_Lac.choices.(task_nm).a < 0.05 &&...
            pval.no_outliers.plasma_Lac.choices.(task_nm).b < 0.05
        pval.no_outliers.signif.plasma_Lac.choices.(task_nm) = max(pval.no_outliers.plasma_Lac.choices.(task_nm).a,...
            pval.no_outliers.plasma_Lac.choices.(task_nm).b);
    end
    
    % store when direct path is significant
    if pval.no_outliers.plasma_Lac.choices.(task_nm).c < 0.05
        pval.direct_path.no_outliers.signif.plasma_Lac.choices.(task_nm) = pval.no_outliers.plasma_Lac.choices.(task_nm).c;
    end
    
    %% same but with boxcox transformation of behavioral parameters
    behavPrm_boxcox = (boxcox(choice_hE.(task_nm)'))';
    [mediation_path.boxcox.plasma_Lac.choices.(task_nm).a,...
        mediation_path.boxcox.plasma_Lac.choices.(task_nm).b,...
        mediation_path.boxcox.plasma_Lac.choices.(task_nm).c,...
        mediation_path.boxcox.plasma_Lac.choices.(task_nm).c_prime,...
        pval.boxcox.plasma_Lac.choices.(task_nm),...
        stats.boxcox.plasma_Lac.choices.(task_nm)] = mediation(plasma_Lac(goodSubs),...
        con_data(goodSubs),...
        behavPrm_boxcox(goodSubs),...
        X_nm, M_nm, Y_nm, dispMed);
    % store how many subjects were kept
    N_goodSubs.boxcox.plasma_Lac.choices = sum(goodSubs);
    
    % store when mediation is significant
    if pval.boxcox.plasma_Lac.choices.(task_nm).a < 0.05 &&...
            pval.boxcox.plasma_Lac.choices.(task_nm).b < 0.05
        pval.boxcox.signif.plasma_Lac.choices.(task_nm) = max(pval.boxcox.plasma_Lac.choices.(task_nm).a,...
            pval.boxcox.plasma_Lac.choices.(task_nm).b);
    end
    
    % store when direct path is significant
    if pval.boxcox.plasma_Lac.choices.(task_nm).c < 0.05
        pval.direct_path.boxcox.signif.plasma_Lac.choices.(task_nm) = pval.boxcox.plasma_Lac.choices.(task_nm).c;
    end
    
    %% perform the same but removing "outliers" (><mean*3SD)
    [~, ~, behavPrm_boxcox_clean] = rmv_outliers_3sd(behavPrm_boxcox);
    goodSubs_ter_choices = ~isnan(plasma_Lac_clean).*~isnan(con_data_clean).*~isnan(behavPrm_boxcox_clean) == 1;
    
    [mediation_path.boxcox.no_outliers.plasma_Lac.choices.(task_nm).a,...
        mediation_path.boxcox.no_outliers.plasma_Lac.choices.(task_nm).b,...
        mediation_path.boxcox.no_outliers.plasma_Lac.choices.(task_nm).c,...
        mediation_path.boxcox.no_outliers.plasma_Lac.choices.(task_nm).c_prime,...
        pval.boxcox.no_outliers.plasma_Lac.choices.(task_nm),...
        stats.boxcox.no_outliers.plasma_Lac.choices.(task_nm)] = mediation(plasma_Lac(goodSubs_ter_choices),...
        con_data(goodSubs_ter_choices),...
        behavPrm_boxcox(goodSubs_ter_choices),...
        X_nm, M_nm, Y_nm, dispMed);
    % store how many subjects were kept
    N_goodSubs.boxcox.no_outliers.plasma_Lac.choices.(task_nm) = sum(goodSubs_ter_choices);
    
    % store when mediation is significant
    if pval.boxcox.no_outliers.plasma_Lac.choices.(task_nm).a < 0.05 &&...
            pval.boxcox.no_outliers.plasma_Lac.choices.(task_nm).b < 0.05
        pval.boxcox.no_outliers.signif.plasma_Lac.choices.(task_nm) = max(pval.boxcox.no_outliers.plasma_Lac.choices.(task_nm).a,...
            pval.boxcox.no_outliers.plasma_Lac.choices.(task_nm).b);
    end
    
    % store when direct path is significant
    if pval.boxcox.no_outliers.plasma_Lac.choices.(task_nm).c < 0.05
        pval.direct_path.boxcox.no_outliers.signif.plasma_Lac.choices.(task_nm) = pval.boxcox.no_outliers.plasma_Lac.choices.(task_nm).c;
    end
end % task loop for choices

%% behavioral parameters
for iPrm = 1:nPrm
    prm_nm = parameter_names{iPrm};
    behavPrm = prm.(prm_nm);
    
    Y_nm = prm_nm;
    [mediation_path.plasma_Lac.prm.(prm_nm).a,...
        mediation_path.plasma_Lac.prm.(prm_nm).b,...
        mediation_path.plasma_Lac.prm.(prm_nm).c,...
        mediation_path.plasma_Lac.prm.(prm_nm).c_prime,...
        pval.plasma_Lac.prm.(prm_nm),...
        stats.plasma_Lac.prm.(prm_nm)] = mediation(plasma_Lac(goodSubs),...
        con_data(goodSubs),...
        behavPrm(goodSubs),...
        X_nm, M_nm, Y_nm, dispMed);
    % store how many subjects were kept
    N_goodSubs.plasma_Lac.prm = sum(goodSubs);
    
    % store when mediation is significant
    if pval.plasma_Lac.prm.(prm_nm).a < 0.05 &&...
            pval.plasma_Lac.prm.(prm_nm).b < 0.05
        pval.signif.plasma_Lac.prm.(prm_nm) = max(pval.plasma_Lac.prm.(prm_nm).a,...
            pval.plasma_Lac.prm.(prm_nm).b);
    end
    
    % store when direct path is significant
    if pval.plasma_Lac.prm.(prm_nm).c < 0.05
        pval.direct_path.signif.plasma_Lac.prm.(prm_nm) = pval.plasma_Lac.prm.(prm_nm).c;
    end
    
    %% perform the same but removing "outliers" (><mean*3SD)
    [~, ~, plasma_Lac_clean] = rmv_outliers_3sd(plasma_Lac);
    [~, ~, con_data_clean] = rmv_outliers_3sd(con_data);
    [~, ~, behavPrm_clean] = rmv_outliers_3sd(behavPrm);
    goodSubs_bis_prm = ~isnan(plasma_Lac_clean).*~isnan(con_data_clean).*~isnan(behavPrm_clean) == 1;
    
    [mediation_path.no_outliers.plasma_Lac.prm.(prm_nm).a,...
        mediation_path.no_outliers.plasma_Lac.prm.(prm_nm).b,...
        mediation_path.no_outliers.plasma_Lac.prm.(prm_nm).c,...
        mediation_path.no_outliers.plasma_Lac.prm.(prm_nm).c_prime,...
        pval.no_outliers.plasma_Lac.prm.(prm_nm),...
        stats.no_outliers.plasma_Lac.prm.(prm_nm)] = mediation(plasma_Lac(goodSubs_bis_prm),...
        con_data(goodSubs_bis_prm),...
        behavPrm(goodSubs_bis_prm),...
        X_nm, M_nm, Y_nm, dispMed);
    % store how many subjects were kept
    N_goodSubs.no_outliers.plasma_Lac.prm.(prm_nm) = sum(goodSubs_bis_prm);
    
    % store when mediation is significant
    if pval.no_outliers.plasma_Lac.prm.(prm_nm).a < 0.05 &&...
            pval.no_outliers.plasma_Lac.prm.(prm_nm).b < 0.05
        pval.no_outliers.signif.plasma_Lac.prm.(prm_nm) = max(pval.no_outliers.plasma_Lac.prm.(prm_nm).a,...
            pval.no_outliers.plasma_Lac.prm.(prm_nm).b);
    end
    
    % store when direct path is significant
    if pval.no_outliers.plasma_Lac.prm.(prm_nm).c < 0.05
        pval.direct_path.no_outliers.signif.plasma_Lac.prm.(prm_nm) = pval.no_outliers.plasma_Lac.prm.(prm_nm).c;
    end
    
    %% same but with boxcox transformation of behavioral parameters
    if ~strcmp(prm_nm,'kBiasM')
        behavPrm_boxcox = (boxcox(prm.(prm_nm)'))';
        [mediation_path.boxcox.plasma_Lac.prm.(prm_nm).a,...
            mediation_path.boxcox.plasma_Lac.prm.(prm_nm).b,...
            mediation_path.boxcox.plasma_Lac.prm.(prm_nm).c,...
            mediation_path.boxcox.plasma_Lac.prm.(prm_nm).c_prime,...
            pval.boxcox.plasma_Lac.prm.(prm_nm),...
            stats.boxcox.plasma_Lac.prm.(prm_nm)] = mediation(plasma_Lac(goodSubs),...
            con_data(goodSubs),...
            behavPrm_boxcox(goodSubs),...
            X_nm, M_nm, Y_nm, dispMed);
        % store how many subjects were kept
        N_goodSubs.boxcox.plasma_Lac.prm = sum(goodSubs);
        
        % store when mediation is significant
        if pval.boxcox.plasma_Lac.prm.(prm_nm).a < 0.05 &&...
                pval.boxcox.plasma_Lac.prm.(prm_nm).b < 0.05
            pval.boxcox.signif.plasma_Lac.prm.(prm_nm) = max(pval.boxcox.plasma_Lac.prm.(prm_nm).a,...
                pval.boxcox.plasma_Lac.prm.(prm_nm).b);
        end
        
        % store when direct path is significant
        if pval.boxcox.plasma_Lac.prm.(prm_nm).c < 0.05
            pval.direct_path.boxcox.signif.plasma_Lac.prm.(prm_nm) = pval.boxcox.plasma_Lac.prm.(prm_nm).c;
        end
        
        %% perform the same but removing "outliers" (><mean*3SD)
        [~, ~, behavPrm_boxcox_clean] = rmv_outliers_3sd(behavPrm_boxcox);
        goodSubs_ter_prm = ~isnan(plasma_Lac_clean).*~isnan(con_data_clean).*~isnan(behavPrm_boxcox_clean) == 1;
        
        [mediation_path.boxcox.no_outliers.plasma_Lac.prm.(prm_nm).a,...
            mediation_path.boxcox.no_outliers.plasma_Lac.prm.(prm_nm).b,...
            mediation_path.boxcox.no_outliers.plasma_Lac.prm.(prm_nm).c,...
            mediation_path.boxcox.no_outliers.plasma_Lac.prm.(prm_nm).c_prime,...
            pval.boxcox.no_outliers.plasma_Lac.prm.(prm_nm),...
            stats.boxcox.no_outliers.plasma_Lac.prm.(prm_nm)] = mediation(plasma_Lac(goodSubs_ter_prm),...
            con_data(goodSubs_ter_prm),...
            behavPrm_boxcox(goodSubs_ter_prm),...
            X_nm, M_nm, Y_nm, dispMed);
        % store how many subjects were kept
        N_goodSubs.boxcox.no_outliers.plasma_Lac.prm.(prm_nm) = sum(goodSubs_ter_prm);
        
        % store when mediation is significant
        if pval.boxcox.no_outliers.plasma_Lac.prm.(prm_nm).a < 0.05 &&...
                pval.boxcox.no_outliers.plasma_Lac.prm.(prm_nm).b < 0.05
            pval.boxcox.no_outliers.signif.plasma_Lac.prm.(prm_nm) = max(pval.boxcox.no_outliers.plasma_Lac.prm.(prm_nm).a,...
                pval.boxcox.no_outliers.plasma_Lac.prm.(prm_nm).b);
        end
        
        % store when direct path is significant
        if pval.boxcox.no_outliers.plasma_Lac.prm.(prm_nm).c < 0.05
            pval.direct_path.boxcox.no_outliers.signif.plasma_Lac.prm.(prm_nm) = pval.boxcox.no_outliers.plasma_Lac.prm.(prm_nm).c;
        end
    end % filter bias parameter because bias is both negative and positive + already follows normal distribution
end % parameter loop

%% display main results
%% choices Ep+Em
disp(['EpEm plasma Lac: ',...
    'c = ',num2str(round(stats.no_outliers.plasma_Lac.choices.EpEm.r.c,3)),...
    '; p(c) =  ',num2str(round(pval.no_outliers.plasma_Lac.choices.EpEm.c,3))]);
disp(['EpEm plasma Lac: ',...
    'a = ',num2str(round(stats.no_outliers.plasma_Lac.choices.EpEm.r.a,3)),...
    '; p(a) =  ',num2str(round(pval.no_outliers.plasma_Lac.choices.EpEm.a,3))]);
disp(['EpEm plasma Lac: ',...
    'b = ',num2str(round(stats.no_outliers.plasma_Lac.choices.EpEm.r.b,3)),...
    '; p(b) =  ',num2str(round(pval.no_outliers.plasma_Lac.choices.EpEm.b,3))]);
disp(['EpEm plasma Lac: ',...
    'c'' = ',num2str(round(stats.no_outliers.plasma_Lac.choices.EpEm.r.c_prime,3)),...
    '; p(c'') =  ',num2str(round(pval.no_outliers.plasma_Lac.choices.EpEm.c_prime,3))]);

%% choices Ep
disp(['Ep plasma Lac: ',...
    'c = ',num2str(round(stats.no_outliers.plasma_Lac.choices.Ep.r.c,3)),...
    '; p(c) =  ',num2str(round(pval.no_outliers.plasma_Lac.choices.Ep.c,3))]);
disp(['Ep plasma Lac: ',...
    'a = ',num2str(round(stats.no_outliers.plasma_Lac.choices.Ep.r.a,3)),...
    '; p(a) =  ',num2str(round(pval.no_outliers.plasma_Lac.choices.Ep.a,3))]);
disp(['Ep plasma Lac: ',...
    'b = ',num2str(round(stats.no_outliers.plasma_Lac.choices.Ep.r.b,3)),...
    '; p(b) =  ',num2str(round(pval.no_outliers.plasma_Lac.choices.Ep.b,3))]);
disp(['Ep plasma Lac: ',...
    'c'' = ',num2str(round(stats.no_outliers.plasma_Lac.choices.Ep.r.c_prime,3)),...
    '; p(c'') =  ',num2str(round(pval.no_outliers.plasma_Lac.choices.Ep.c_prime,3))]);

%% choices Em
disp(['Em plasma Lac: ',...
    'r(c) = ',num2str(round(stats.no_outliers.plasma_Lac.choices.Em.r.c,3)),...
    '; p(c) =  ',num2str(round(pval.no_outliers.plasma_Lac.choices.Em.c,3))]);
disp(['Em plasma Lac: ',...
    'r(a) = ',num2str(round(stats.no_outliers.plasma_Lac.choices.Em.r.a,3)),...
    '; p(a) =  ',num2str(round(pval.no_outliers.plasma_Lac.choices.Em.a,3))]);
disp(['Em plasma Lac: ',...
    'r(b) = ',num2str(round(stats.no_outliers.plasma_Lac.choices.Em.r.b,3)),...
    '; p(b) =  ',num2str(round(pval.no_outliers.plasma_Lac.choices.Em.b,3))]);
disp(['Em plasma Lac: ',...
    'r(c'') = ',num2str(round(stats.no_outliers.plasma_Lac.choices.Em.r.c_prime,3)),...
    '; p(c'') =  ',num2str(round(pval.no_outliers.plasma_Lac.choices.Em.c_prime,3))]);

%% kEp
disp(['kEp plasma Lac: ',...
    'r(c) = ',num2str(round(stats.no_outliers.plasma_Lac.prm.kEp.r.c,3)),...
    '; p(c) =  ',num2str(round(pval.no_outliers.plasma_Lac.prm.kEp.c,3))]);
disp(['kEp plasma Lac: ',...
    'r(a) = ',num2str(round(stats.no_outliers.plasma_Lac.prm.kEp.r.a,3)),...
    '; p(a) =  ',num2str(round(pval.no_outliers.plasma_Lac.prm.kEp.a,3))]);
disp(['kEp plasma Lac: ',...
    'r(b) = ',num2str(round(stats.no_outliers.plasma_Lac.prm.kEp.r.b,3)),...
    '; p(b) =  ',num2str(round(pval.no_outliers.plasma_Lac.prm.kEp.b,3))]);
disp(['kEp plasma Lac: ',...
    'r(c'') = ',num2str(round(stats.no_outliers.plasma_Lac.prm.kEp.r.c_prime,3)),...
    '; p(c'') =  ',num2str(round(pval.no_outliers.plasma_Lac.prm.kEp.c_prime,3))]);

%% kEm
disp(['kEm plasma Lac: ',...
    'r(c) = ',num2str(round(stats.no_outliers.plasma_Lac.prm.kEm.r.c,3)),...
    '; p(c) =  ',num2str(round(pval.no_outliers.plasma_Lac.prm.kEm.c,3))]);
disp(['kEm plasma Lac: ',...
    'r(a) = ',num2str(round(stats.no_outliers.plasma_Lac.prm.kEm.r.a,3)),...
    '; p(a) =  ',num2str(round(pval.no_outliers.plasma_Lac.prm.kEm.a,3))]);
disp(['kEm plasma Lac: ',...
    'r(b) = ',num2str(round(stats.no_outliers.plasma_Lac.prm.kEm.r.b,3)),...
    '; p(b) =  ',num2str(round(pval.no_outliers.plasma_Lac.prm.kEm.b,3))]);
disp(['kEm plasma Lac: ',...
    'r(c'') = ',num2str(round(stats.no_outliers.plasma_Lac.prm.kEm.r.c_prime,3)),...
    '; p(c'') =  ',num2str(round(pval.no_outliers.plasma_Lac.prm.kEm.c_prime,3))]);

%% save results automatically
savePath = fullfile('P:','boulot','postdoc_CarmenSandi','results','mediation','plasma_Lac');
file_nm = [savePath,filesep,'GLM',num2str(GLM),'_',...
    ROI_coords.ROI_nm.ROI_1_shortName,'_ROI_',num2str(NS),'subs.mat'];
if ~exist(file_nm,'file')
    save(file_nm);
end



