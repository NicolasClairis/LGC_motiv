% script to correlate behavioral parameters from the model to the blood
% variables related to NAD.

%% working directory
pcRoot = LGCM_root_paths;
switch pcRoot
    case 'E:\'
        gitPath = fullfile('C:','Users','clairis','Desktop','GitHub',...
            'LGC_motiv');
    case 'L:\human_data_private\raw_data_subject\'
        gitPath = fullfile('C:','Users','Loco','Documents','GitHub',...
            'LGC_motiv');
end

%% select study and participants
study_nm = 'study1';
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm,condition);
dataResultsPath = fullfile(gitPath, 'LGC_Motiv_results',study_nm);

%% extract all blood data
[bloodTable, blood_NAD_sub_List] = load_blood_NAD(study_nm);
blood_metabolites = fieldnames(bloodTable);
nNADHmetab = length(blood_metabolites);

%% load behavioral parameters
[prm_all, mdlType, mdlN] = prm_extraction(study_nm, subject_id);
parameters = fieldnames(prm_all);
n_prm = length(parameters);
for iPrm = 1:n_prm
    prm_nm = parameters{iPrm};
    if ~strcmp(prm_nm,'CID')
        prm.(prm_nm) = NaN(1,NS);
    end
end

%% loop through subjects
for iBloodM = 1:nNADHmetab
    blood.(blood_metabolites{iBloodM}) = NaN(1,NS);
end
for iS = 1:NS
    sub_nm = subject_id{iS};
    sub_nm_bis = ['CID',sub_nm];
    % identify index for the current subject on the different files
    sub_blood_idx = find(strcmp(blood_NAD_sub_List, sub_nm_bis));
    %% extract blood values
    for iBloodM = 1:nNADHmetab
        blood_m_nm = blood_metabolites{iBloodM};
        if ~isempty(sub_blood_idx)
            blood.(blood_m_nm)(iS) = bloodTable.(blood_m_nm)(sub_blood_idx);
        end
    end % metabolite loop
    %% load parameters
    jS = strcmp(sub_nm,prm_all.CID);
    for iPrm = 1:n_prm
        prm_nm = parameters{iPrm};
        if ~strcmp(prm_nm,'CID')
            prm.(prm_nm)(iS) = prm_all.(prm_nm)(jS);
        end
    end
end % subject loop

%% test the correlations
for iPrm = 1:n_prm
    prm_nm = parameters{iPrm};
    if ~strcmp(prm_nm,'CID')
        for iBlood = 1:nNADHmetab
            blood_m_nm = blood_metabolites{iBlood};
            % filter bad subjects
            subs_ok = (~isnan(blood.(blood_m_nm)).*~isnan(prm.(prm_nm))) == 1;
            if sum(subs_ok) > 0
                % perform glm brain GSH = b0 + b1*nutrition variable
                [betas.(prm_nm).(['f_',blood_m_nm]),~,stats_tmp] = glmfit(blood.(blood_m_nm)(subs_ok), prm.(prm_nm)(subs_ok),'normal');
                % extract p.value
                pval.(prm_nm).(['f_',blood_m_nm]) = stats_tmp.p;
                % note which are significant
                if stats_tmp.p(2) < 0.05
                    pval.signif.([prm_nm,'_f_',blood_m_nm]) = stats_tmp.p(2);
                end
                % extract the fit
                blood_bis.(prm_nm).(blood_m_nm) = sort(blood.(blood_m_nm)(subs_ok));
                fittedData.(prm_nm).(['f_',blood_m_nm]) = glmval(betas.(prm_nm).(['f_',blood_m_nm]), blood_bis.(prm_nm).(blood_m_nm), 'identity');
            end
        end % nutrition component
    end % CID field filter
end % parameter

%% figure
figDisp = 1;
if figDisp == 1
    
    lWidth = 3;
    pSize = 25;
    rawDataCol = 'k';
    fitCol = [255 143 143]./255;
    
    %% show results
    for iPrm = 1:n_prm
        prm_nm = parameters{iPrm};
        
        if ~strcmp(prm_nm,'CID')
            % initiliaze figures
            % one figure for all measures and one for integrated measures
            fig1 = fig; j_fig1 = 0;
            fig2 = fig; j_fig2 = 0;
            for iBlood = 1:nNADHmetab
                blood_m_nm = blood_metabolites{iBlood};
                %% raw food
                switch blood_m_nm
                    case {'Nam','NMN','NR','NAD','NADH','NADP','NADPH','MeNam','MeXPY'}
                        figure(fig1);
                        j_fig1 = j_fig1 + 1;
                        subplot(3,4,j_fig1);
                    case {'NAD_div_NADH','NADP_div_NADPH',...
                            'total_NAD_precursors','total_NAD',...
                            'total_NAD_with_precursors',...
                            'total_NAD_with_byproducts','total_NAD_byproducts'}
                        figure(fig2);
                        j_fig2 = j_fig2 + 1;
                        subplot(3,3,j_fig2);
                    otherwise
                        error([blood_m_nm,' not in the list']);
                end
                hold on;
                % show raw data
                curve_hdl = scatter(blood.(blood_m_nm), prm.(prm_nm));
                curve_hdl.LineWidth = lWidth;
                curve_hdl.MarkerEdgeColor = rawDataCol;
                % add the fit
                curve_fit_hdl = plot(blood_bis.(prm_nm).(blood_m_nm),...
                    fittedData.(prm_nm).(['f_',(blood_m_nm)]));
                curve_fit_hdl.LineStyle = '--';
                curve_fit_hdl.LineWidth = lWidth;
                curve_fit_hdl.Color = fitCol;
                xlabel(blood_m_nm);
                ylabel(prm_nm);
                legend_size(pSize);
                
            end % whole-blood metabolite
        end % CID field filter
    end % parameter loop
end % figure display