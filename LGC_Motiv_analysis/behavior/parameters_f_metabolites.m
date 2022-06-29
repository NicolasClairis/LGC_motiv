%% script to correlate brain metabolites with behavioral scores
%
% designed by N.Clairis - 2022

%% define metabolite and ROI you want to focus on
% ROI
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

%% define all subjects
condition = subject_condition();
[subject_id, NS] = LGCM_subject_selection('study1', condition);

%% extract all metabolites
[metabolites] = metabolite_load(subject_id);
% focus on metabolite and brain area selected
metabolite_allSubs = metabolites.(ROI_nm).(metabolite_nm);

%% extract behavioral parameters
listPossibleConditions = {'bayesian','simple'};
mdlType_idx = listdlg('promptstring','Which model type?','ListString',listPossibleConditions);
mdlType = listPossibleConditions{mdlType_idx};


switch mdlType
    case 'bayesian'
        bayesian_mdl = getfield(load('behavioral_prm_tmp.mat','bayesian_mdl3'),...
            'bayesian_mdl3');
        [prm.kR, prm.kP,...
            prm.kEp, prm.kFp,...
            prm.kEm, prm.kFm] = deal(NaN(1,NS));
        parameters = {'kR','kP',...
            'kEp','kFp',...
            'kEm','kFm'};
        for iS = 1:NS
            sub_nm = subject_id{iS};
            % extract parameters
            sub_idx = strcmp(bayesian_mdl.subject_id, sub_nm);
            if sum(sub_idx == 1)
                prm.kR(iS) = bayesian_mdl.kR(sub_idx);
                prm.kP(iS) = bayesian_mdl.kP(sub_idx);
                prm.kEp(iS) = bayesian_mdl.kEp(sub_idx);
                prm.kFp(iS) = bayesian_mdl.kFp(sub_idx);
                prm.kEm(iS) = bayesian_mdl.kEm(sub_idx);
                prm.kFm(iS) = bayesian_mdl.kFm(sub_idx);
            end % filter if subject extracted by Arthur
        end % subject list
    case 'simple'
        %% perform behavioral model
        figDispGroup = 0;
        dispMoneyOrLevels = 'levels';
        [betas_fullList, pvalues_fullList] = logitfit_choices_group(figDispGroup, dispMoneyOrLevels);

        [prm.kRp, prm.kPp, prm.kEp, prm.kFp,...
            prm.kRm, prm.kPm, prm.kEm, prm.kFm] = deal(NaN(1,NS));
        parameters = {'kRp','kPp','kEp','kFp',...
            'kRm','kPm','kEm','kFm'};
        mdl_nm = 'mdl_3';
        for iS = 1:NS
            sub_nm = subject_id{iS};
            % extract physical task parameters
            sub_idx = strcmp(betas_fullList.subList,sub_nm);
            prm.kRp(iS) = betas_fullList.Ep.(mdl_nm).kR(sub_idx);
            prm.kPp(iS) = betas_fullList.Ep.(mdl_nm).kP(sub_idx);
            prm.kEp(iS) = betas_fullList.Ep.(mdl_nm).kEffort(sub_idx);
            prm.kFp(iS) = betas_fullList.Ep.(mdl_nm).kFatigue(sub_idx);
            % extract mental task parameters
            prm.kRm(iS) = betas_fullList.Em.(mdl_nm).kR(sub_idx);
            prm.kPm(iS) = betas_fullList.Em.(mdl_nm).kP(sub_idx);
            prm.kEm(iS) = betas_fullList.Em.(mdl_nm).kEffort(sub_idx);
            prm.kFm(iS) = betas_fullList.Em.(mdl_nm).kFatigue(sub_idx);
        end % subject list
end % bayesian/behavioral

%% perform GLM and correlation between behavioral parameters and spectroscopy
n_prm = length(parameters);
% remove NaN subjects
for iPrm = 1:n_prm
    prm_nm = parameters{iPrm};
    goodSubs = ~isnan(metabolite_allSubs).*~isnan(prm.(prm_nm));
    var_nm = [prm_nm,'_f_',metabolite_nm];
    [betas.(var_nm), ~,stats_tmp] = glmfit(metabolite_allSubs(goodSubs), prm.(prm_nm)(goodSubs),'normal');
    pval.(var_nm) = stats_tmp.p;
    prm_f_metabolite_fit.(var_nm) = glmval(betas.(var_nm),metabolite_allSubs(goodSubs),'identity');
end % parameters loop

%% display results
pSize = 50;
for iPrm = 1:n_prm
    prm_nm = parameters{iPrm};
    var_nm = [prm_nm,'_f_',metabolite_nm];
    fig;
    hdl=scatter(metabolite_allSubs, prm.(prm_nm),...
        'MarkerEdgeColor','k','SizeData',100,'LineWidth',3);
    hold on;
    [metabolite_allSubs_sorted, idx_sorted] = sort(metabolite_allSubs(goodSubs));
    plot(metabolite_allSubs_sorted, prm_f_metabolite_fit.(var_nm)(idx_sorted),...
        'LineStyle','--','LineWidth',3);
    xlabel([ROI_nm,' ',metabolite_nm]);
    ylabel(prm_nm);
    legend_size(pSize);

end % parameter loop