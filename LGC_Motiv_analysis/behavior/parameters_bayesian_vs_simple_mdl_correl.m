%% check correlation between parameters extracted with simple frequentist
% model vs bayesian modelling approach

%% define all subjects
study_nm = 'study1';
condition = subject_condition();
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%%
bayesian_root = fullfile('C:','Users','clairis','Desktop','GitHub',...
    'LGC_motiv','LGC_Motiv_analysis','behavior');

%% extract behavioral parameters
%% extract bayesian model
bayesian_mdl = getfield(load([bayesian_root,filesep,...
    'behavioral_prm_tmp.mat'],'bayesian_mdl3'),...
    'bayesian_mdl3');
[bayesian_prm.kR, bayesian_prm.kP,...
    bayesian_prm.kEp, bayesian_prm.kFp,...
    bayesian_prm.kEm, bayesian_prm.kFm] = deal(NaN(1,NS));
bayesian_parameters = fieldnames(bayesian_prm);
for iS = 1:NS
    sub_nm = subject_id{iS};
    % extract parameters
    sub_idx = strcmp(bayesian_mdl.subject_id, sub_nm);
    if sum(sub_idx == 1)
        for iBPrm = 1:length(bayesian_parameters)
            bayesian_prm_nm = bayesian_parameters{iBPrm};
            bayesian_prm.(bayesian_prm_nm)(iS) = bayesian_mdl.(bayesian_prm_nm)(sub_idx);
        end
    end % filter if subject extracted by Arthur
end % subject list

%% perform behavioral model
figDispGroup = 0;
dispMoneyOrLevels = 'levels';
[betas_fullList, pvalues_fullList] = logitfit_choices_group(figDispGroup, dispMoneyOrLevels);

mdl_nm = 'mdl_3';
switch mdl_nm
    case 'mdl_3'
        [simpleMdl_prm.kMp, simpleMdl_prm.kEp, simpleMdl_prm.kFp,...
            simpleMdl_prm.kMm, simpleMdl_prm.kEm, simpleMdl_prm.kFm] = deal(NaN(1,NS));
    case 'mdl_4'
        [simpleMdl_prm.kRp, simpleMdl_prm.kPp, simpleMdl_prm.kEp, simpleMdl_prm.kFp,...
            simpleMdl_prm.kRm, simpleMdl_prm.kPm, simpleMdl_prm.kEm, simpleMdl_prm.kFm] = deal(NaN(1,NS));
end
simple_parameters = fieldnames(simpleMdl_prm);
for iS = 1:NS
    sub_nm = subject_id{iS};
    % extract physical task parameters
    sub_idx = strcmp(betas_fullList.subList,sub_nm);
    switch mdl_nm
        case 'mdl_3'
            simpleMdl_prm.kMp(iS) = betas_fullList.Ep.(mdl_nm).kMoney(sub_idx);
            simpleMdl_prm.kMm(iS) = betas_fullList.Em.(mdl_nm).kMoney(sub_idx);
        case 'mdl_4'
            simpleMdl_prm.kRp(iS) = betas_fullList.Ep.(mdl_nm).kR(sub_idx);
            simpleMdl_prm.kPp(iS) = betas_fullList.Ep.(mdl_nm).kP(sub_idx);
            simpleMdl_prm.kRm(iS) = betas_fullList.Em.(mdl_nm).kR(sub_idx);
            simpleMdl_prm.kPm(iS) = betas_fullList.Em.(mdl_nm).kP(sub_idx);
        otherwise
            error('case not ready yet');
    end
    simpleMdl_prm.kEp(iS) = betas_fullList.Ep.(mdl_nm).kEffort(sub_idx);
    simpleMdl_prm.kEm(iS) = betas_fullList.Em.(mdl_nm).kEffort(sub_idx);
    simpleMdl_prm.kFp(iS) = betas_fullList.Ep.(mdl_nm).kFatigue(sub_idx);
    simpleMdl_prm.kFm(iS) = betas_fullList.Em.(mdl_nm).kFatigue(sub_idx);
end % subject list

%% perform correlation tests
goodSubs = ~isnan(bayesian_prm.kEp);
R_Rp = corrcoef(simpleMdl_prm.kMp(goodSubs), bayesian_prm.kR(goodSubs));
R_Rm = corrcoef(simpleMdl_prm.kMm(goodSubs), bayesian_prm.kR(goodSubs));
R_Ep = corrcoef(simpleMdl_prm.kEp(goodSubs), bayesian_prm.kEp(goodSubs));
R_Em = corrcoef(simpleMdl_prm.kEm(goodSubs), bayesian_prm.kEm(goodSubs));
R_Fp = corrcoef(simpleMdl_prm.kFp(goodSubs), bayesian_prm.kFp(goodSubs));
R_Fm = corrcoef(simpleMdl_prm.kFm(goodSubs), bayesian_prm.kFm(goodSubs));