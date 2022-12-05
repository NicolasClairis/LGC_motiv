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
bayesian_mdlN = '3';
[bayesian_prm] = prm_extraction(study_nm, subject_id, 'bayesian', bayesian_mdlN);
bayesian_parameters = fieldnames(bayesian_prm);

%% perform simple behavioral model
simple_mdlN = '3';
[simpleMdl_prm] = prm_extraction(study_nm, subject_id, 'simple', simple_mdlN);
simple_parameters = fieldnames(simpleMdl_prm);

%% perform correlation tests
goodSubs = ~isnan(bayesian_prm.kEp);
switch simple_mdlN
    case {'1','3'} % money/R-P
        R_sMp_bR = corrcoef(simpleMdl_prm.kMp(goodSubs), bayesian_prm.kR(goodSubs));
        R_sMm_bR = corrcoef(simpleMdl_prm.kMm(goodSubs), bayesian_prm.kR(goodSubs));
        R_sMp_bP = corrcoef(simpleMdl_prm.kMp(goodSubs), bayesian_prm.kR(goodSubs));
        R_sMm_bP = corrcoef(simpleMdl_prm.kMm(goodSubs), bayesian_prm.kR(goodSubs));
        % store into output
        corr.sMp_bR = R_sMp_bR(2,1);
        corr.sMm_bR = R_sMm_bR(2,1);
        corr.sMp_bP = R_sMp_bP(2,1);
        corr.sMm_bP = R_sMm_bP(2,1);
    case {'2','4'} % R/P
        R_sRp_bR = corrcoef(simpleMdl_prm.kRp(goodSubs), bayesian_prm.kR(goodSubs));
        R_sPp_bP = corrcoef(simpleMdl_prm.kPp(goodSubs), bayesian_prm.kP(goodSubs));
        R_sRm_bR = corrcoef(simpleMdl_prm.kRm(goodSubs), bayesian_prm.kR(goodSubs));
        R_sPm_bR = corrcoef(simpleMdl_prm.kPm(goodSubs), bayesian_prm.kP(goodSubs));
        % store into output
        corr.sRp_bR = R_sRp_bR(2,1);
        corr.sPp_bP = R_sPp_bP(2,1);
        corr.sRm_bR = R_sRm_bR(2,1);
        corr.sPm_bR = R_sPm_bR(2,1);
end
% effort
R_sEp_bEp = corrcoef(simpleMdl_prm.kEp(goodSubs), bayesian_prm.kEp(goodSubs));
R_sEm_bEm = corrcoef(simpleMdl_prm.kEm(goodSubs), bayesian_prm.kEm(goodSubs));
% store into output
corr.sEp_bEp = R_sEp_bEp(2,1);
corr.sEm_bEm = R_sEm_bEm(2,1);

% fatigue
switch simple_mdlN
    case {'3','4'}
        R_sFp_bFp = corrcoef(simpleMdl_prm.kFp(goodSubs), bayesian_prm.kFp(goodSubs));
        R_sFm_bFm = corrcoef(simpleMdl_prm.kFm(goodSubs), bayesian_prm.kLm(goodSubs));
        % store into output
        corr.sFp_bFp = R_sFp_bFp(2,1);
        corr.sFm_bFm = R_sFm_bFm(2,1);
end
