% check if dmPFC slope can be explained by several metabolites => kEp

%% launch mediation
mediation_allMetabolites_fMRI_allBehavPrm;

%% try combining everyone
% prm_nm='kEm';
% behavPrm = prm.(prm_nm);
dmPFC_GSH_allSubs = metabolites.dmPFC.GSH;
dmPFC_Tau_allSubs = metabolites.dmPFC.Tau;
dmPFC_Lac_allSubs = metabolites.dmPFC.Lac;
okSubs = ~isnan(dmPFC_GSH_allSubs).*~isnan(dmPFC_Tau_allSubs).*...
    ~isnan(dmPFC_Lac_allSubs).*~isnan(con_data) == 1; %.*~isnan(behavPrm)

X_mb = [dmPFC_GSH_allSubs', dmPFC_Tau_allSubs', dmPFC_Lac_allSubs'];
[beta_mb_fMRI,~,stats_mb_fMRI] = glmfit(X_mb(okSubs,:), con_data(okSubs),'normal');


%% same model without Taurine
dmPFC_GSH_allSubs = metabolites.dmPFC.GSH;
dmPFC_Lac_allSubs = metabolites.dmPFC.Lac;
okSubs2 = ~isnan(dmPFC_GSH_allSubs).*~isnan(dmPFC_Lac_allSubs).*~isnan(con_data) == 1; %.*~isnan(behavPrm)

X_mb2 = [dmPFC_GSH_allSubs', dmPFC_Lac_allSubs'];
[beta_mb_fMRI2,~,stats_mb_fMRI2] = glmfit(X_mb2(okSubs2,:), con_data(okSubs2),'normal');

%% now trying Glu, GSH and Gly
dmPFC_GSH_allSubs = metabolites.dmPFC.GSH;
dmPFC_Glu_allSubs = metabolites.dmPFC.Glu;
dmPFC_Gly_allSubs = metabolites.dmPFC.Gly;
okSubs3 = ~isnan(dmPFC_GSH_allSubs).*~isnan(dmPFC_Glu_allSubs).*~isnan(dmPFC_Gly_allSubs).*~isnan(con_data) == 1; %.*~isnan(behavPrm)

X_mb3 = [dmPFC_GSH_allSubs', dmPFC_Glu_allSubs', dmPFC_Gly_allSubs'];
[beta_mb_fMRI3,~,stats_mb_fMRI3] = glmfit(X_mb3(okSubs3,:), con_data(okSubs3),'normal');