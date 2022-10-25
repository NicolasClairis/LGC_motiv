function[low_met_subs, high_met_subs, metabolite_nm, ROI_nm] = medSplit_metabolites(subject_id)
% [low_met_subs, high_met_subs, metabolite_nm, ROI_nm] = medSplit_metabolites(subject_id)
% medSplit_metabolites will allow you to select the metabolite and the area
% that you want in order to perform a median split based on it.
%
% INPUTS
% subject_id: cell with list of subjects for which you will extract the
% metabolites
%
% OUTPUTS
% low_met_subs: index of the subjects with low levels of the selected
% metabolite
%
% high_met_subs: index of the subjects with high levels of the selected
% metabolite
%
% metabolite_nm: name of the selected metabolite
%
% ROI_nm: name of the brain area for which the metabolite has been selected
%
% Developed by N. Clairis - june 2022


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

%% extract all metabolites
[metabolites] = metabolite_load(subject_id);
% focus on metabolite and brain area selected
metabolite_allSubs = metabolites.(ROI_nm).(metabolite_nm);

%% perform a median split based on the metabolite selected in the ROI
% selected
med_metabolite_allSubs = median(metabolite_allSubs,'omitnan');
% extract index of participants with low or high level of metabolites
low_met_subs = metabolite_allSubs <= med_metabolite_allSubs;
high_met_subs = metabolite_allSubs > med_metabolite_allSubs;

end % function