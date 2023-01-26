function[low_met_subs, high_met_subs, metabolite_nm, MRS_ROI_nm, metabolite_allSubs] = medSplit_metabolites(study_nm, subject_id)
% [low_met_subs, high_met_subs, metabolite_nm, MRS_ROI_nm, metabolite_allSubs] = medSplit_metabolites(study_nm, subject_id)
% medSplit_metabolites will allow you to select the metabolite and the area
% that you want in order to perform a median split based on it.
%
% INPUTS
% study_nm: study name ('study1'/'study2')
%
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
% MRS_ROI_nm: name of the brain area for which the metabolite has been selected
%
% metabolite_allSubs: level of metabolite selected for each subject
%
% Developed by N. Clairis - june 2022


%% define metabolite and ROI you want to focus on
[metabolite_allSubs, MRS_ROI_nm, metabolite_nm] = metabolite_extraction(study_nm, subject_id);

%% perform a median split based on the metabolite selected in the ROI
% selected
med_metabolite_allSubs = median(metabolite_allSubs,'omitnan');
% extract index of participants with low or high level of metabolites
low_met_subs = metabolite_allSubs <= med_metabolite_allSubs;
high_met_subs = metabolite_allSubs > med_metabolite_allSubs;

end % function