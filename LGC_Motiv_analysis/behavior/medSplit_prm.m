function[low_prm_subs, high_prm_subs, prm_allSubs] = medSplit_prm(study_nm, subject_id, prm_nm)
% [low_prm_subs, high_prm_subs, prm_allSubs] = medSplit_prm(study_nm, subject_id, prm_nm)
% medSplit_prm will allow you to select the parameter and the area
% that you want in order to perform a median split based on it.
%
% INPUTS
% study_nm: study name ('study1'/'study2')
%
% subject_id: cell with list of subjects for which you will extract the
% parameters
%
% prm_nm: parameter name (allows you to select which parameter to extract)
%
% OUTPUTS
% low_prm_subs: index of the subjects with low levels of the selected
% parameter
%
% high_prm_subs: index of the subjects with high levels of the selected
% parameter
%
% prm_allSubs: level of parameter selected for each subject
%
% Developed by N. Clairis - march 2023


%% define parameter and ROI you want to focus on
prm = prm_extraction(study_nm, subject_id, 'bayesian', '3');
prm_allSubs = prm.(prm_nm);

%% perform a median split based on the parameter selected in the ROI
% selected
med_prm_allSubs = median(prm_allSubs,'omitnan');
% extract index of participants with low or high level of behavioral
% parameter
low_prm_subs = prm_allSubs <= med_prm_allSubs;
high_prm_subs = prm_allSubs > med_prm_allSubs;

end % function