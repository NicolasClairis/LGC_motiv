function[fatigue,...
    preMRS_fatigue, postMRS_fatigue,...
    prefMRI_fatigue, postfMRI_fatigue,...
    deltaFatiguePrePostExp, deltaFatiguePrePostfMRI] = extract_subjective_fatigue_ratings(study_nm, subject_id, NS)
% [fatigue, preMRS_fatigue, postMRS_fatigue,...
%     prefMRI_fatigue, postfMRI_fatigue,...
%     deltaFatiguePrePostExp, deltaFatiguePrePostfMRI] = extract_subjective_fatigue_ratings(study_nm, subject_id, NS)
% extract_subjective_fatigue_ratings will extract the subjective fatigue
% ratings (on a visual analogue scale from 1 to 10) for the subjecte
% entered in input.
%
% INPUTS
% study_nm: study name (eg 'study1')
% 
% subject_id: list of subjects
%
% NS: number of subjects
%
% OUTPUTS
% fatigue: structure containing all the following variables:
%
% preMRS_fatigue: subjective fatigue when subject arrives (before MRS)
%
% postMRS_fatigue: subjective fatigue when subject came out of MRS
%
% prefMRI_fatigue: subjective fatigue after training for the behavioral
% task and before behavioral task in fMRI
%
% postfMRI_fatigue: subjective fatigue after behavioral task in fMRI
%
% deltaFatiguePrePostExp: difference between initial subjective fatigue
% rating and final level of fatigue
%
% deltaFatiguePrePostfMRI: difference between subjective fatigue rating
% before vs after the task in fMRI
%
% Designed by N.Clairis - november 2023


%% extract fatigue levels
[preMRS_fatigue, postMRS_fatigue,...
    prefMRI_fatigue, postfMRI_fatigue,...
    deltaFatiguePrePostExp, deltaFatiguePrePostfMRI] = deal(NaN(1,NS));

% load fatigue
[excelReadGeneralFile] = load_gal_data_bis(study_nm);
%% loop through subjects
for iS = 1:NS
    sub_nm = ['CID',subject_id{iS}];
    fatigue_sub_idx = strcmp(sub_nm, excelReadGeneralFile.CID);
    
    preMRS_fatigue(iS)          = excelReadGeneralFile.FatiguePr__MRS(fatigue_sub_idx);
    postMRS_fatigue(iS)          = excelReadGeneralFile.FatiguePost_MRS(fatigue_sub_idx);
    prefMRI_fatigue(iS)          = excelReadGeneralFile.FatiguePr__fMRI(fatigue_sub_idx);
    postfMRI_fatigue(iS)          = excelReadGeneralFile.FatiguePost_fMRI(fatigue_sub_idx);
    deltaFatiguePrePostExp(iS)  = postfMRI_fatigue(iS) - preMRS_fatigue(iS);
    deltaFatiguePrePostfMRI(iS) = postfMRI_fatigue(iS) - prefMRI_fatigue(iS);
end

%% create structure containing all the previous data
fatigue.preMRS = preMRS_fatigue;
fatigue.postMRS = postMRS_fatigue;
fatigue.prefMRI = prefMRI_fatigue;
fatigue.postfMRI = postfMRI_fatigue;
fatigue.deltaPrePostExp = deltaFatiguePrePostExp;
fatigue.deltaFatiguePrePostfMRI = deltaFatiguePrePostfMRI;

end % function