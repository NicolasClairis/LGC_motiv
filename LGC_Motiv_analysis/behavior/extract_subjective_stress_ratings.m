function[stress,...
    preMRS_stress, postMRS_stress,...
    prefMRI_stress, postfMRI_stress,...
    deltaStressPrePostExp] = extract_subjective_stress_ratings(study_nm, subject_id, NS)
% [stress, preMRS_stress, postMRS_stress,...
%     prefMRI_stress, postfMRI_stress,...
%     deltaStressPrePostExp] = extract_subjective_stress_ratings(study_nm, subject_id, NS)
% extract_subjective_stress_ratings will extract the subjective stress
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
% stress: structure containing all the following variables:
%
% preMRS_stress: subjective stress when subject arrives (before MRS)
%
% postMRS_stress: subjective stress when subject came out of MRS
%
% prefMRI_stress: subjective stress after training for the behavioral
% task and before behavioral task in fMRI
%
% postfMRI_stress: subjective stress after behavioral task in fMRI
%
% deltaStressPrePostExp: difference between initial subjective stress
% rating and final level of stress
%
% Designed by N.Clairis - november 2023


%% extract stress levels
[preMRS_stress, postMRS_stress,...
    prefMRI_stress, postfMRI_stress,...
    deltaStressPrePostExp] = deal(NaN(1,NS));

% load stress
[excelReadGeneralFile] = load_gal_data_bis(study_nm);
%% loop through subjects
for iS = 1:NS
    sub_nm = ['CID',subject_id{iS}];
    stress_sub_idx = strcmp(sub_nm, excelReadGeneralFile.CID);
    
    preMRS_stress(iS)          = excelReadGeneralFile.StressPr__MRS(stress_sub_idx);
    postMRS_stress(iS)          = excelReadGeneralFile.StressPost_MRS(stress_sub_idx);
    prefMRI_stress(iS)          = excelReadGeneralFile.StressPr__fMRI(stress_sub_idx);
    postfMRI_stress(iS)          = excelReadGeneralFile.StressPost_fMRI(stress_sub_idx);
    deltaStressPrePostExp(iS)  = postfMRI_stress(iS) - preMRS_stress(iS);
end

%% create structure containing all the previous data
stress.preMRS = preMRS_stress;
stress.postMRS = postMRS_stress;
stress.prefMRI = prefMRI_stress;
stress.postfMRI = postfMRI_stress;
stress.deltaPrePostExp = deltaStressPrePostExp;

end % function