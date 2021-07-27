function[subject_id, NS] = LGCM_subject_selection(study_nm)
% [subject_id, NS] = LGCM_subject_selection(study_nm)
% LGCM_subject_selection will select the subject names and number of
% subjects for the current study.
%
% INPUTS
% study_nm: definition of the study on which you want to analyze the data
%
% OUTPUTS
% subject_id: list of subject names
%
% NS: number of subjects

switch study_nm
    case 'fMRI_pilots'
%         subject_id = {'pilot_s1','pilot_s2'};
        subject_id = {'pilot_s1'};
    case 'study1'
        subject_id = {};
    case 'study2_clinical'
        subject_id = {};
end

NS = length(subject_id);

end % function