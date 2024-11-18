function[stress_measures] = stress_pool(study_nm, condition, subject_id, NS, genderFilter)
% [stress_measures] = stress_pool(study_nm, condition, subject_id, NS, genderFilter)
% stress_pool will extract and pool all the variables of the experiment
% that somehow relate to stress.
%
% INPUTS
% study_nm: study name
%
% condition: condition to use
%
% subject_id: subject list
%
% NS: number of subjects
%
% genderFilter: use of gender filter?
%
% OUTPUTS
% stress_measures: structure with all the relevant information. It will
% include:
% - the list of subjects and condition used for the selection
% - the stress questionnaires (STAI-T, PSS14 and SIAS)
% - the stress ratings the day of the experiment
% - cortisol levels before/after MRS and before/after fMRI

%% subject selection
if ~exist('study_nm','var') && ~exist('condition','var') &&...
        ~exist('subject_id','var') && ~exist('NS','var') &&...
        ~exist('genderFilter','var')
    [study_nm, condition, subject_id, NS, genderFilter] = subject_selection;
else
    % study
    if ~exist('study_nm','var')
        study_nm = 'study1';
    end
    % condition
    if ~exist('condition','var')
        condition = subject_condition;
    end
    % subject list
    if ~exist('subject_id','var') || ~exist('NS','var')
        [subject_id, NS] = LGCM_subject_selection(study_nm, condition);
    end
end

%% prepare the data of interest (to force all variables to be with the same dimensions)
[stress_measures.STAI_T,...
    stress_measures.SIAS,...
    stress_measures.PSS14,...
    stress_measures.S1,...
    stress_measures.S2,...
    stress_measures.S3,...
    stress_measures.S4,...
    stress_measures.S4_S1,...
    stress_measures.S4_S3,...
    stress_measures.CORT1,...
    stress_measures.CORT2,...
    stress_measures.CORT3,...
    stress_measures.CORT4] = deal(NaN(1,NS));
stress_varnames = fieldnames(stress_measures);
n_stress_vars = length(stress_varnames);
stress_mtrx = NaN(n_stress_vars, NS);

%% store subject-relevant information
stress_measures.sub_selection.condition = condition;
stress_measures.sub_selection.subject_id = subject_id;
stress_measures.sub_selection.NS = NS;
stress_measures.sub_selection.genderFilter = genderFilter;

%% load the data
% questionnaires related to fatigue (MPSTEFS + JPI-R) + nb of covid
% infections
[questionnaires] = extract_questionnaires(study_nm, subject_id, NS);

%% stress-related questionnaires
stress_measures.STAI_T(:) = questionnaires.stress_anxiety.STAI_T;
stress_measures.SIAS(:) = questionnaires.stress_anxiety.SIAS;
stress_measures.PSS14(:) = questionnaires.stress_anxiety.PSS14;

%% stress ratings on the day of the experiment
[~,...
    preMRS_stress, postMRS_stress,...
    prefMRI_stress, postfMRI_stress,...
    deltaStressPrePostExp, deltaStressPrePostfMRI] = extract_subjective_stress_ratings(study_nm, subject_id, NS);
stress_measures.S1(:) = preMRS_stress;
stress_measures.S2(:) = postMRS_stress;
stress_measures.S3(:) = prefMRI_stress;
stress_measures.S4(:) = postfMRI_stress;
stress_measures.S4_S1(:) = deltaStressPrePostExp;
stress_measures.S4_S3(:) = deltaStressPrePostfMRI;

%% cortisol levels on the day of the experiment
[CORT_data] = load_CORT(study_nm, subject_id);
stress_measures.CORT1(:) = CORT_data.CORT(1,:);
stress_measures.CORT2(:) = CORT_data.CORT(2,:);
stress_measures.CORT3(:) = CORT_data.CORT(3,:);
stress_measures.CORT4(:) = CORT_data.CORT(4,:);

%% load all into fatigue_mtrx for future tests
for iVar = 1:n_stress_vars
    var_nm = stress_varnames{iVar};
    stress_mtrx(iVar, :) = stress_measures.(var_nm);
end

end % function