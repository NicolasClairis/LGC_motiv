function[fatigue_measures] = fatigue_pool()
% [fatigue_measures] = fatigue_pool()
% fatigue_pool will extract and pool all the variables of the experiment
% that somehow relate to physical, mental or central fatigue.
%
% OUTPUTS
% fatigue_measures: structure with all the relevant information. It will
% include:
% - the list of subjects and condition used for the selection
% - the fatigue MPSTEFS and energy JPI-R questionnaire
% - the fatigue ratings the day of the experiment
% - task fatigue parameter kFp
% - the decrease in performance pre/post each run due to fatigue
% - the decrease in performance per effort level within the task
% - the number of covid infections

%% subject selection
study_nm = 'study1';
[study_nm, condition, subject_id, NS, genderFilter] = subject_selection(study_nm);

%% prepare the data of interest (to force all variables to be with the same dimensions)
[fatigue_measures.MPSTEFS_physical,...
    fatigue_measures.MPSTEFS_mental,...
    fatigue_measures.JPIR,...
    fatigue_measures.n_covid,...
    fatigue_measures.F1,...
    fatigue_measures.F2,...
    fatigue_measures.F3,...
    fatigue_measures.F4,...
    fatigue_measures.F4_F1,...
    fatigue_measures.F4_F3,...
    fatigue_measures.MVC_R1a,...
    fatigue_measures.MVC_R1b,...
    fatigue_measures.MVC_R2a,...
    fatigue_measures.MVC_R2b,...
    fatigue_measures.MVC_dR1,...
    fatigue_measures.MVC_dR2,...
    fatigue_measures.MVC_dTask,...
    fatigue_measures.perf_decrease_slope_E0,...
    fatigue_measures.perf_decrease_slope_E1,...
    fatigue_measures.perf_decrease_slope_E2,...
    fatigue_measures.perf_decrease_slope_E3,...
    fatigue_measures.choices_kFp] = deal(NaN(1,NS));
fatigue_varnames = fieldnames(fatigue_measures);
n_fatigue_vars = length(fatigue_varnames);
fatigue_mtrx = NaN(n_fatigue_vars, NS);

%% store subject-relevant information
fatigue_measures.sub_selection.condition = condition;
fatigue_measures.sub_selection.subject_id = subject_id;
fatigue_measures.sub_selection.NS = NS;
fatigue_measures.sub_selection.genderFilter = genderFilter;

%% load the data
% questionnaires related to fatigue (MPSTEFS + JPI-R) + nb of covid
% infections
[questionnaires] = extract_questionnaires(study_nm, subject_id, NS);
%% number of covid infections
fatigue_measures.n_covid = questionnaires.general.n_covid;

%% fatigue-related questionnaires
fatigue_measures.MPSTEFS_physical = questionnaires.Motivation.MPSTEFS_physical;
fatigue_measures.MPSTEFS_mental = questionnaires.Motivation.MPSTEFS_mental;
fatigue_measures.JPIR = questionnaires.Motivation.JPIR;

%% fatigue ratings on the day of the experiment
[~,...
    preMRS_fatigue, postMRS_fatigue,...
    prefMRI_fatigue, postfMRI_fatigue,...
    deltaFatiguePrePostExp, deltaFatiguePrePostfMRI] = extract_subjective_fatigue_ratings(study_nm, subject_id, NS);
fatigue_measures.F1 = preMRS_fatigue;
fatigue_measures.F2 = postMRS_fatigue;
fatigue_measures.F3 = prefMRI_fatigue;
fatigue_measures.F4 = postfMRI_fatigue;
fatigue_measures.F4_F1 = deltaFatiguePrePostExp;
fatigue_measures.F4_F3 = deltaFatiguePrePostfMRI;

%% performance: change in maximum voluntary contraction (MVC) due to fatigue (ie decrease in physical capacity per block and across blocks)
maxPerf = maxPerfEvolutionAcrossRuns_group([], 0, 0, study_nm, condition, subject_id, NS);
% physical task
% max perf
fatigue_measures.MVC_R1a = maxPerf.Ep.allData(1,:);
fatigue_measures.MVC_R1b = maxPerf.Ep.allData(2,:);
fatigue_measures.MVC_R2a = maxPerf.Ep.allData(3,:);
fatigue_measures.MVC_R2b = maxPerf.Ep.allData(4,:);
% difference after-before each session
fatigue_measures.MVC_dR1 = maxPerf.Ep.allData(2,:) - maxPerf.Ep.allData(1,:);
fatigue_measures.MVC_dR2 = maxPerf.Ep.allData(4,:) - maxPerf.Ep.allData(3,:);
% difference between final and initial maxPerf
fatigue_measures.MVC_dTask = maxPerf.Ep.allData(4,:) - maxPerf.Ep.allData(1,:);

%% performance: change in performance due to fatigue within each block
b_peakF_perEch_f_time = peakF_f_time(study_nm, subject_id, condition, 0); % slope of decrease in peak force per effort level
fatigue_measures.perf_decrease_slope_E0 = b_peakF_perEch_f_time(1,:);
fatigue_measures.perf_decrease_slope_E1 = b_peakF_perEch_f_time(2,:);
fatigue_measures.perf_decrease_slope_E2 = b_peakF_perEch_f_time(3,:);
fatigue_measures.perf_decrease_slope_E3 = b_peakF_perEch_f_time(4,:);

%% choice-derived behavioral fatigue parameter kFp
prm = prm_extraction(study_nm, subject_id);
fatigue_measures.choices_kFp = prm.kFp;

%% load all into fatigue_mtrx for future tests
for iVar = 1:n_fatigue_vars
    var_nm = fatigue_varnames{iVar};
    fatigue_mtrx(iVar, :) = fatigue_measures.(var_nm);
end

end % function