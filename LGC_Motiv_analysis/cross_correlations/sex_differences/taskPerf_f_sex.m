function[] = taskPerf_f_sex()
% perf_f_sex will compare different variables related to performance
% (maximal performance during calibration, average RT, difference in
% accuracy, overshoot, etc.) between males and females


%% subject selection
study_nm = 'study1';
[male_CIDS, female_CIDS, male_NS, female_NS, condition] = subject_selection_per_sex;

%% extract performance variables
% extract male data
[male_MVC, male_PCSA] = extract_grip_MVC_PCSA(study_nm, male_CIDS, male_NS);
[male_NMP] = extract_Nback_NMP(study_nm, male_CIDS, male_NS);
% extract female data
[female_MVC, female_PCSA] = extract_grip_MVC_PCSA(study_nm, female_CIDS, female_NS);
[female_NMP] = extract_Nback_NMP(study_nm, female_CIDS, female_NS);

%% perform comparisons
% compare maximal performance
[~,pval.MVC] = ttest2(male_MVC, female_MVC);
[~,pval.PCSA] = ttest2(male_PCSA, female_PCSA);
[~,pval.NMP] = ttest2(male_NMP, female_NMP);

end % function