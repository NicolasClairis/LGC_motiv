function[male_CIDS, female_CIDS, male_NS, female_NS, condition] = subject_selection_per_sex()
% [male_CIDS, female_CIDS, male_NS, female_NS, condition] = subject_selection_per_sex()
 % subject_selection_per_sex will extract the male and female subject ids
 % independently according to the condition selected.
 %
 % OUTPUTS
 % male_CIDS, female_CIDS: cell with list of male and female, respectively,
 % subject names, depending on the condition selected
 %
 % male_NS, female_NS: number of subjects for male and female,
 % respectively, conditions
 %
 % condition: string indicating the condition selected in input


%% study name
study_nm = 'study1';
%% extract general condition
condition = subject_condition;
%% load males
[male_CIDS, male_NS] = LGCM_subject_selection(study_nm, condition, 'males');
%% load females
[female_CIDS, female_NS] = LGCM_subject_selection(study_nm, condition, 'females');

end % function