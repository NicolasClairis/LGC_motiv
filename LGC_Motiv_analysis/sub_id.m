function[study_nm, condition, gender, subject_id, NS] = sub_id()
% [study_nm, condition, gender, subject_id, NS] = sub_id
% sub_id will provide you the list of subjects to look for.
%
% OUTPUTS
% study_nm: string with study name ('study1' by default now)
%
% condition: string with condition selected for this analysis
%
% gender: string with gender selected for this analysis
%
% subject_id: list of subjects to be included
%
% NS: number of subjects included

%% selection of the study
% study_names = {'study1','study2'};
% which_study = listdlg('PromptString','Which study?',...
%     'SelectionMode','single','ListString',study_names);
% study_nm = study_names{which_study};
study_nm = 'study1';
%% selection of the condition (which subjects to include
condition = subject_condition;

%% selection of the gender
% gender_names = {'all','males','females'};
% which_gender = listdlg('PromptString','Which gender?',...
%     'SelectionMode','single','ListString',gender_names);
% gender = gender_names{which_gender};
gender = 'all';

%% extract list of subjects
[subject_id, NS] = LGCM_subject_selection(study_nm, condition, gender);

end % function