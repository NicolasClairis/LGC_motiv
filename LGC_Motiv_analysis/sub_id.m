function[study_nm, condition, sex, subject_id, NS] = sub_id()
% [study_nm, condition, sex, subject_id, NS] = sub_id
% sub_id will provide you the list of subjects to look for.
%
% OUTPUTS
% study_nm: string with study name ('study1' by default now)
%
% condition: string with condition selected for this analysis
%
% sex: string with sex selected for this analysis
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

%% selection of the sex
% sex_names = {'all','males','females'};
% which_sex = listdlg('PromptString','Which sex?',...
%     'SelectionMode','single','ListString',sex_names);
% sex = sex_names{which_sex};
sex = 'all';

%% extract list of subjects
[subject_id, NS] = LGCM_subject_selection(study_nm, condition, sex);

end % function