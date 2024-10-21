function[study_nm, condition, subject_id, NS, genderFilter] = subject_selection(study_nm, condition, subject_id, NS, genderFilter)
% [study_nm, condition, subject_id, NS, genderFilter] = subject_selection(study_nm, condition, subject_id, NS, genderFilter)
%subject_selection allows to select the subjects and conditions of
%interest.
%
% INPUTS (any input can be left empty and will be selected accordingly)
% study_nm: study name ('study1' by default if left empty)
%
% condition: condition for subjects and runs to include (list will be
% proposed and script will ask you which you want if left empty)
%
% subject_id: list of subjects to include (will be determined automatically
% if left empty)
%
% NS: number of subjects (will be determined automatically if left empty)
%
% genderFilter: include only subjects of one sex ('males'/'females') or
% everybody ('all'). The latter will be done by default.
%
% OUTPUTS
% same variables as inputs, ie:
% study_nm: study name
%
% condition: condition for subjects and runs to include
%
% subject_id: list of subjects to include
%
% NS: number of subjects
%
% genderFilter: include only subjects of one sex ('males'/'females') or
% everybody ('all').

%% study selection
if ~exist('study_nm','var') || isempty(study_nm) || ~ismember(study_nm,{'study1','study2'})
   study_nm = 'study1'; 
end

%% condition
if ~exist('condition','var') || isempty(condition)
    condition = subject_condition;
end

%% gender filter
if ~exist('genderFilter','var') || isempty(genderFilter) ||...
        ~ismember(genderFilter,{'females','males','all'})
   genderFilter = 'all'; 
end

%% list of subjects
if ~exist('subject_id','var') || isempty(subject_id)
    [subject_id, NS] = LGCM_subject_selection(study_nm, condition, genderFilter);
else
    if ~exist('NS','var') || isempty(NS) || NS <= 0
        NS = length(subject_id);
    end
end

end % function