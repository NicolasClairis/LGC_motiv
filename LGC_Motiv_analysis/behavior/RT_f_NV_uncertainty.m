function[RT, b_RT, pval] = RT_f_NV_uncertainty(study_nm, subject_id, condition, figDisp)
% [RT, b_RT, pval] = RT_f_NV_uncertainty(study_nm, subject_id, condition, figDisp)
% RT_f_NV_uncertainty will look at how reaction times (RT) during choice
% vary with net value depending on the uncertainty rating and uncertainty
% derived from the model, looking at all tasks together and each task
% separately.
%
% INPUTS
% study_nm: study name
%
% subject_id: list of subjects (can be left empty)
%
% condition: condition used
%
% figDisp: display figure (1) or not (0)
%
% OUTPUTS
% RT: structure with RT
%
% b_RT: beta for linear regressions
%
% pval: p.value for different linear regressions

%% subject selection
if ~exist('study_nm','var') || isempty(study_nm) ||...
        ~ismember(study_nm,{'study1','study2'})
    study_nm = 'study1';
end
if ~exist('condition','var') || isempty(condition)
    condition = subject_condition;
end
if ~exist('subject_id','var') || isempty(subject_id)
    [subject_id, NS] = LGCM_subject_selection(study_nm, condition);
else
    NS = length(subject_id);
end
if ~exist('figDisp','var') || isempty(figDisp)
    figDisp = 1;
end

%% working directories
computerRoot = LGCM_root_paths;
dataRoot = [computerRoot, filesep, study_nm, filesep];

for iS = 1:NS
    sub_nm = subject_id{iS};
    subBehaviorFolder = [dataRoot, filesep, 'CID',sub_nm, filesep, 'behavior',filesep];
    runs = runs_definition(study_nm, sub_nm, condition);
    for iRun = 1:runs.nb_runs.Ep
        jRun = runs.Ep.runsToKeep(iRun);
        run_nm = num2str(jRun);
        run_nm_bis = ['run',num2str(iRun)];
    end % run loop
end % subject loop
end % function