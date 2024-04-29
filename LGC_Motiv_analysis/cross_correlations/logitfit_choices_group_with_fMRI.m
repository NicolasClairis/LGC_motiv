function[betas, pvalues] = logitfit_choices_group_with_fMRI(study_nm, subject_id, computerRoot)
% [betas, pvalues] = logitfit_choices_group_with_fMRI(study_nm, subject_id, computerRoot, figDispGroup)
%
% INPUTS
% study_nm: study name
%
% subject_id: list of subjects to extract (will be asked if left empty)
%
% computerRoot: root where data is stored (will be asked if left empty)
%
% OUTPUT
% betas: structure with betas
%
% pvalues: structure with p.values

%% if root not defined => ask for it
if ~exist('computerRoot','var') || isempty(computerRoot)
    computerRoot = LGCM_root_paths;
end

%% study name
if ~exist('study_nm','var') || isempty(study_nm)
    study_nm = 'study1';
end

%% working directories
studyBehaviorFolder = [computerRoot, filesep, study_nm, filesep];
resultFolder_a = [studyBehaviorFolder,'results',filesep];
if ~exist(resultFolder_a,'dir')
    mkdir(resultFolder_a);
end
resultFolder = [resultFolder_a,'figures',filesep];
if ~exist(resultFolder,'dir')
    mkdir(resultFolder);
end

%% subject selection
if ~exist('subject_id','var') || isempty(subject_id)
    [condition] = subject_condition;
    [subject_id, NS] = LGCM_subject_selection(study_nm, condition);
else
    NS = length(subject_id);
end
% store subject list to know which beta corresponds to which subject
betas.subList = subject_id;
%% initialize variables of interest
for iPM = 1:2
    switch iPM
        case 1
            task_id = 'Ep';
        case 2
            task_id = 'Em';
    end
    [betas.(task_id).bias,...
        betas.(task_id).bR,...
        betas.(task_id).b_VS_R,...
        betas.(task_id).bP,...
        betas.(task_id).b_aINS_P,...
        betas.(task_id).bE,...
        betas.(task_id).b_dmPFC_E,...
        betas.(task_id).bF,...
        betas.(task_id).b_dmPFC_F] = deal(NaN(1,NS));
    prm_names.(task_id) = fieldnames(betas.(task_id));
    nPrm.(task_id) = length(prm_names.(task_id));
end % physical/mental loop

%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    % load individual data
    [betas_tmp] = logitfit_choices_with_fMRI(computerRoot, study_nm, sub_nm, condition);
    
    % pool data across subjects
    for iPM = 1:2
        switch iPM
            case 1
                task_id = 'Ep';
            case 2
                task_id = 'Em';
        end
        % extract betas
        for iPrm = 1:nPrm.(task_id)
            prm_nm = prm_names.(task_id){iPrm};
            betas.(task_id).(prm_nm)(iS) = betas_tmp.(task_id).(prm_nm);
        end
    end % physical/mental loop
    disp(['Subject ',sub_nm,' - ',num2str(iS),'/',num2str(NS),' done']);
end % subject loop

%% average data
for iPM = 1:2
    switch iPM
        case 1
            task_id = 'Ep';
        case 2
            task_id = 'Em';
    end
    
    % average betas
    for iPrm = 1:nPrm.(task_id)
        prm_nm = prm_names.(task_id){iPrm};
        [betas.mean.(task_id).(prm_nm),...
            betas.sem.(task_id).(prm_nm),...
            betas.sd.(task_id).(prm_nm)] = mean_sem_sd(betas.(task_id).(prm_nm),2);
        
        % test how significant betas are
        [~,pvalues.(task_id).(prm_nm)] = ttest(betas.(task_id).(prm_nm));
    end % parameter loop
end % physical/mental loop

end % function