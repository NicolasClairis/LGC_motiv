%% script performing the correlation between fMRI betas in a given ROI and 
% behavioral parameters extracted from the model

%% define subjects to include
study_nm = 'study1';
condition = subject_condition();
[subject_id, NS] = LGCM_subject_selection(study_nm, condition, 'all');

%% define GLM to work on
GLM_str = inputdlg('Which GLM?');
GLM = str2double(GLM_str);

%% define fMRI ROI to use
[con_vec_all,...
    ~, ~, ~,...
    con_names,...
    ROI_coords, ttest_ROI] = ROI_extraction_group('study1', GLM,...
    subject_id, condition, 0);
n_cons = size(con_vec_all, 1);
n_ROIs = size(con_vec_all,3);
if n_ROIs > 1
    error(['more than 1 ROI selected, mediation script cannot work that way',...
        'please focus on one and do it again.']);
end

%% define regression estimate to look for in the fMRI GLM
con_idx = listdlg('promptstring','Which contrast?',...
    'ListString',con_names);
% con_nm = con_names{con_idx};
con_nm = inputdlg('Contrast short name?');
con_data = NaN(1,NS);
con_data(:) = con_vec_all(con_idx, :, 1);

%% extract behavioural parameters
[prm] = prm_extraction(subject_id);
parameters = fieldnames(prm);
behavPrm_CID = prm.CID;
parameter_names = parameters;
parameter_names(strcmp(parameter_names,'CID'))=[]; % remove indication of subject ID


%% perform the linear regression
pval.signif = struct;
for iPrm = 1:length(parameters)
    prm_nm = parameters{iPrm};
    if ~strcmp(prm_nm,'CID')
        goodSubs = ~isnan(prm.(prm_nm));
        [b_tmp,~,stats_tmp] = glmfit(con_data(goodSubs), prm.(prm_nm)(goodSubs), 'normal');
        
        % name to store
        fMRI_f_prm_nm = [prm_nm,'_f_',ROI_coords.ROI_nm.ROI_1_shortName,'_',con_nm{1}];
        betas.(fMRI_f_prm_nm) = b_tmp;
        pval.(fMRI_f_prm_nm) = stats_tmp.p;
        
        % fit
        fitted.(fMRI_f_prm_nm) = glmval(b_tmp,con_data(goodSubs),'identity');
        
        %% store it if p is significant
        if stats_tmp.p(2) < 0.05
            pval.signif.(fMRI_f_prm_nm) = stats_tmp.p(2);
        end
    end % CID
end % behavioral parameter loop