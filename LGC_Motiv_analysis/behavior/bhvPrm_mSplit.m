function[low_prm_CID, high_prm_CID, behavioralPrm_mSplit] = bhvPrm_mSplit(study_nm, subject_id)
% [low_prm_CID, high_prm_CID, behavioralPrm_mSplit] = bhvPrm_mSplit(study_nm, subject_id)
% bhvPrm_mSplit will extract the subjects with low vs high levels of the
% selected parameter.
%
% INPUTS
% study_nm: study name ('study1'/'study2')
%
% subject_id: list of subjects to look at
%
% OUTPUTS
% low_prm_CID: name of the subjects with low parameter
%
% high_prm_CID: name of the subjects with high parameter
%
% behavioralPrm_mSplit: structure containing mean, SEM and SD for the parameter selected
% in the low vs high group (and also for everyone)

%% extract behavioural parameters
[prm] = prm_extraction(study_nm, subject_id);
parameters = fieldnames(prm);
behavPrm_CID = prm.CID;
parameter_names = parameters;
% remove indication of subject ID
parameter_names(strcmp(parameter_names,'CID')) = [];

%% ask which parameter to use for median split
which_prm = listdlg('PromptString','Which parameter?',...
    'ListString',parameter_names,...
    'SelectionMode','single');
bhvPrm = prm.(parameter_names{which_prm});

%% median split
med_prm = median(bhvPrm, 'omitnan');
low_prm = bhvPrm <= med_prm;
high_prm = bhvPrm > med_prm;
low_prm_CID = behavPrm_CID(low_prm);
high_prm_CID = behavPrm_CID(high_prm);

% extract corresponding values of the parameter
[behavioralPrm_mSplit.low.mean,...
    behavioralPrm_mSplit.low.sem,...
    behavioralPrm_mSplit.low.sd] = mean_sem_sd(bhvPrm(low_prm),2);
[behavioralPrm_mSplit.high.mean,...
    behavioralPrm_mSplit.high.sem,...
    behavioralPrm_mSplit.high.sd] = mean_sem_sd(bhvPrm(high_prm),2);
[behavioralPrm_mSplit.all.mean,...
    behavioralPrm_mSplit.all.sem,...
    behavioralPrm_mSplit.all.sd] = mean_sem_sd(bhvPrm,2);

end % function