% convert delta net value data sent by Arthur to homogenize with previous data
%
% CID_NV_PROBA:
% line 1 = subject_id
% lines 2-217 = deltaNV (SV high E - SV low E) across all trials and runs
% lines 218-433 = p(choice = high E) across all trials and runs

%% clear variables to avoid interference
clear;
%% general parameters
study_nm = 'study1';
nTrialsPerRun = 54;
nRuns = 4;

%% working directories
whichPc = 'Lab';
switch whichPc
    case 'Lab'
        gitRoot = fullfile('C:','Users','clairis','Desktop');
    case 'home'
        gitRoot = fullfile('C:','Users','Loco','Documents');
end
gitPath = [gitRoot, filesep, 'GitHub',filesep,...
    'LGC_motiv',filesep,'LGC_Motiv_results',filesep,...
    study_nm, filesep,...
    'bayesian_modeling', filesep];

%% load model to modify
if exist([gitPath, 'bayesian_deltaNV_data.mat'],'file')
    bayesian_deltaNV = getfield(load([gitPath, 'bayesian_deltaNV_data.mat'],...
        'bayesian_deltaNV'),'bayesian_deltaNV');
end
%% extraction model 3
NVmodel3_data = getfield(load('CID_NV_Proba.mat'),'CID_NV_PROBA');
NS = size(NVmodel3_data, 2);
bayesian_deltaNV.mdl_3.subject_id = convert_sub_id_from_num_to_cell(NVmodel3_data(1,:));
for iS = 1:NS
    sub_nm = bayesian_deltaNV.mdl_3.subject_id{iS};
    for iRun = 1:nRuns
        run_nm = ['run',num2str(iRun)];
        trial_idx = (1:nTrialsPerRun) + nTrialsPerRun*(iRun - 1) + 1; % first line is subject_id
        bayesian_deltaNV.mdl_3.(['CID',sub_nm]).(run_nm) = NVmodel3_data(trial_idx,iS);
    end % run loop
end % subject loop

%% save data updated
save([gitPath, 'bayesian_deltaNV_data.mat'],'bayesian_deltaNV');