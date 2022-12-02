% convert delta net value data sent by Arthur to homogenize with previous data

%% clear variables to avoid interference
clear;
%% general parameters
study_nm = 'study1';
nTrialsPerRun = 54;
nRuns = 4;

%% working directories
whichPc = 'home';
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
if exist([gitPath, 'bayesian_pChoice_data.mat'],'file')
    bayesian_pChoice = getfield(load([gitPath, 'bayesian_pChoice_data.mat'],...
        'bayesian_pChoice'),'bayesian_pChoice');
end
%% extraction model 3
NVmodel3_data = getfield(load('CID_NV_Proba.mat'),'CID_NV_Proba');
NS = size(NVmodel3_data, 2);
bayesian_pChoice.mdl3.subject_id = deal(cell(1,NS));
for iS = 1:NS
    sub_n = NVmodel3_data(1,iS);
    if sub_n < 10
        sub_nm = ['00',num2str(sub_n)];
    elseif sub_n >= 10 && sub_n < 100
        sub_nm = ['0',num2str(sub_n)];
    elseif sub_n >= 100
        sub_nm = num2str(sub_n);
    end
    bayesian_pChoice.mdl3.subject_id{iS} = sub_nm;
    for iRun = 1:nRuns
        run_nm = ['run',num2str(iRun)];
        trial_idx = (1:nTrialsPerRun) + nTrialsPerRun*(iRun - 1) +...
            1 + 216; % first line is subject_id and 216 next are delta net value
        bayesian_pChoice.mdl3.(['CID',sub_nm]).(run_nm) = NVmodel3_data(trial_idx,iS);
    end % run loop
end % subject loop

%% save data updated
save([gitPath, 'bayesian_pChoice_data.mat'],'bayesian_pChoice');