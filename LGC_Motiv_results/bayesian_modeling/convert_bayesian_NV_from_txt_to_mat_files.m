%% extract net value based on bayesian models for each subject, each run 
% and each trial to be able to perform other analysis later on

%% working directories
whichPc = 'Lab';
switch whichPc
    case 'Lab'
        gitRoot = fullfile('C:','Users','clairis','Desktop');
    case 'Home'
        gitRoot = fullfile('C:','Users','Loco','Documents');
end
gitPath = [gitRoot, filesep, 'GitHub',filesep,...
    'LGC_motiv',filesep,'LGC_Motiv_results',...
    filesep,'bayesian_modeling', filesep];

%% general parameters
study_nm = 'study1';
nTrialsPerRun = 54;
nRuns = 4;

%% extraction model 1
model_NV_raw.mdl1 = load([gitPath,'net_values_biais.txt']);
sub_id_line1 = model_NV_raw.mdl1(1,:);
NV_data_mdl1 = model_NV_raw.mdl1(2:end,:);
[bayesian_deltaNV.mdl1.subject_id] = convert_sub_id_from_num_to_cell(sub_id_line1);
NS = length(bayesian_deltaNV.mdl1.subject_id);
for iS = 1:NS
    sub_nm = bayesian_deltaNV.mdl1.subject_id{iS};
    for iRun = 1:nRuns
        run_nm = ['run',num2str(iRun)];
        trial_idx = (1:nTrialsPerRun) + nTrialsPerRun*(iRun - 1);
        bayesian_deltaNV.mdl1.(['CID',sub_nm]).(run_nm) = NV_data_mdl1(trial_idx);
    end % run loop
end % subject loop

%% extraction model 2
model_NV_raw.mdl2 = load([gitPath,'net_values_classic.txt']);
sub_id_line2 = model_NV_raw.mdl2(1,:);
NV_data_mdl2 = model_NV_raw.mdl2(2:end,:);
[bayesian_deltaNV.mdl2.subject_id] = convert_sub_id_from_num_to_cell(sub_id_line2);
NS = length(bayesian_deltaNV.mdl2.subject_id);
for iS = 1:NS
    sub_nm = bayesian_deltaNV.mdl2.subject_id{iS};
    for iRun = 1:nRuns
        run_nm = ['run',num2str(iRun)];
        trial_idx = (1:nTrialsPerRun) + nTrialsPerRun*(iRun - 1);
        bayesian_deltaNV.mdl2.(run_nm) = NV_data_mdl2(trial_idx);
    end % run loop
end % subject loop

%% make one global save
save([gitPath, 'bayesian_deltaNV_data.mat'],'bayesian_deltaNV');