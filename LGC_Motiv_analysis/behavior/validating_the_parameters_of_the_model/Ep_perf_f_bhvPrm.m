% performance in physical task (latency, overshoot, peak force, etc.)
% depending on behavioral parameters (kEp, kFp) for the selected subjects.

%% subject selection
study_nm = 'study1';
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% working directories
computerRoot = LGCM_root_paths;
dataRoot = [computerRoot, filesep, study_nm, filesep];

%% load parameters
[mdlType, mdlN] = behavioral_model_selection;
prm = prm_extraction(study_nm, subject_id, mdlType, mdlN);
% extract parameter of interest
prmToTest = {'kFp','kEp'};
nPrmToTest = length(prmToTest);

%% figure
pSize = 30;
lWidth = 3;
grey = [143 143 143]./255;


for iPrm = 1:nPrmToTest
    prm_nm = prmToTest{iPrm};
    prm_of_interest = prm.(prm_nm);
    
end % parameter loop