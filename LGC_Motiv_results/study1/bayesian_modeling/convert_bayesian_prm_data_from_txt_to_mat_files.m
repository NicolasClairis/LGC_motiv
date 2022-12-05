%% script to convert subject_id data from Arthur's txt files into cell


% model 1 with bias
bayesian_mdl.mdl1.kR = model_raw.mdl1(1,:);
bayesian_mdl.mdl1.kP = model_raw.mdl1(2,:);
bayesian_mdl.mdl1.kEp = model_raw.mdl1(3,:);
bayesian_mdl.mdl1.kEm = model_raw.mdl1(4,:);
bayesian_mdl.mdl1.kBiasM = model_raw.mdl1(5,:);
bayesian_mdl.mdl1.kFp = model_raw.mdl1(6,:);
bayesian_mdl.mdl1.kLm = model_raw.mdl1(7,:);
lSub = 8;
% model 2 (no bias)
bayesian_mdl.mdl2.kR = model_raw.mdl2(1,:);
bayesian_mdl.mdl2.kP = model_raw.mdl2(2,:);
bayesian_mdl.mdl2.kEp = model_raw.mdl2(3,:);
bayesian_mdl.mdl2.kEm = model_raw.mdl2(4,:);
bayesian_mdl.mdl2.kFp = model_raw.mdl2(5,:);
bayesian_mdl.mdl2.kLm = model_raw.mdl2(6,:);
lSub = 7;

% model 3 (bias)
bayesian_mdl.mdl3.kR = model_raw.mdl3(2,:);
bayesian_mdl.mdl3.kP = model_raw.mdl3(3,:);
bayesian_mdl.mdl3.kEp = model_raw.mdl3(4,:);
bayesian_mdl.mdl3.kEm = model_raw.mdl3(5,:);
bayesian_mdl.mdl3.kBiasM = model_raw.mdl3(6,:);
bayesian_mdl.mdl3.kFp = model_raw.mdl3(7,:);
bayesian_mdl.mdl3.kLm = model_raw.mdl3(8,:);
lSub = 1;

lSub = input('Which line for subject_id?');
% define the model data
model_n = input('Which model number?');
mdl_nm = ['mdl',num2str(model_n)];
model_data = model_raw.(mdl_nm)(lSub,:);
bayesian_mdl.(mdl_nm).subject_id = convert_sub_id_from_num_to_cell(model_data);

%% save the data
save('behavioral_prm_tmp.mat','bayesian_mdl','model_raw');