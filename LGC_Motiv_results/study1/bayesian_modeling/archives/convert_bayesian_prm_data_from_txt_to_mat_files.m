%% script to convert subject_id data from Arthur's txt files into cell


% model 1 with bias
bayesian_mdl.mdl_1.kR = model_raw.mdl_1(1,:);
bayesian_mdl.mdl_1.kP = model_raw.mdl_1(2,:);
bayesian_mdl.mdl_1.kEp = model_raw.mdl_1(3,:);
bayesian_mdl.mdl_1.kEm = model_raw.mdl_1(4,:);
bayesian_mdl.mdl_1.kBiasM = model_raw.mdl_1(5,:);
bayesian_mdl.mdl_1.kFp = model_raw.mdl_1(6,:);
bayesian_mdl.mdl_1.kLm = model_raw.mdl_1(7,:);
lSub_f_mdl.mdl_1 = 8;

% model 2 (no bias)
bayesian_mdl.mdl_2.kR = model_raw.mdl_2(1,:);
bayesian_mdl.mdl_2.kP = model_raw.mdl_2(2,:);
bayesian_mdl.mdl_2.kEp = model_raw.mdl_2(3,:);
bayesian_mdl.mdl_2.kEm = model_raw.mdl_2(4,:);
bayesian_mdl.mdl_2.kFp = model_raw.mdl_2(5,:);
bayesian_mdl.mdl_2.kLm = model_raw.mdl_2(6,:);
lSub_f_mdl.mdl_2 = 7;

% model 3 (bias)
bayesian_mdl.mdl_3.kR = model_raw.mdl_3(2,:);
bayesian_mdl.mdl_3.kP = model_raw.mdl_3(3,:);
bayesian_mdl.mdl_3.kEp = model_raw.mdl_3(4,:);
bayesian_mdl.mdl_3.kEm = model_raw.mdl_3(5,:);
bayesian_mdl.mdl_3.kBiasM = model_raw.mdl_3(6,:);
bayesian_mdl.mdl_3.kFp = model_raw.mdl_3(7,:);
bayesian_mdl.mdl_3.kLm = model_raw.mdl_3(8,:);
lSub_f_mdl.mdl_3 = 1;

% model 4 (bias + no confidence in the model)
bayesian_mdl.mdl_4.kR = model_raw.mdl_4(2,:);
bayesian_mdl.mdl_4.kP = model_raw.mdl_4(3,:);
bayesian_mdl.mdl_4.kEp = model_raw.mdl_4(4,:);
bayesian_mdl.mdl_4.kEm = model_raw.mdl_4(5,:);
bayesian_mdl.mdl_4.kBiasM = model_raw.mdl_4(6,:);
bayesian_mdl.mdl_4.kFp = model_raw.mdl_4(7,:);
bayesian_mdl.mdl_4.kLm = model_raw.mdl_4(8,:);
lSub_f_mdl.mdl_4 = 1;

% define the model data
model_n = input('Which model number?');
mdl_nm = ['mdl_',num2str(model_n)];
lSub = lSub_f_mdl.(mdl_nm);
model_data = model_raw.(mdl_nm)(lSub,:);
bayesian_mdl.(mdl_nm).subject_id = convert_sub_id_from_num_to_cell(model_data);

%% save the data
save('behavioral_prm_tmp.mat','bayesian_mdl','model_raw');