function[NV_chosen_min_unch, deltaNV_hE_min_lE, confidence_highE, pChoice,...
    deltaNV_hE_min_lE_plus_bias, NV_ch_min_unch_with_bias] = extract_bayesian_mdl(resultsFolder, subBehaviorFolder,...
    sub_nm, run_nm, task_fullName, mdl_nm)
% [NV_chosen_min_unch, deltaNV_hE_min_lE, confidence_highE, pChoice,...
%     deltaNV_hE_min_lE_plus_bias, NV_ch_min_unch_with_bias] = extract_bayesian_mdl(resultsFolder, subBehaviorFolder,...
%     sub_nm, run_nm, task_fullName, mdl_nm)
% extract_bayesian_mdl will extract the different variables inferred thanks
% to the bayesian modeling approach for the model and subjects and run
% defined in the inputs.
%
% INPUTS
% resultsFolder: folder where bayesian data is stored
%
% subBehaviorFolder: folder where subject data is stored
%
% sub_nm: string with subject CID name
%
% run_nm: string with run name
%
% task_fullName: task full name 'mental'/'physical'
%
% mdl_nm: model name (should be a string with the form 'mdl_X' with X being
% the model number)
%
% OUTPUTS
% NV_chosen_min_unch: 1*nTrials vector with information about the net value of the 
% chosen - unchosen option for the current study, subject and run.
%
% deltaNV_hE_min_lE: 1*nTrials vector with information about the net value of the 
% high effort option - low effort option for the current study, subject and run.
%
% confidence_highE: 1*nTrials vector with information about the confidence
% inferred by the model for the current study, subject and run. Confidence
% is computed in reference of the high effort option as [p(high E)-0.5]²
%
% pChoice: 1*nTrials vector with information about the probability of
% choosing the high effort option for the current study, subject and run.
%
% deltaNV_hE_min_lE_plus_bias: 1*nTrials vector with deltaNV_hE_min_lE + bias 
% (ie everything used for p(choice) before softmax transformation.
%
% NV_ch_min_unch_with_bias: 1*nTrials vector with information about the net value of the 
% chosen - unchosen option for the current study, subject and run. Same as
% NV_chosen_min_unch, but also including bias (ie same as pChoice before
% the sigmoid transformation)

%% load data
mdlN = strrep(mdl_nm,'mdl_','');
model_data_struct = load([resultsFolder,'bayesian_model_',mdlN,'_results.mat']);

%% choice of low or high effort option?
[choice_highE] = extract_choice_hE(subBehaviorFolder, sub_nm, run_nm, task_fullName);

%% load model bias (if present in the model)
if ismember('kBias',fieldnames(model_data_struct.prm))
    kBias = prm.kBias;
else
    kBias = 0;
end

%% identify relevant subject and trials
sub_nm_bis = ['CID',sub_nm];
run_nm_bis = ['run',run_nm];

% verify fields exist
% subject
if ~ismember(sub_nm_bis, fieldnames(model_data_struct.dV_pred_perSub_perRun))
    error([sub_nm_bis,' missing for bayesian model ',mdlN,' extraction']);
end
% run
if ~ismember(run_nm_bis, fieldnames(model_data_struct.dV_pred_perSub_perRun.(sub_nm_bis)))
    error([run_nm_bis,' missing for ',sub_nm_bis,' with bayesian model ',mdlN]);
end

%% extract net value for current subject and run
deltaNV_hE_min_lE = model_data_struct.dV_pred_perSub_perRun.(sub_nm_bis).(run_nm_bis);

% in case of infinite values, replace with values giving a high probability
% but not equal to +Inf or -Inf
plus_inf_idx = find(deltaNV_hE_min_lE == Inf); % p(choice=hE) = 1
minus_inf_idx = find(deltaNV_hE_min_lE == -Inf); % p(choice=hE) = 0
if ~isempty(plus_inf_idx) || ~isempty(minus_inf_idx)
    
    % replace Inf value by value giving similar probability
    if ~isempty(plus_inf_idx) % p(choice=hE) = 1
        pc_high = 0.9999999999999999;
        deltaNV_hE_min_lE(plus_inf_idx) = -log((1-pc_high)/pc_high) - kBias;
    end
    if ~isempty(minus_inf_idx) % p(choice=hE) = 0
        pc_low = 0.0000000000000001;
        deltaNV_hE_min_lE(plus_inf_idx) = -log((1-pc_low)/pc_low) - kBias;
    end
end % filter Inf values
% extract net value for chosen option
choice_vector = (choice_highE' == 1) - (choice_highE' == 0);
NV_chosen_min_unch = deltaNV_hE_min_lE.*choice_vector;
NV_chosen_min_unch(isnan(choice_highE)) = NaN;

%% delta NV + bias (ie everything used to compute p(choice) before the softmax transformation)
deltaNV_hE_min_lE_plus_bias = deltaNV_hE_min_lE + kBias;
NV_ch_min_unch_with_bias = deltaNV_hE_min_lE_plus_bias.*choice_vector;

%% extract p(choice) predicted by the model
pChoice = model_data_struct.choices_pred_perSub_perRun.(sub_nm_bis).(run_nm_bis);

%% confidence
confidence_highE = (pChoice - 0.5).^2;
% normalize confidence to make it between 0 and 1
confidence_highE = confidence_highE./0.25;
end % function