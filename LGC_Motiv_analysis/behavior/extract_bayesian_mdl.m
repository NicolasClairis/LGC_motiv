function[NV_chosen, deltaNV_hE_min_lE, confidence_highE, pChoice] = extract_bayesian_mdl(resultsFolder, subBehaviorFolder,...
    sub_nm, run_nm, task_fullName, mdl_nm)
% [NV_chosen, deltaNV_hE_min_lE, confidence_highE, pChoice] = extract_bayesian_mdl(resultsFolder, sub_nm, run_nm, mdl_nm)
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
% NV_chosen: 1*nTrials vector with information about the net value of the 
% chosen option for the current study, subject and run.
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

study_nm = 'study1';
%% load data
deltaNVstruct = getfield(load([resultsFolder,'bayesian_deltaNV_data.mat']),'bayesian_deltaNV');
pChoiceStruct = getfield(load([resultsFolder,'bayesian_pChoice_data.mat']),'bayesian_pChoice');
run_nm_bis = ['run',run_nm];
%% choice of low or high effort option?
[choice_highE] = extract_choice_hE(subBehaviorFolder, sub_nm, run_nm, task_fullName);
%% extract net value
deltaNV_hE_min_lE = deltaNVstruct.(mdl_nm).(['CID',sub_nm]).(run_nm_bis);
% in case of infinite values, replace with values giving a high probability
% but not equal to +Inf or -Inf
plus_inf_idx = find(deltaNV_hE_min_lE == Inf); % p(choice=hE) = 1
minus_inf_idx = find(deltaNV_hE_min_lE == -Inf); % p(choice=hE) = 0
if ~isempty(plus_inf_idx) || ~isempty(minus_inf_idx)
    % load model bias (if present in the model)
    mdl_nm_start = strfind(mdl_nm,'_') + 1;
    mdlN = mdl_nm(mdl_nm_start:end);
    prm = prm_extraction(study_nm, {sub_nm},...
        'bayesian', mdlN);
    if ~isempty(strcmp(fieldnames(prm),'kBiasM'))
        kBias = prm.kBiasM;
    else
        kBias = 0;
    end
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
NV_chosen = deltaNV_hE_min_lE.*(choice_highE' == 1) +...
    -deltaNV_hE_min_lE.*(choice_highE' == 0);
NV_chosen(isnan(choice_highE)) = NaN;

%% confidence
pChoice = pChoiceStruct.(mdl_nm).(['CID',sub_nm]).(run_nm_bis);
confidence_highE = (pChoice - 0.5).^2;
% normalize confidence to make it between 0 and 1
confidence_highE = confidence_highE./0.25;
end % function