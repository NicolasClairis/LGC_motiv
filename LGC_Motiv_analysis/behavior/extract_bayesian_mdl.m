function[NV_chosen, deltaNV_hE_min_lE, confidence, pChoice] = extract_bayesian_mdl(resultsFolder, subBehaviorFolder,...
    sub_nm, run_nm, task_fullName, mdl_nm)
% [NV_chosen, deltaNV_hE_min_lE, confidence, pChoice] = extract_bayesian_mdl(resultsFolder, sub_nm, run_nm, mdl_nm)
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
% mdl_nm: model name
%
% OUTPUTS
% NV_chosen: 1*nTrials vector with information about the net value of the 
% chosen option for the current study, subject and run.
%
% deltaNV_hE_min_lE: 1*nTrials vector with information about the net value of the 
% high effort option - low effort option for the current study, subject and run.
%
% confidence: 1*nTrials vector with information about the confidence
% inferred by the model for the current study, subject and run.
%
% pChoice: 1*nTrials vector with information about the probability of
% choosing the high effort option for the current study, subject and run.

%% load data
deltaNVstruct = getfield(load([resultsFolder,'bayesian_deltaNV_data.mat']),'bayesian_deltaNV');
pChoiceStruct = getfield(load([resultsFolder,'bayesian_pChoice_data.mat']),'bayesian_pChoice');
run_nm_bis = ['run',run_nm];
%% choice of low or high effort option?
[choice_highE] = extract_choice_hE(subBehaviorFolder, sub_nm, run_nm, task_fullName);
%% extract high effort level
deltaNV_hE_min_lE = deltaNVstruct.(mdl_nm).(['CID',sub_nm]).(run_nm_bis);
NV_chosen = deltaNV_hE_min_lE.*(choice_highE == 1) -deltaNV_hE_min_lE.*(choice_highE == 0);

%% confidence
pChoice = pChoiceStruct.(mdl_nm).(['CID',sub_nm]).(run_nm_bis);
confidence = (pChoice - 0.5).^2;
% normalize confidence to make it between 0 and 1
confidence = confidence./0.25;
end % function