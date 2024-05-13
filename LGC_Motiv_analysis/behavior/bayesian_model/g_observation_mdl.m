function [gx] = g_observation_mdl( x, phi, var, inG )
% [gx] = g_observation_mdl( x, phi, var, inG )
%
% INPUTS
% x: hidden state (not used currently)
%
% phi: phi observation parameters
%
% var: 1*nParameters vector with inputs to be used for the function
%
% inG: structure with relevant information to know what parameters to use
% for the current model
%
% OUTPUTS
% gx: predicted output for the VBA_NL_StateSpaceModel.m function


%% load all input variables
% Difference of incentive between the two choices
deltaR = var(strcmp(inG.var_names,'dR'));
deltaP = var(strcmp(inG.var_names,'dP'));
% Difference of effort between the two choices
deltaE = var(strcmp(inG.var_names,'dE'));
% Are we in a physical (1) or mental (0) block
Ep_or_Em_trials = var(strcmp(inG.var_names,'EpEm'));
Ep_trial = Ep_or_Em_trials == 1;
Em_trial = Ep_or_Em_trials == 0;
% physical fatigue
Fp = var(strcmp(inG.var_names,'Fp'));
% current trial or previous trial efficiency in solving mental effort
currEff = var(strcmp(inG.var_names,'currEff'));
prevEff = var(strcmp(inG.var_names,'prevEff'));

%% load model parameters
mdl_prm = inG.mdl_prm;
mdl_n = mdl_prm.mdl_n;
pos = mdl_prm.pos;

%% load parameters
switch mdl_n
    case 1
        [kR] = prm_input(phi(1), pos, 'kR'); % reward sensitivity
        [kP] = prm_input(phi(2), pos, 'kP'); % punishment sensitivity
    case 2
        [kR] = prm_input(phi(1), pos, 'kR'); % reward sensitivity
        [kP] = prm_input(phi(2), pos, 'kP'); % punishment sensitivity
        [kEp] = prm_input(phi(3), pos, 'kEp'); % physical effort sensitivity
        [kEm] = prm_input(phi(4), pos, 'kEm'); % mental effort sensitivity
    case {3,4,5,6} % 7 parameters: kR, kP, kEp, kEm, kBias, kFp, kLm
        [kR] = prm_input(phi(1), pos, 'kR'); % reward sensitivity
        [kP] = prm_input(phi(2), pos, 'kP'); % punishment sensitivity
        [kEp] = prm_input(phi(3), pos, 'kEp'); % physical effort sensitivity
        [kEm] = prm_input(phi(4), pos, 'kEm'); % mental effort sensitivity
        [kBias] = prm_input(phi(5), pos, 'kBias'); % bias for high effort
        [kFp] = prm_input(phi(6), pos, 'kFp'); % physical fatigue sensitivity
        [kLm] = prm_input(phi(7), pos, 'kLm'); % mental learning sensitivity
    otherwise
        error(['model ',num2str(mdl_n),' not ready yet']);
end

%% Subjective value computation
switch mdl_n
    case 1
        dV = kR.*deltaR + kP.*deltaP;
    case 2
        dV = (kR.*deltaR + kP.*deltaP) - deltaE.*(Ep_trial.*kEp + Em_trial.*kEm);
    case {3,4}
        dV = (kR.*deltaR + kP.*deltaP) - deltaE.*(Ep_trial.*(kEp + kFp.*Fp) + Em_trial.*(kEm - kLm.*currEff));
    case {5,6}
        dV = (kR.*deltaR + kP.*deltaP) - deltaE.*(Ep_trial.*(kEp + kFp.*Fp) + Em_trial.*(kEm - kLm.*prevEff));
    otherwise
        error(['model ',num2str(mdl_n),' not ready yet']);
end

%% sigmoid transformation
switch mdl_n
    case {1,2} % no bias
        gx = 1./(1 + exp(-dV));
    case {3,4,5,6} % include the bias
        gx = 1./(1 + exp(-(dV + kBias)));
    otherwise
        error(['model ',num2str(mdl_n),' not ready yet']);
end

end

function[prm] = prm_input(prm_raw, pos, prm_nm)
% function to transform parameter with positivity constraint when necessary
%
% INPUTS
% prm_raw: raw value of parameter (before any transformation)
%
% pos: structure containing a binary variable indicating if there is a positivity constraint
% (true) or not (false) for each parameter
%
% prm_nm: name of the current parameter to check if there is a positivity
% constraint for this parameter or not
%
% OUTPUTS
% prm: parameter after transformation with or without positivity constraint

switch pos.(prm_nm)
    case false % no transformation
        prm = prm_raw;
    case true % positivity constraint transformation
        prm = log(1+exp(prm_raw));
end

end % subfunction