function [ fx ] = LGC_RL_f_evolution( x, P, u, inF )
%[ fx ] = LGC_RL_f_evolution( x, P, u, inF )
%
% INPUTS
% x: hidden states = Q.values
%
% P: parameters to estimate
%
% u: task parameters
%
% inF: indication about which model to use
%
% See also MS2_RL_Fabien_version & MS2_RL_g_evolution.m

%% load model parameters
alpha_prm = inF.RL_mdl_prm.alpha_prm;
Qvalues_relation = inF.RL_mdl_prm.Qvalues_relation;
RP_weights = inF.RL_mdl_prm.RP_weights;
% side_bias = inF.RL_mdl_prm.side_bias;
fbk_relative_coding = inF.RL_mdl_prm.fbk_relative_coding;

%% load Q.values
Q_GP_gain   = x(1);
Q_GP_ntal   = x(2);
Q_LP_ntal   = x(3);
Q_LP_loss   = x(4);

%% load parameters to estimate
iP = 1;
switch alpha_prm
    case 'one_learningRate' % same learning rate for gain and loss
        alpha = P(iP);
        alpha_G = alpha;
        alpha_L = alpha;
    case 'GL_learningRates' % different learning rate for gain and for loss
        alpha_G = P(iP);
        iP = iP + 1;
        alpha_L = P(iP);
end

switch RP_weights
    case 0 % no scaling
        R_weight = 1;
        L_weight = 1;
    case 1
        iP = iP + 1;
        R_weight = P(iP);
        iP = iP + 1;
        L_weight = P(iP);
end

%% load task variables
% trialN             = u(1);
%  pairValence        = u(2);
%  goodSide           = u(3);
lastPairValence    = u(4);
lastChoice         = u(5);
lastOutcome        = u(6);

% adapt feedback values if relative coding
if strcmp(fbk_relative_coding,'yes')
    switch lastPairValence % transform gain neutral and loss neutral for relative coding
        case 1 % gain pair
            if lastOutcome == 0 % neutral => like a loss (relative coding)
                lastOutcome = -1;
            end
        case -1 % loss pair
            if lastOutcome == 0 % neutral => like a gain (relative coding)
                lastOutcome = 1;
            end
    end
end

% which option needs updating for the current trial?
switch lastPairValence
    case 1 % gain
        switch lastChoice
            case 0 % gain neutral
                cue_to_update = 'GP_ntal';
            case 1 % gain gain
                cue_to_update = 'GP_gain';
        end
    case 0 % neutral
        switch lastChoice
            case 0 % neutral A
                cue_to_update = 'NP_A';
            case 1 % neutral B
                cue_to_update = 'NP_B';
        end
    case -1 % loss
        switch lastChoice
            case 0 % loss loss
                cue_to_update = 'LP_loss';
            case 1 % loss neutral
                cue_to_update = 'LP_ntal';
        end
end

%% update Q.values according to the model used
switch Qvalues_relation
    case 'counterfactual_Q' % update Q.values for both options
        
        switch fbk_relative_coding
            case 'no'
                GP_counterfactual = 1;
                LP_counterfactual = -1;
            case 'yes'
                GP_counterfactual = 0;
                LP_counterfactual = 0;
        end
        
        switch cue_to_update
            case 'GP_ntal'
                Q_GP_ntal = Q_GP_ntal + alpha_G*(lastOutcome*R_weight - Q_GP_ntal);
                Q_GP_gain = Q_GP_gain + alpha_G*( (GP_counterfactual - lastOutcome)*R_weight - Q_GP_gain); %as if received feedback +1 for gain option when neutral option outcome is 0, and reciprocally
            case 'GP_gain'
                Q_GP_gain = Q_GP_gain + alpha_G*(lastOutcome*R_weight - Q_GP_gain);
                Q_GP_ntal = Q_GP_ntal + alpha_G*( (GP_counterfactual - lastOutcome)*R_weight - Q_GP_ntal);
            case 'LP_loss'
                Q_LP_loss = Q_LP_loss + alpha_L*(lastOutcome*L_weight - Q_LP_loss);
                Q_LP_ntal = Q_LP_ntal + alpha_L*( (LP_counterfactual - lastOutcome)*L_weight - Q_LP_ntal);% as if received feedback -1 for loss option when neutral option outcome is 0 and reciprocally
            case 'LP_ntal'
                Q_LP_ntal = Q_LP_ntal + alpha_L*(lastOutcome*L_weight - Q_LP_ntal);
                Q_LP_loss = Q_LP_loss + alpha_L*( (LP_counterfactual - lastOutcome)*L_weight - Q_LP_loss);
        end
    case 'indpdt_Q' % update only the chosen option which received a feedback
        switch cue_to_update
            case 'GP_ntal'
                Q_GP_ntal = Q_GP_ntal + alpha_G*(lastOutcome*R_weight - Q_GP_ntal);
            case 'GP_gain'
                Q_GP_gain = Q_GP_gain + alpha_G*(lastOutcome*R_weight - Q_GP_gain);
            case 'LP_loss'
                Q_LP_loss = Q_LP_loss + alpha_L*(lastOutcome*L_weight - Q_LP_loss);
            case 'LP_ntal'
                Q_LP_ntal = Q_LP_ntal + alpha_L*(lastOutcome*L_weight - Q_LP_ntal);
        end
    case 'counterfactual_Q_bis' % Qunchosen = -Qchosen always
        switch cue_to_update
            case 'GP_ntal'
                Q_GP_ntal = Q_GP_ntal + alpha_G*(lastOutcome*R_weight - Q_GP_ntal);
                Q_GP_gain = -Q_GP_ntal; %as if received feedback +1 for gain option when neutral option outcome is 0, and reciprocally
            case 'GP_gain'
                Q_GP_gain = Q_GP_gain + alpha_G*(lastOutcome*R_weight - Q_GP_gain);
                Q_GP_ntal = -Q_GP_gain;
            case 'LP_loss'
                Q_LP_loss = Q_LP_loss + alpha_L*(lastOutcome*L_weight - Q_LP_loss);
                Q_LP_ntal = -Q_LP_loss;% as if received feedback -1 for loss option when neutral option outcome is 0 and reciprocally
            case 'LP_ntal'
                Q_LP_ntal = Q_LP_ntal + alpha_L*(lastOutcome*L_weight - Q_LP_ntal);
                Q_LP_loss = -Q_LP_ntal;
        end
end

%% output
fx = [Q_GP_gain, Q_GP_ntal, Q_LP_ntal, Q_LP_loss];


end % function