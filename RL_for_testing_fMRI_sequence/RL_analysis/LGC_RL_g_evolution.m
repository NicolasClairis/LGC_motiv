function [ gx ] = LGC_RL_g_evolution( x, P, u, inG )
%[ gx ] = LGC_RL_g_evolution( x, P, u, inG )
%
% INPUTS
% x: hidden states = Q.values
%
% P: parameters to estimate
%
% u: task parameters
%
% inG: indication about which model to use
%
% See also MS2_RL_Fabien_version & MS2_RL_f_evolution.m

%% load model parameters
side_bias   = inG.RL_mdl_prm.side_bias;
opt_bias    = inG.RL_mdl_prm.opt_bias;

iP = 1;
% inverse_temperature_prm = P(iP); % beta inverse temperature
temperature_prm = exp(P(iP)); % beta temperature
switch side_bias
    case 0
        side_bias_prm = 0;
    case 1
        iP = iP + 1;
        side_bias_prm = P(iP);
end
switch opt_bias
    case 'ntal'
        iP = iP + 1;
        opt_bias_N = P(iP);
    case 'all_pairs'
        iP = iP + 1;
        opt_bias_G = P(iP);
        iP = iP + 1;
        opt_bias_N = P(iP);
        iP = iP + 1;
        opt_bias_L = P(iP);
end

%% load hidden states
Q_GP_gain   = x(1);
Q_GP_ntal   = x(2);
Q_LP_ntal   = x(3);
Q_LP_loss   = x(4);

%% load task variables
% trialN             = u(1);
 pairValence        = u(2);
 goodSide           = u(3); % side of best option (1 right/-1 left)
%  lastPairValence    = u(4);
%  lastChoice         = u(5);
%  lastOutcome        = u(6);

%% determine DV for choice model
switch pairValence
    case 1
        DV = Q_GP_gain - Q_GP_ntal;
    case 0
        DV = 0; % both Q.values should be equal to zero always...
    case -1
        DV = Q_LP_ntal - Q_LP_loss;
end

%% motor bias congruent or not with best option
mBias = goodSide*side_bias_prm;

%% option bias
switch opt_bias
    case ''
        opt_bias_prm = 0;
    case 'ntal'
        switch pairValence
            case {-1,1} % no bias for gain and loss pairs
                opt_bias_prm = 0;
            case 0
                opt_bias_prm = opt_bias_N;
        end
    case 'all_pairs'
        switch pairValence
            case 1
                opt_bias_prm = opt_bias_G;
            case 0
                opt_bias_prm = opt_bias_N;
            case -1
                opt_bias_prm = opt_bias_L;
        end
end

%% output = p(choice = best option)
% gx = sigmo( (DV + mBias + opt_bias_prm), inverse_temperature_prm);
gx = sigmo( (DV + mBias + opt_bias_prm), (1/temperature_prm) );

end % function