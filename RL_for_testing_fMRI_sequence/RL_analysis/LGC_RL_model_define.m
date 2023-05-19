function[RL_mdl_prm] = LGC_RL_model_define(RL_model_n)
%[RL_mdl_prm] = LGC_RL_model_define(RL_model_n)
% MS2_RL_model_define sets the parameters to use depending on the model
% number
%
% INPUTS
% RL_model_n: model number
%
% OUTPUTS
% RL_mdl_prm: structure with model parameters
%       .alpha_prm
%           'one_learningRate': one single learning rate
%           'GL_learningRates': one learning rate per condition
%           (gain/loss)
%
%       .RP_weights
%           (0): no difference in the weighting between gains and losses
%           update: Q(t+1) = Q(t) + a*[outcome - Q(t)] where
%           (1): Q(t+1) = Q(t) + a*[R/P*outcome - Q(t)] where R/P is
%           different according to valence
%
%       .sigmaQ_prm: value for prior on initial Q.values
%           (0): start with no SD
%           (X): start with X as SD for each hidden state
%
%       .Qvalues_relation
%           'indpdt_Q': independent Q.values
%           'counterfactual_Q': for any given pair, Q.values evolve in
%           parallel as if the feedback for both options had been seen
%           'counterfactual_Q_bis': for any given pair, Qunchosen is always
%           equal to -Qchosen
%
%       .fbk_relative_coding
%           'no': feedback = +1 for gain (+10€)/ 0 for neutral/-1 for loss
%           (-10€) = 3 feedbacks = -1/0/+1
%           'yes': feedback = +1 for gain (+10€) in gain pair or neutral
%           (0€) in loss pair/-1 for neutral (0€) in gain pair or loss
%           (-10€) in loss pair = only feedbacks = +1/-1
%
%       .side_bias
%           (0): no side bias included
%           (1): add a side bias
%
%       .opt_bias: option bias = bias towards one option, independent of
%       the learned value
%           '': no bias
%           'all_pairs': for all pairs: might interfere with Q.values
%           learning...
%           'ntal': for neutral pair only (since Q.values should not change
%           and remain at zero)
% See also MS2_RL_model_bis.m

switch RL_model_n
    case 1
        RL_mdl_prm.alpha_prm = 'one_learningRate';
        RL_mdl_prm.sigmaQ_prm = 0.5;
        RL_mdl_prm.Qvalues_relation = 'counterfactual_Q';
        RL_mdl_prm.RP_weights = 0;
        RL_mdl_prm.side_bias = 0;
        RL_mdl_prm.opt_bias = 'all_pairs';
        RL_mdl_prm.fbk_relative_coding = 'no';
    case 2
        RL_mdl_prm.alpha_prm = 'one_learningRate';
        RL_mdl_prm.sigmaQ_prm = 0.5;
        RL_mdl_prm.Qvalues_relation = 'indpdt_Q';
        RL_mdl_prm.RP_weights = 0;
        RL_mdl_prm.side_bias = 0;
        RL_mdl_prm.opt_bias = 'ntal';
        RL_mdl_prm.fbk_relative_coding = 'no';
    case 3 % like model 1 but option bias only for neutral pair
        RL_mdl_prm.alpha_prm = 'one_learningRate';
        RL_mdl_prm.sigmaQ_prm = 0.5;
        RL_mdl_prm.Qvalues_relation = 'counterfactual_Q';
        RL_mdl_prm.RP_weights = 0;
        RL_mdl_prm.side_bias = 0;
        RL_mdl_prm.opt_bias = 'ntal';
        RL_mdl_prm.fbk_relative_coding = 'no';
    case 4 % 2 learning rates alphaGain and alphaLoss
        RL_mdl_prm.alpha_prm = 'GL_learningRates';
        RL_mdl_prm.sigmaQ_prm = 0.5;
        RL_mdl_prm.Qvalues_relation = 'counterfactual_Q';
        RL_mdl_prm.RP_weights = 0;
        RL_mdl_prm.side_bias = 0;
        RL_mdl_prm.opt_bias = 'ntal';
        RL_mdl_prm.fbk_relative_coding = 'no';
    case 5
        RL_mdl_prm.alpha_prm = 'GL_learningRates';
        RL_mdl_prm.sigmaQ_prm = 0.5;
        RL_mdl_prm.Qvalues_relation = 'counterfactual_Q';
        RL_mdl_prm.RP_weights = 0;
        RL_mdl_prm.side_bias = 0;
        RL_mdl_prm.opt_bias = 'all_pairs';
        RL_mdl_prm.fbk_relative_coding = 'no';
    case 6 % like model 3 but initial SD = 0
        RL_mdl_prm.alpha_prm = 'one_learningRate';
        RL_mdl_prm.sigmaQ_prm = 0;
        RL_mdl_prm.Qvalues_relation = 'counterfactual_Q';
        RL_mdl_prm.RP_weights = 0;
        RL_mdl_prm.side_bias = 0;
        RL_mdl_prm.opt_bias = 'ntal';
        RL_mdl_prm.fbk_relative_coding = 'no';
    case 7
        RL_mdl_prm.alpha_prm = 'one_learningRate';
        RL_mdl_prm.sigmaQ_prm = 0.5;
        RL_mdl_prm.Qvalues_relation = 'counterfactual_Q';
        RL_mdl_prm.RP_weights = 0;
        RL_mdl_prm.side_bias = 0;
        RL_mdl_prm.opt_bias = 'ntal';
        RL_mdl_prm.fbk_relative_coding = 'no';
    case 8
        RL_mdl_prm.alpha_prm = 'one_learningRate';
        RL_mdl_prm.sigmaQ_prm = 0;
        RL_mdl_prm.Qvalues_relation = 'counterfactual_Q';
        RL_mdl_prm.RP_weights = 0;
        RL_mdl_prm.side_bias = 0;
        RL_mdl_prm.opt_bias = 'all_pairs';
        RL_mdl_prm.fbk_relative_coding = 'no';
    case 9
        RL_mdl_prm.alpha_prm = 'one_learningRate';
        RL_mdl_prm.sigmaQ_prm = 0;
        RL_mdl_prm.Qvalues_relation = 'counterfactual_Q_bis';
        RL_mdl_prm.RP_weights = 0;
        RL_mdl_prm.side_bias = 0;
        RL_mdl_prm.opt_bias = 'ntal';
        RL_mdl_prm.fbk_relative_coding = 'no';
    case 10 % like model 6 but independent Q.values
        RL_mdl_prm.alpha_prm = 'one_learningRate';
        RL_mdl_prm.sigmaQ_prm = 0;
        RL_mdl_prm.Qvalues_relation = 'indpdt_Q';
        RL_mdl_prm.RP_weights = 0;
        RL_mdl_prm.side_bias = 0;
        RL_mdl_prm.opt_bias = 'ntal';
        RL_mdl_prm.fbk_relative_coding = 'no';
    case 11 % like model 6 but 2 different learning rates (gain/loss)
        RL_mdl_prm.alpha_prm = 'GL_learningRates';
        RL_mdl_prm.sigmaQ_prm = 0;
        RL_mdl_prm.Qvalues_relation = 'counterfactual_Q';
        RL_mdl_prm.RP_weights = 0;
        RL_mdl_prm.side_bias = 0;
        RL_mdl_prm.opt_bias = 'ntal';
        RL_mdl_prm.fbk_relative_coding = 'no';
    case 12 % like model 6 but with relative coding
        RL_mdl_prm.alpha_prm = 'one_learningRate';
        RL_mdl_prm.sigmaQ_prm = 0;
        RL_mdl_prm.Qvalues_relation = 'counterfactual_Q';
        RL_mdl_prm.RP_weights = 0;
        RL_mdl_prm.side_bias = 0;
        RL_mdl_prm.opt_bias = 'ntal';
        RL_mdl_prm.fbk_relative_coding = 'yes';
end

end % function