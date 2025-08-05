function[mdl_prm] = computational_mdl_prm(mdl_n)
% [mdl_prm] = computational_mdl_prm(mdl_n)
% script defining model features for each model number
% 
% INPUTS
% mdl_n: model number
%
% OUTPUTS
% mdl_prm: structure containing the different model parameters to be used

switch mdl_n
    case 1 % simple model (only kR and kP)
        F_prm_names = {};
        G_prm_names = {'kR','kP'};
        include_Inc = false;
        include_R = true;
        include_P = true;
        include_Ep = false;
        include_Em = false;
        include_Fp = false;
        include_currEff = false;
        include_prevEff = false;
        binary_answers = false;
        % positivity constraints
        pos.kR = true;
        pos.kP = true;
    case 2 % simple model (no time effect + no bias)
        F_prm_names = {};
        G_prm_names = {'kR','kP','kEp','kEm'};
        include_Inc = false;
        include_R = true;
        include_P = true;
        include_Ep = true;
        include_Em = true;
        include_Fp = false;
        include_currEff = false;
        include_prevEff = false;
        binary_answers = false;
        % positivity constraints
        pos.kR = true;
        pos.kP = true;
        pos.kEp = true;
        pos.kEm = true;
    case 3 % Arthur's model: include current efficiency instead of previous
        F_prm_names = {};
        G_prm_names = {'kR','kP','kEp','kEm','kBias','kFp','kLm'};
        include_Inc = false;
        include_R = true;
        include_P = true;
        include_Ep = true;
        include_Em = true;
        include_Fp = true;
        include_currEff = true;
        include_prevEff = false;
        binary_answers = false;
        % positivity constraints
        pos.kR = true;
        pos.kP = true;
        pos.kEp = true;
        pos.kEm = true;
        pos.kFp = true;
        pos.kLm = true;
        pos.kBias = false;
    case 4 % modified Arthur's model: binary output for choices (0/1) instead of 4-levels (0/0.25/0.75/1)
        F_prm_names = {};
        G_prm_names = {'kR','kP','kEp','kEm','kBias','kFp','kLm'};
        include_Inc = false;
        include_R = true;
        include_P = true;
        include_Ep = true;
        include_Em = true;
        include_Fp = true;
        include_currEff = true;
        include_prevEff = false;
        binary_answers = true;
        % positivity constraints
        pos.kR = true;
        pos.kP = true;
        pos.kEp = true;
        pos.kEm = true;
        pos.kFp = true;
        pos.kLm = true;
        pos.kBias = false;
    case 5 % main model (used in the 2 papers)
        F_prm_names = {};
        G_prm_names = {'kR','kP','kEp','kEm','kBias','kFp','kLm'};
        include_Inc = false;
        include_R = true;
        include_P = true;
        include_Ep = true;
        include_Em = true;
        include_Fp = true;
        include_currEff = false;
        include_prevEff = true;
        binary_answers = false;
        % positivity constraints
        pos.kR = true;
        pos.kP = true;
        pos.kEp = true;
        pos.kEm = true;
        pos.kFp = true;
        pos.kLm = true;
        pos.kBias = false;
    case 6 % model with binary output for choices (0/1) instead of 4-levels (0/0.25/0.75/1) (for applying oMCD)
        F_prm_names = {};
        G_prm_names = {'kR','kP','kEp','kEm','kBias','kFp','kLm'};
        include_Inc = false;
        include_R = true;
        include_P = true;
        include_Ep = true;
        include_Em = true;
        include_Fp = true;
        include_currEff = false;
        include_prevEff = true;
        binary_answers = true;
        % positivity constraints
        pos.kR = true;
        pos.kP = true;
        pos.kEp = true;
        pos.kEm = true;
        pos.kFp = true;
        pos.kLm = true;
        pos.kBias = false;
    case 7  % simple model like model 2 but including bias (no time effect)
        F_prm_names = {};
        G_prm_names = {'kR','kP','kEp','kEm','kBias'};
        include_Inc = false;
        include_R = true;
        include_P = true;
        include_Ep = true;
        include_Em = true;
        include_Fp = false;
        include_currEff = false;
        include_prevEff = false;
        binary_answers = false;
        % positivity constraints
        pos.kR = true;
        pos.kP = true;
        pos.kEp = true;
        pos.kEm = true;
        pos.kBias = false;
    otherwise
        error(['no model ',num2str(mdl_n)]);
end
% extract number of parameters
n_F_prm = length(F_prm_names);
n_G_prm = length(G_prm_names);

%% pool all parameters in output structure
% model number
mdl_prm.mdl_n = mdl_n;
% global number of parameters
mdl_prm.F_prm_names = F_prm_names;
mdl_prm.n_F_prm = n_F_prm;
mdl_prm.G_prm_names = G_prm_names;
mdl_prm.n_G_prm = n_G_prm;

% which parameters to include?
% monetary incentives
mdl_prm.include_Inc = include_Inc;
mdl_prm.include_R = include_R;
mdl_prm.include_P = include_P;
% physical E-specific parameters
mdl_prm.include_Ep = include_Ep;
mdl_prm.include_Fp = include_Fp; % physical fatigue
% mental E-specific parameters
mdl_prm.include_Em = include_Em;
mdl_prm.include_currEff = include_currEff; % current efficiency
mdl_prm.include_prevEff = include_prevEff; % previous efficiency
% choices: model as binary output (0=low E; 1=high E) or with confidence
% (0/0.25/0.75/1, i.e. lowE high Conf, low E low conf, high E low conf, high E high conf)
mdl_prm.binary_answers = binary_answers;
% positivity constraints on the parameters?
mdl_prm.pos = pos;

end % function