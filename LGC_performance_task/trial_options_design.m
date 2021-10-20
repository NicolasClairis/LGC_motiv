function[choiceOptions] = trial_options_design(n_R_levels, n_E_levels, punishment_yn, nTrials)
% [choiceOptions] = trial_options_design(n_R_levels, n_E_levels, punishment_yn, nTrials)
% trial_options_design will prepare a design matrix repeating small
% mini-blocks
%
% INPUTS
%
% OUTPUTS
% choiceOptions: structure with the design matrix

nOptions = n_R_levels*n_E_levels + n_R_levels*n_E_levels*strcmp(punishment_yn,'yes');

%% check that the number of trials is compatible with the number of options
if mod(nTrials, nOptions) ~= 0
    error(['number of trials = ',num2str(nTrials),' while number of possible options = ',num2str(nOptions),'. Something needs to be fixed.']);
end

nMiniBlocksRepeats = nTrials/nOptions;
R_options = [1, 1, 1, 2, 2, 2, 3, 3, 3];
E_options = [1, 2, 3, 1, 2, 3, 1, 2, 3];
% add punishment options eventually
if strcmp(punishment_yn,'yes')
    P_options = -R_options;
    moneyOptions = [R_options, P_options];
    E_options = [E_options, E_options];
else
    moneyOptions = R_options;
end

if length(moneyOptions) ~= nOptions || length(E_options) ~= nOptions
   error('problem in number of options somewhere');
end

%% prepare mini-blocks
[choiceOptions.moneyLevel,...
    choiceOptions.effortLevel] = deal(NaN(1, nTrials));

for iMiniB = 1:nMiniBlocksRepeats
    % randomize the order of the options within a mini-block
    miniBlockOrderTmp = randperm(nOptions);
    % extract the index for the current mini-block
    trial_idx = (1:nOptions) + nOptions*(iMiniB - 1);
    
    choiceOptions.moneyLevel(trial_idx)    = moneyOptions(miniBlockOrderTmp);
    choiceOptions.effortLevel(trial_idx)   = E_options(miniBlockOrderTmp);
end % mini-blocks for the design matrix

end % function end