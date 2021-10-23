function[trialOptions] = design_trialOptions(n_R_levels, n_E_levels, punishment_yn, nTrials)
% 
% design_trialOptions will create a potential design matrix for the version
% with a default option.
% 
% INPUTS
% n_R_levels: number of money levels
%
% n_E_levels: number of effort levels
%
% punishment_yn:
% 'yes': add punishments
% 'no': only rewards
%
% nTrials: number of trials in total
%
% OUTPUTS
% trialOptions: structure with information for each trial


%% define potential options for 1 mini-block
switch punishment_yn
    case 'no'
        options = 1:(n_R_levels*n_Elevels);
        [R_optionsPerBlock,...
            E_optionsPerBlock] = deal(NaN(1,(n_R_levels*n_E_levels)));
        jOption = 0;
        for iRoption = 1:n_R_levels
            for iEoption = 1:n_E_levels
                jOption = jOption + 1;
                R_optionsPerBlock(jOption) = iRoption;
                E_optionsPerBlock(jOption) = iEoption;
            end % effort loop
        end % reward loop
    case 'yes'
        options = 1:((n_R_levels*n_E_levels)*2);
        [R_optionsPerBlock,...
            E_optionsPerBlock] = deal(NaN(1,(n_R_levels*n_E_levels)*2));
        jOption = 0;
        for iRoption = [(-n_R_levels):(-1), 1:n_R_levels]
            for iEoption = 1:n_E_levels
                jOption = jOption + 1;
                R_optionsPerBlock(jOption) = iRoption;
                E_optionsPerBlock(jOption) = iEoption;
            end % effort loop
        end % reward loop
end
nOptions = length(options);

%% prepare the blocks
nBlocks = nTrials/nOptions;
if floor(nBlocks) < nBlocks
    error(['problem in the number of trials: you have ',...
        num2str(nTrials),' trials while there are ',...
        num2str(nOptions),' options. It cannot work like that.']);
end

R_options = repmat(R_optionsPerBlock,1,nBlocks);
E_options = repmat(E_optionsPerBlock,1,nBlocks);

%% randomize the order within each block
for iBlock = 1:nBlocks
    block_trials_idx = (1:nOptions) + nOptions*(iBlock - 1);
    block_rdm = randperm(nOptions);
    R_options(block_trials_idx) = R_options( block_trials_idx(block_rdm) );
    E_options(block_trials_idx) = E_options( block_trials_idx(block_rdm) );
end % loop on blocks

%% randomize the left/right side of the default option within each block
if floor(nOptions/2) < (nOptions/2)
    error(['problem for splitting left/right half-half in every block because number of options/2 is equal to ',num2str(nOptions/2)]);
end
default_LR = repmat([zeros(1,nOptions/2), ones(1,nOptions/2)],1,nBlocks);
for iBlock = 1:nBlocks
    block_trials_idx = (1:nOptions) + nOptions*(iBlock - 1);
    block_rdm_bis = randperm(nOptions);
    default_LR(block_trials_idx) = default_LR( block_trials_idx(block_rdm_bis) );
end % loop on blocks

%% extract information of left/right options + reward or punishment trial
trialOptions.default_LR = default_LR;
trialOptions.R_or_P = cell(1,nTrials);
[trialOptions.R.left,...
    trialOptions.R.right,...
    trialOptions.E.left,...
    trialOptions.E.right] = deal(NaN());
for iTrial = 1:nTrials
    % extract reward or punishment trial
    if R_options(iTrial) < 0
        trialOptions.R_or_P{iTrial} = 'P';
    elseif R_options(iTrial) > 0
        trialOptions.R_or_P{iTrial} = 'R';
    end
    % extract reward/effort per trial (R/E=0 for default option)
    switch default_LR(iTrial)
        case 0
            trialOptions.R.left(iTrial)     = 0;
            trialOptions.R.right(iTrial)    = abs(R_options(iTrial));
            trialOptions.E.left(iTrial)     = 0;
            trialOptions.E.right(iTrial)    = E_options(iTrial);
        case 1
            trialOptions.R.left(iTrial)     = abs(R_options(iTrial));
            trialOptions.R.right(iTrial)    = 0;
            trialOptions.E.left(iTrial)     = E_options(iTrial);
            trialOptions.E.right(iTrial)    = 0;
    end
end

end % function