function[choiceOptions] = design_choiceOptions(n_R_levels_withDefault, n_E_levels_withDefault, punishment_yn, nTrials)
% [choiceOptions] = design_choiceOptions(n_R_levels_withDefault, n_E_levels_withDefault, punishment_yn, nTrials)
% design_choiceOptions will create a potential design matrix for the version
% with a default option.
% 
% INPUTS
% n_R_levels_withIP: number of money levels (including the default option)
%
% n_E_levels_withIP: number of effort levels (including the default option)
%
% punishment_yn:
% 'yes': add punishments
% 'no': only rewards
%
% nTrials: number of trials in total
%
% OUTPUTS
% choiceOptions: structure with information for each choice trial

%% remove IP
n_R_levels_withoutDefault = n_R_levels_withDefault - 1;
n_E_levels_withoutDefault = n_E_levels_withDefault - 1;

%% define potential options for 1 mini-block
switch punishment_yn
    case 'no'
        options = 1:(n_R_levels_withoutDefault*n_E_levels_withoutDefault);
        [R_optionsPerBlock,...
            E_optionsPerBlock] = deal(NaN(1,(n_R_levels_withoutDefault*n_E_levels_withoutDefault)));
        jOption = 0;
        for iRoption = 1:n_R_levels_withoutDefault
            for iEoption = 1:n_E_levels_withoutDefault
                jOption = jOption + 1;
                R_optionsPerBlock(jOption) = iRoption;
                E_optionsPerBlock(jOption) = iEoption;
            end % effort loop
        end % reward loop
    case 'yes'
        options = 1:((n_R_levels_withoutDefault*n_E_levels_withoutDefault)*2);
        [R_optionsPerBlock,...
            E_optionsPerBlock] = deal(NaN(1,(n_R_levels_withoutDefault*n_E_levels_withoutDefault)*2));
        jOption = 0;
        for iRoption = [(-n_R_levels_withoutDefault):(-1), 1:n_R_levels_withoutDefault]
            for iEoption = 1:n_E_levels_withoutDefault
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
matrxOk = 0;
nRP_repeatThreshold = 4; % script will not accept a matrix where consecutive reward or punishment trials are more than this number
while matrxOk == 0
    % randomize the blocks
    for iBlock = 1:nBlocks
        block_trials_idx = (1:nOptions) + nOptions*(iBlock - 1);
        block_rdm = randperm(nOptions);
        R_options(block_trials_idx) = R_options( block_trials_idx(block_rdm) );
        E_options(block_trials_idx) = E_options( block_trials_idx(block_rdm) );
    end % loop on blocks

    % verify that the reward and punishment options are alternating well
    % => if not perform another randomization
    jR = 0; % index for consecutive reward trials
    jP = 0; % index for consecutive punishment trials
    [maxConsecutiveR, maxConsecutiveP] = deal(0);
    for iTrial = 1:nTrials
        if R_options(iTrial) < 0 % punishment trial
            jR = 0;
            jP = jP + 1;
            maxConsecutiveP = max(maxConsecutiveP, jP);
        elseif R_options(iTrial) > 0 % reward trial
            jR = jR + 1;
            jP = 0;
            maxConsecutiveR = max(maxConsecutiveR, jR);
        end
    end % loop through trials
    
    if (maxConsecutiveR <= nRP_repeatThreshold) && (maxConsecutiveP <= nRP_repeatThreshold)
        matrxOk = 1;
    end
end

%% randomize the left/right side of the default option within each block
if floor(nOptions/2) < (nOptions/2)
    error(['problem for splitting left/right half-half in every block because number of options/2 is equal to ',num2str(nOptions/2)]);
end
default_LR = repmat([-ones(1,nOptions/2), ones(1,nOptions/2)],1,nBlocks);
for iBlock = 1:nBlocks
    block_trials_idx = (1:nOptions) + nOptions*(iBlock - 1);
    block_rdm_bis = randperm(nOptions);
    default_LR(block_trials_idx) = default_LR( block_trials_idx(block_rdm_bis) );
end % loop on blocks

%% extract information of left/right options + reward or punishment trial
choiceOptions.default_LR = default_LR;
choiceOptions.R_or_P = cell(1,nTrials);
[choiceOptions.R.left,...
    choiceOptions.R.right,...
    choiceOptions.E.left,...
    choiceOptions.E.right] = deal(NaN());
for iTrial = 1:nTrials
    % extract reward or punishment trial
    if R_options(iTrial) < 0
        choiceOptions.R_or_P{iTrial} = 'P';
    elseif R_options(iTrial) > 0
        choiceOptions.R_or_P{iTrial} = 'R';
    end
    % extract reward/effort per trial (R/E=0 for default option)
    switch default_LR(iTrial)
        case -1 % default on the left
            choiceOptions.R.left(iTrial)     = 0;
            choiceOptions.R.right(iTrial)    = abs(R_options(iTrial));
            choiceOptions.E.left(iTrial)     = 0;
            choiceOptions.E.right(iTrial)    = E_options(iTrial);
        case 1 % default on the right
            choiceOptions.R.left(iTrial)     = abs(R_options(iTrial));
            choiceOptions.R.right(iTrial)    = 0;
            choiceOptions.E.left(iTrial)     = E_options(iTrial);
            choiceOptions.E.right(iTrial)    = 0;
    end
end

end % function