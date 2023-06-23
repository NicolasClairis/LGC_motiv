function [side, npair, lottery] = designLearning75_bis(ntrialsGain, ntrialsLoss, ntrialsNeutral, indpProb)
% design matrix for taskLearning75
%
% Inputs:
% ntrialsGain: number of gain/neutral pair trials
% ntrialsLoss: number of loss/neutral pair trials
% ntrialsNeutral: number of neutral/neutral pair trials
% indpProb: independent probabilities for two items of a pair (1) or dependent probabilities (0)?
% In the latter, if one is winning, the other is not and vice-versa
% (winning pair, same logic for loss pair). In the former, you could have
% two items giving a winning (or neutral or losing) outcome in the same
% trial.
%
% Outputs:
% side: side for the "best" cue, half on the left, half on the right across
% trials
% npair: different types of pairs (gain/neutral/loss) ordered across trials
% lottery: probabilities associated to each cue for each trial: 1 = most
% likely outcome (big Gain/small gain for corresponding item for gains,
% neutral/neutral for loss and big loss/small loss for corresponding item),
% -1: less likely outcome = small gain for item generally associated with
% big gain, neutral/neutral and small loss for item of big loss
%
% nunlikely: probability for each item of the pair (here 75% => should be a total number compatible
% with this value)
%
% Script made by Nicolas Clairis on february 2017 for fMRI task
% Initially based on designLearning75.m by Nicolas Borderies
%
% See also taskLearning75.m

% total nber of trials
totalTrial = ntrialsGain + ntrialsNeutral + ntrialsLoss;

%% divide trials into smaller sequences (to spread probabilities more equally)
nSmallSeq = 6; % number of smaller sequences in which you want to divide the bigger one to better spread the probabilities
nsmallGainSeq = ntrialsGain/nSmallSeq;
nsmallLossSeq = ntrialsLoss/nSmallSeq;
% check whether mini-block parameters are ok
if nsmallGainSeq ~= floor(nsmallGainSeq) || nsmallLossSeq ~= floor(nsmallLossSeq)
    warning('problem in nSmallSeq: you have to change mini-blocks parameters or number of gain/loss trials please');
    return;
end

%% select randmoly order of 25% trials where winning item loses
% and also 25% trials where neutral item wins
% (do the same for losses)
% apply this to small blocks and then add it to global list of items
lottery_positivePairs = NaN(2,ntrialsGain); % First line for winning item (75%), second line for less winning item (25%)
lottery_lossPairs = NaN(2,ntrialsLoss); % First line for losing item (75%), second line for less losing item (25%)
lottery_neutralPairs = ones(2,ntrialsNeutral); % all ones since probabilities are always of 100% to get a neutral outcome for neutral pairs
nunlikely = nsmallGainSeq/4; % 1/4(=25%) of unlikely outcomes per item
for iBlock = 1:nSmallSeq
    % win pair
    winSeq_idx = (1:nsmallGainSeq) + (iBlock - 1)*nsmallGainSeq;
    % 75% winning item
    winningItemSmallSeqTrials = ones(1,nsmallGainSeq); % 1 = most likely outcome (75%)
    nonGain25Trials = randperm(nsmallGainSeq,nunlikely); % 25% of non-winning trials for the winning item
    winningItemSmallSeqTrials(nonGain25Trials) = -1; % these are the trials with the less likely outcome (-1 = 25%)
    lottery_positivePairs(1,winSeq_idx) = winningItemSmallSeqTrials;
    % 25% winning item
    nonwinningItemSmallSeqTrials = ones(1,nsmallGainSeq); % 1 = most likely outcome (75%)
    gain25Trials = randperm(nsmallGainSeq,nunlikely); % 25% of winning trials for the non-winning item
    nonwinningItemSmallSeqTrials(gain25Trials) = -1; % these are the trials with the less likely outcome (-1 = 25%)
    lottery_positivePairs(2,winSeq_idx) = nonwinningItemSmallSeqTrials;
    
    % lose pair
    loseSeq_idx = (1:nsmallLossSeq) + (iBlock - 1)*nsmallLossSeq;
    % 75% lose item
    losingItemSmallSeqTrials = ones(1,nsmallLossSeq); % 1 = most likely outcome (75%)
    nonGain25Trials = randperm(nsmallLossSeq,nunlikely); % 25% of non-losing trials for the losing item
    losingItemSmallSeqTrials(nonGain25Trials) = -1; % these are the trials with the less likely outcome (-1 = 25%)
    lottery_lossPairs(1,loseSeq_idx) = losingItemSmallSeqTrials;
    % 25% lose item
    nonlosingItemSmallSeqTrials = ones(1,nsmallLossSeq); % 1 = most likely outcome (75%)
    gain25Trials = randperm(nsmallLossSeq,nunlikely); % 25% of losing trials for the non-losing item
    nonlosingItemSmallSeqTrials(gain25Trials) = -1; % these are the trials with the less likely outcome (-1 = 25%)
    lottery_lossPairs(2,loseSeq_idx) = nonlosingItemSmallSeqTrials;
end

%% randomize side of cues (50/50) for each pair (side)
% win pair
winSide = [ repmat(-1,1,ntrialsGain/2) ones(1,ntrialsGain/2)]; % 50/50 for sides
winSidePerm = randperm(ntrialsGain); % randomize order for each side
winSide = winSide(winSidePerm);
% neutral pair
neutralSide = [ repmat(-1,1,ntrialsNeutral/2) ones(1,ntrialsNeutral/2)]; % 50/50 for sides
neutralSidePerm = randperm(ntrialsNeutral);% randomize order for each side
neutralSide = neutralSide(neutralSidePerm);
% loss pair
lossSide = [ repmat(-1,1,ntrialsLoss/2) ones(1,ntrialsLoss/2)]; % 50/50 for sides
lossSidePerm = randperm(ntrialsLoss);% randomize order for each side
lossSide = lossSide(lossSidePerm);

%% randomize type of pair (gain/neutral/loss) (npair) in mini-blocks 
miniBlockSize = 5; % 2 R/r trials, 2 L/l trials and 1 neutral/neutral trial
pairType = [ones(1,2), 2*ones(1,1), 3*ones(1,2)]; % one id for each type of pair (gain=1, neutral=2, loss=3)
nPairMiniBlocks = totalTrial/miniBlockSize;
npair = NaN(1,totalTrial);
if (ntrialsGain == ntrialsLoss) && (ntrialsGain == 2*ntrialsNeutral)
    for iPairBlock = 1:nPairMiniBlocks
        randomizePairs = randperm(miniBlockSize); % randomize order of pairs in the PairMiniBlock sequence
        npair_idx = (1 + miniBlockSize*(iPairBlock - 1)):(miniBlockSize*iPairBlock);
        npair(npair_idx) = pairType(randomizePairs); % extract miniblock in npair
    end
else
    warning(['Script not ready for case where you don''t have the following attributes:',...
        ' number of gains = number of loss trials',...
        'AND number of neutral trials = half of the gain and loss trials']);
    return;
end

%% create lottery and side variable with all pairs in the correct order now
% find trials index for each condition
gainTrials = find(npair == 1);
neutralTrials = find(npair == 2);
lossTrials = find(npair == 3);
% inject each condition side variable into global side variable
side = NaN(1,totalTrial);
side(gainTrials) = winSide;
side(neutralTrials) = neutralSide;
side(lossTrials) = lossSide;
% inject each condition lottery variable into global lottery variable
if indpProb == 1 % independent variables => first line is for the best item of a pair while the second line is for the worst
    lottery = NaN(2,totalTrial);
    lottery(1:2,gainTrials) = lottery_positivePairs;
    lottery(1:2,neutralTrials) = lottery_neutralPairs;
    lottery(1:2,lossTrials) = lottery_lossPairs;
elseif indpProb == 0 % dependent variables => only consider the best item of a pair
    lottery = NaN(1,totalTrial);
    lottery(1,gainTrials) = lottery_positivePairs(1,:);
    lottery(1,neutralTrials) = lottery_neutralPairs(1,:);
    lottery(1,lossTrials) = lottery_lossPairs(1,:);
end
end