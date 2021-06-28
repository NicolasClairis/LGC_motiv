function [choice_design] = choice_design_RE_levels(n_RP_levels, n_E_levels, n_trials, RP_condition)
%[choice_design] = choice_design_RE_levels(n_RP_levels, n_E_levels, n_trials, RP_condition)
% choice_design_RE_levels will extract the levels of difficulty for
% each choice trial based on the number of trials and number of reward and
% effort conditions required.
%
% INPUTS
% n_RP_levels: number of reward (or punishment) levels
%
% n_E_levels: number of effort levels
%
% n_trials: total number of trials of the current session for the current
% condition (reward or punishment)
%
% RP_condition: state if reward or punishment condition => adapt the design
% accordingly (the relevant parts of the design matrix are not the same
% depending on the condition)
% 'R': reward condition
% 'P': punishment condition
%
% OUTPUTS
% choice_design: structure with level of reward and effort for left and
% right options
%   .allOptions (trials*RL/RR/EL/ER columns)
%   .reward/punishment
%       .left/.right: reward (or punishment) level for left vs right option
%   .effort
%       .left/right: effort level for left vs right option
%
% See also choice_option_design.m

%% define possible values for R and E
possible_E = 1:n_E_levels;
possible_RP = 1:n_RP_levels;

%% define design matrix size
design_mtrx_potential_options = cell(n_E_levels*2 - 1, n_RP_levels*2 - 1);

%% draw all possible option from the design matrix (based on RL-RR and EL-ER difference)
% (ignore case where RL=RR or where EL=ER because completely
% non-informative (=1D trivial choices based on pure Reward or pure
% Effort difference)
jIncentiveLeft_IncentiveRight = 0;
for iRleft_Rright = -(n_RP_levels - 1):(n_RP_levels - 1) % loop through all possible delta between reward left (RL) and reward right (RR) options
    jIncentiveLeft_IncentiveRight = jIncentiveLeft_IncentiveRight + 1;
    jEleft_Eright = 0;
    for iEL_ER = -(n_E_levels - 1):(n_E_levels - 1) % loop through all possible delta between effort left (EL) and effort right (ER) options
        jEleft_Eright = jEleft_Eright + 1;
        
        % define all possible combinations of left and right option for each box of the matrix
        kOption = 0;
        for iRL_tmp = possible_RP
            iRR_tmp = iRL_tmp - iRleft_Rright; % corresponding R right value for the R left value selected
            if ismember(iRR_tmp, possible_RP)
                for iEL_tmp = possible_E
                    iER_tmp = iEL_tmp - iEL_ER; % corresponding R right value for the R left value selected
                    if ismember(iER_tmp, possible_E)
                        kOption = kOption + 1;
                         % record all possible combinations for this RL-RR
                         % and EL-ER difference
                        design_mtrx_potential_options{jEleft_Eright, jIncentiveLeft_IncentiveRight}(kOption,:) = [iRL_tmp, iRR_tmp, iEL_tmp, iER_tmp];
                    end % if there is a possible E right for the selected E left
                end % E left loop
            end % if there is a possible R right for the selected R left
        end % R left loop
        
    end % EL-ER loop
end % RL-RR loop

%% now define your sample matrix (how many samples you will get from each box of the design matrix)
sample_mtrx = zeros(n_E_levels*2 - 1, n_RP_levels*2 - 1);

% no sample for when RL-RR or EL-ER = 0
sample_mtrx(n_E_levels,:) = 0;
sample_mtrx(:,n_RP_levels) = 0;

% sample uniformly the relevant parts (where RL-RR and EL-ER have same
% sign)
n_samples_per_relevantBox = 2;
% loop through all boxes and put uniform distribution for these ones
switch RP_condition
    case 'R'
        jRleft_Rright = 0;
        for iRleft_Rright = -(n_RP_levels - 1):(n_RP_levels - 1)
            jRleft_Rright = jRleft_Rright + 1;
            jEleft_Eright = 0;
            for iEL_ER = -(n_E_levels - 1):(n_E_levels - 1)
                jEleft_Eright = jEleft_Eright + 1;
                if iRleft_Rright*iEL_ER > 0 % same sign only
                    sample_mtrx(jEleft_Eright, jRleft_Rright) = n_samples_per_relevantBox;
                end % sign RL-RR vs EL-ER
            end % EL-ER loop
        end % RL-RR loop
        
    case 'P'
        jPleft_Pright = 0;
        for iPleft_Pright = -(n_RP_levels - 1):(n_RP_levels - 1)
            jPleft_Pright = jPleft_Pright + 1;
            jEleft_Eright = 0;
            for iEL_ER = -(n_E_levels - 1):(n_E_levels - 1)
                jEleft_Eright = jEleft_Eright + 1;
                if iPleft_Pright*iEL_ER < 0 % opposite signs only
                    sample_mtrx(jEleft_Eright, jPleft_Pright) = n_samples_per_relevantBox;
                end % sign RL-RR vs EL-ER
            end % EL-ER loop
        end % RL-RR loop
end

% take also a few samples from the trivial parts of the matrix (where RL-RR
% and EL-ER have opposite signs = when more rewarded option is also the
% less effortful)
n_samples_matrix_irrelevantParts = 4; % = number of samples for each subpart of the matrix where data is less relevant for the estimation
% check number of boxes per subpart of the matrix
n_irrelevant_boxes = (n_E_levels - 1)*(n_RP_levels - 1);

if n_samples_matrix_irrelevantParts > 0
    switch RP_condition
        case 'R'
            if n_samples_matrix_irrelevantParts == n_irrelevant_boxes
                % if you have a match => extract uniformly the whole space of
                % irrelevant options
                jRleft_Rright = 0;
                for iRleft_Rright = -(n_RP_levels - 1):(n_RP_levels - 1)
                    jRleft_Rright = jRleft_Rright + 1;
                    jEleft_Eright = 0;
                    for iEL_ER = -(n_E_levels - 1):(n_E_levels - 1)
                        jEleft_Eright = jEleft_Eright + 1;
                        if iRleft_Rright*iEL_ER < 0 % opposite sign only
                            sample_mtrx(jEleft_Eright, jRleft_Rright) = 1;
                        end % sign RL-RR vs EL-ER
                    end % EL-ER loop
                end % RL-RR loop
                
            elseif n_samples_matrix_irrelevantParts < n_irrelevant_boxes
                % select randomly which parts of the matrix will be included (to do
                % twice, once one each part of the graph)
                
                % do for RL-RR>0 and EL-ER<0
                vec_irrel1 = randperm(n_irrelevant_boxes, n_samples_matrix_irrelevantParts); % pick 1 box randomly over all the possible boxes
                [E_idx, R_idx] = ind2sub([(n_E_levels - 1), (n_RP_levels - 1)], vec_irrel1); % extract the indexes for the selected box
                for iSample = 1:n_samples_matrix_irrelevantParts
                    sample_mtrx(E_idx(iSample) + n_E_levels, R_idx(iSample)) = 1; % set 1 sample for each box selected (in the correct space)
                end
                
                % same for RL-RR<0 and EL-ER>0
                vec_irrel2 = randperm(n_irrelevant_boxes, n_samples_matrix_irrelevantParts); % pick 1 box randomly over all the possible boxes
                [E_idx, R_idx] = ind2sub([(n_E_levels - 1), (n_RP_levels - 1)], vec_irrel2); % extract the indexes for the selected box
                for iSample = 1:n_samples_matrix_irrelevantParts
                    sample_mtrx(E_idx(iSample), R_idx(iSample) + n_RP_levels) = 1; % set 1 sample for each box selected (in the correct space)
                end
                
            elseif n_samples_matrix_irrelevantParts > n_irrelevant_boxes
                error('case not ready yet: why so many useless trials?');
                % if you can, just distribute samples uniformly, if not possible
                % set uniform distribution as much as possible and then select
                % randomly the other samples
            end
            
        case 'P' % punishment
            if n_samples_matrix_irrelevantParts == n_irrelevant_boxes
                % if you have a match => extract uniformly the whole space of
                % irrelevant options
                jPleft_Pright = 0;
                for iPleft_Pright = -(n_RP_levels - 1):(n_RP_levels - 1)
                    jPleft_Pright = jPleft_Pright + 1;
                    jEleft_Eright = 0;
                    for iEL_ER = -(n_E_levels - 1):(n_E_levels - 1)
                        jEleft_Eright = jEleft_Eright + 1;
                        if iPleft_Pright*iEL_ER > 0 % same sign only
                            sample_mtrx(jEleft_Eright, jPleft_Pright) = 1;
                        end % sign RL-RR vs EL-ER
                    end % EL-ER loop
                end % PL-PR loop
                
            elseif n_samples_matrix_irrelevantParts < n_irrelevant_boxes
                % select randomly which parts of the matrix will be included (to do
                % twice, once one each part of the graph)
                
                % do for PL-PR<0 and EL-ER<0
                vec_irrel1 = randperm(n_irrelevant_boxes, n_samples_matrix_irrelevantParts); % pick 1 box randomly over all the possible boxes
                [E_idx, P_idx] = ind2sub([(n_E_levels - 1), (n_RP_levels - 1)], vec_irrel1); % extract the indexes for the selected box
                for iSample = 1:n_samples_matrix_irrelevantParts
                    sample_mtrx(E_idx(iSample), P_idx(iSample)) = 1; % set 1 sample for each box selected (in the correct space)
                end
                
                % same for PL-PR>0 and EL-ER>0
                vec_irrel2 = randperm(n_irrelevant_boxes, n_samples_matrix_irrelevantParts); % pick 1 box randomly over all the possible boxes
                [E_idx, P_idx] = ind2sub([(n_E_levels - 1), (n_RP_levels - 1)], vec_irrel2); % extract the indexes for the selected box
                for iSample = 1:n_samples_matrix_irrelevantParts
                    sample_mtrx(E_idx(iSample) + n_E_levels, P_idx(iSample) + n_RP_levels) = 1; % set 1 sample for each box selected (in the correct space)
                end
                
            elseif n_samples_matrix_irrelevantParts > n_irrelevant_boxes
                error('case not ready yet: why so many useless trials?');
                % if you can, just distribute samples uniformly, if not possible
                % set uniform distribution as much as possible and then select
                % randomly the other samples
            end
    end
end

%% check that the number of samples matches the total number of trials
n_trials_to_sample = sum(sum(sample_mtrx));
if n_trials_to_sample ~= n_trials
    error(['number of trials = ',num2str(n_trials),' while design matrix only includes ',num2str(n_trials_to_sample)]);
end

%% extract the relevant trials
% -sample randomly in each box if number of samples required is lower than the number of options
% - if the number of samples for a given box is higher than the number of
% options: repeat the same trial
choice_design.allOptions = NaN(n_trials, 4); % RL/RR/EL/ER
jTrial = 0;
for iE = 1:size(design_mtrx_potential_options, 1)
    for iIncentive = 1:size(design_mtrx_potential_options, 2)
        % extract number of samples required for the current box of the
        % matrix
        n_samples_tmp = sample_mtrx(iE, iIncentive);
        
        % see how many combinations match with the corresponding incentiveLeft-incentiveRight and
        % effortLeft-effortRight difference
        n_possible_options_tmp = size(design_mtrx_potential_options{iE, iIncentive}, 1);
        
        %% repeat the same trials several times
        % repeat all combinations equally as much as
        % possible
        n_repeats = floor(n_samples_tmp/n_possible_options_tmp);
        for iRepeat = 1:n_repeats
            
            for iCombination = 1:n_possible_options_tmp
                jTrial = jTrial + 1;
                choice_design.allOptions(jTrial,:) = design_mtrx_potential_options{iE, iIncentive}(iCombination,:);
            end
            
        end % repeat extraction
        
        %% if you can't sample uniformly all possible combinations
        % => peak some combinations randomly
            
        % is there something remaining when you divide number of
        % samples by possible options for the corresponding box?
        residual_samples = mod(n_samples_tmp, n_possible_options_tmp);
        if residual_samples ~= 0
            
            % select randomly which combinations you will extract
            samples_idx = randperm(n_possible_options_tmp, residual_samples);
            samples_to_extract_tmp = design_mtrx_potential_options{iE, iIncentive}(samples_idx,:);
            % store the selected combinations in the output
            for iSample = 1:length(samples_idx)
                jTrial = jTrial + 1;
                choice_design.allOptions(jTrial,:) = samples_to_extract_tmp(iSample,:);
            end
            
        end % if you can't sample uniformly all possible combinations
    end % RL-RR loop
end % EL-ER loop


%% extract the value for each option for the output to be clearer
switch RP_condition
    case 'R'
        choice_design.reward.left   = choice_design.allOptions(:,1);
        choice_design.reward.right  = choice_design.allOptions(:,2);
    case 'P'
        choice_design.punishment.left   = choice_design.allOptions(:,1);
        choice_design.punishment.right  = choice_design.allOptions(:,2);
end
choice_design.effort.left   = choice_design.allOptions(:,3);
choice_design.effort.right  = choice_design.allOptions(:,4);

end % function