function [choice_design] = LGCM_choice_design_RE_levels(n_R_levels, n_E_levels, n_trials)
%[choice_design] = LGCM_choice_design_RE_levels(n_R_levels, n_E_levels, n_trials)
%
% INPUTS
% n_R_levels: number of reward levels
%
% n_E_levels: number of effort levels
%
% n_trials: total number of trial of the current session
%
% OUTPUTS
% choice_design: structure with level of reward and effort for left and
% right options
%   .allOptions (trials*RL/RR/EL/ER columns)
%   .reward
%       .left/.right: reward level for left vs right option
%   .effort
%       .left/right: effort level for left vs right option
%
% See also LGCM_choice_option_design.m

%% define possible values for R and E
possible_E = 1:n_E_levels;
possible_R = 1:n_R_levels;

%% define design matrix size
design_mtrx_potential_options = cell(n_E_levels*2 - 1, n_R_levels*2 - 1);

%% draw all possible option from the design matrix (based on RL-RR and EL-ER difference)
% (ignore case where RL=RR or where EL=ER because completely
% non-informative (=1D trivial choices based on pure Reward or pure
% Effort difference)
jRL_RR = 0;
for iRL_RR = -(n_R_levels - 1):(n_R_levels - 1) % loop through all possible delta between reward left (RL) and reward right (RR) options
    jRL_RR = jRL_RR + 1;
    jEL_ER = 0;
    for iEL_ER = -(n_E_levels - 1):(n_E_levels - 1) % loop through all possible delta between effort left (EL) and effort right (ER) options
        jEL_ER = jEL_ER + 1;
        
        % define all possible combinations of left and right option for each box of the matrix
        kOption = 0;
        for iRL_tmp = possible_R
            iRR_tmp = iRL_tmp - iRL_RR; % corresponding R right value for the R left value selected
            if ismember(iRR_tmp, possible_R)
                for iEL_tmp = possible_E
                    iER_tmp = iEL_tmp - iEL_ER; % corresponding R right value for the R left value selected
                    if ismember(iER_tmp, possible_E)
                        kOption = kOption + 1;
                         % record all possible combinations for this RL-RR
                         % and EL-ER difference
                        design_mtrx_potential_options{jEL_ER, jRL_RR}(kOption,:) = [iRL_tmp, iRR_tmp, iEL_tmp, iER_tmp];
                    end % if there is a possible E right for the selected E left
                end % E left loop
            end % if there is a possible R right for the selected R left
        end % R left loop
        
    end % EL-ER loop
end % RL-RR loop

%% now define your sample matrix (how many samples you will get from each box of the design matrix)
sample_mtrx = zeros(n_E_levels*2 - 1, n_R_levels*2 - 1);

% no sample for when RL-RR or EL-ER = 0
sample_mtrx(n_E_levels,:) = 0;
sample_mtrx(:,n_R_levels) = 0;

% sample uniformly the relevant parts (where RL-RR and EL-ER have same
% sign)
n_samples_per_relevantBox = 2;
% loop through all boxes and put uniform distribution for these ones
jRL_RR = 0;
for iRL_RR = -(n_R_levels - 1):(n_R_levels - 1)
    jRL_RR = jRL_RR + 1;
    jEL_ER = 0;
    for iEL_ER = -(n_E_levels - 1):(n_E_levels - 1)
        jEL_ER = jEL_ER + 1;
        if iRL_RR*iEL_ER > 0 % same sign only
            sample_mtrx(jEL_ER, jRL_RR) = n_samples_per_relevantBox;
        end % sign RL-RR vs EL-ER
    end % EL-ER loop
end % RL-RR loop

% take also a few samples from the trivial parts of the matrix (where RL-RR
% and EL-ER have opposite signs = when more rewarded option is also the
% less effortful)
n_samples_matrix_irrelevantParts = 4; % = number of samples for each subpart of the matrix where data is less relevant for the estimation
% check number of boxes per subpart of the matrix
n_irrelevant_boxes = (n_E_levels - 1)*(n_R_levels - 1);

if n_samples_matrix_irrelevantParts == n_irrelevant_boxes
    % if you have a match => extract uniformly the whole space
    jRL_RR = 0;
    for iRL_RR = -(n_R_levels - 1):(n_R_levels - 1)
        jRL_RR = jRL_RR + 1;
        jEL_ER = 0;
        for iEL_ER = -(n_E_levels - 1):(n_E_levels - 1)
            jEL_ER = jEL_ER + 1;
            if iRL_RR*iEL_ER < 0 % same sign only
                sample_mtrx(jEL_ER, jRL_RR) = 1;
            end % sign RL-RR vs EL-ER
        end % EL-ER loop
    end % RL-RR loop
    
elseif n_samples_matrix_irrelevantParts < n_irrelevant_boxes
    % select randomly which parts of the matrix will be included (to do
    % twice, once one each part of the graph)
    
    % do for RL-RR>0 and EL-ER<0
    vec_irrel1 = randperm(n_irrelevant_boxes, n_samples_matrix_irrelevantParts); % pick 1 box randomly over all the possible boxes
    [E_idx, R_idx] = ind2sub([(n_E_levels - 1), (n_R_levels - 1)], vec_irrel1); % extract the indexes for the selected box
    for iSample = 1:n_samples_matrix_irrelevantParts
        sample_mtrx(E_idx(iSample) + n_E_levels, R_idx(iSample)) = 1; % set 1 sample for each box selected (in the correct space)
    end
    
    % same for RL-RR<0 and EL-ER>0
    vec_irrel2 = randperm(n_irrelevant_boxes, n_samples_matrix_irrelevantParts); % pick 1 box randomly over all the possible boxes
    [E_idx, R_idx] = ind2sub([(n_E_levels - 1), (n_R_levels - 1)], vec_irrel2); % extract the indexes for the selected box
    for iSample = 1:n_samples_matrix_irrelevantParts
        sample_mtrx(E_idx(iSample), R_idx(iSample) + n_R_levels) = 1; % set 1 sample for each box selected (in the correct space)
    end
    
    elseif n_samples_matrix_irrelevantParts > n_irrelevant_boxes
        error('case not ready yet: why so many useless trials?');
        % if you can, just distribute samples uniformly, if not possible
        % set uniform distribution as much as possible and then select
        % randomly the other samples
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
    for iR = 1:size(design_mtrx_potential_options, 2)
        % extract number of samples required for the current box of the
        % matrix
        n_samples_tmp = sample_mtrx(iE, iR);
        
        % see how many combinations match with the corresponding RL-RR and
        % EL-ER difference
        n_possible_options_tmp = size(design_mtrx_potential_options{iE, iR}, 1);
        
        %% repeat the same trials several times
        % repeat all combinations equally as much as
        % possible
        n_repeats = floor(n_samples_tmp/n_possible_options_tmp);
        for iRepeat = 1:n_repeats
            
            for iCombination = 1:n_possible_options_tmp
                jTrial = jTrial + 1;
                choice_design.allOptions(jTrial,:) = design_mtrx_potential_options{iE, iR}(iCombination,:);
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
            samples_to_extract_tmp = design_mtrx_potential_options{iE, iR}(samples_idx,:);
            % store the selected combinations in the output
            for iSample = 1:length(samples_idx)
                jTrial = jTrial + 1;
                choice_design.allOptions(jTrial,:) = samples_to_extract_tmp(iSample,:);
            end
            
        end % if you can't sample uniformly all possible combinations
    end % RL-RR loop
end % EL-ER loop


%% extract the value for each option for the output to be clearer
choice_design.reward.left   = choice_design.allOptions(:,1);
choice_design.reward.right  = choice_design.allOptions(:,2);
choice_design.effort.left   = choice_design.allOptions(:,3);
choice_design.effort.right  = choice_design.allOptions(:,4);

end % function