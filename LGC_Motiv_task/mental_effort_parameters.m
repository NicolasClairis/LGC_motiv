function[mentalE_prm] = mental_effort_parameters()
%[mentalE_prm] = mental_effort_parameters()
% mental_effort_parameters definition of the main parameters for the
% mental effort task.
%
% INPUTS
%
% OUTPUTS
% mentalE_prm: structure with main mental effort parameters
%   .sideQuestion: structure with the side corresponding to each answer
%
%   .switchPerc: percentage of switches to implement in a given sequence
%
%   .mental_n_col: structure with colour to use for each number
%
%   .n_max_to_reach_withInstructions: number of subsequent correct answers
%   to reach in the case where instructions are provided
%
%   .n_max_to_reach_withoutInstructions: number of subsequent correct answers
%   to reach in the case where no instructions are provided
%

%% main parameters for mental effort task
% define side of each expected answer
sideQuestion.oE.pair = -1;
sideQuestion.oE.impair = 1;
sideQuestion.hL.low = -1;
sideQuestion.hL.high = 1;

%% define colours to use for numbers font
white = [255 255 255];
col1 = [233 163 201];
% col2 = [161 215 106];
col2 = white;
NbackLastQuestionCol = white;

%% switch percentage = percentage of questions with a switch per
% trial
switchPerc = 1/2;

%% N-back: define how many answers before you need to answer
mentalE_prm.Nback = 2;

%% task switching version: define colours to use for the font of the numbers according to
%  subject number to alternate the type of colour used

% if mod(i_sub,2) == 0
%     mental_n_col.oddEven = col1;
%     mental_n_col.lowHigh = col2;
%     % record mapping with name also (for learning phase)
%     mental_n_col.col1 = 'oddEven';
%     mental_n_col.col2 = 'lowHigh';
% elseif mod(i_sub,2) == 1
%     mental_n_col.oddEven = col2;
%     mental_n_col.lowHigh = col1;
%     % record mapping with name also (for learning phase)
%     mental_n_col.col1 = 'lowHigh';
%     mental_n_col.col2 = 'oddEven';
% end

% NO task switching
mental_n_col.lowHigh = col2;
mental_n_col.oddEven = col1;
% record mapping with name also (for learning phase)
mental_n_col.col1 = 'lowHigh';
mental_n_col.lastQuestion = NbackLastQuestionCol;

%% store all in output
mentalE_prm.sideQuestion    = sideQuestion;
mentalE_prm.switchPerc      = switchPerc;
mentalE_prm.mental_n_col    = mental_n_col;

end % function