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
sideQuestion.hL.low = -1;
sideQuestion.hL.high = 1;

%% define colours to use for numbers font
white = [255 255 255];

%% N-back: define how many answers before you need to answer
mentalE_prm.Nback = 2;

%% task switching version: define colours to use for the font of the numbers according to
%  subject number to alternate the type of colour used

% NO task switching
mental_n_col.lowHigh = white;
% record mapping with name also (for learning phase)
mental_n_col.col1 = 'lowHigh';
mental_n_col.lastQuestion = white;

%% store all in output
mentalE_prm.sideQuestion    = sideQuestion;
mentalE_prm.mental_n_col    = mental_n_col;

end % function