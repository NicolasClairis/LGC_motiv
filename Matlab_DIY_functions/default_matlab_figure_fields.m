% Function just to remind how to obtain all the default fields of matlab
% for figures to facilitate their identification and modification

%% extract all fields
get(groot,'factory')
%% if you want to change any of these fields as default, just replace 'factory'
% by 'default' in the name and call:
% set(0,'defaultXXX', defaultValue);
%
% For example:
% set(0,'defaultLineLineWidth',1); will set the default line width to 1