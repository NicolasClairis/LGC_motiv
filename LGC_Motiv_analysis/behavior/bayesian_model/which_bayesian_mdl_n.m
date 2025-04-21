function[mdl_n, mdl_n_nm] = which_bayesian_mdl_n(mdl_n)
% [mdl_n, mdl_n_nm] = which_bayesian_mdl_n(mdl_n)
% which_bayesian_mdl_n will ask which computational model number you want to select
% if the input is empty.
%
% INPUTS
% mdl_n: model number (if empty or not matching one of the possible models,
% the script will ask the number, otherwise, it will just output the same
% variable)
%
% OUTPUTS
% mdl_n: model number selected (as a numeric variable)
%
% mdl_n_nm: model number as a string
%
% See also computational_mdl and computational_mdl_prm
%
% Function written by N.Clairis - may 2024


list_potential_models = 1:7;
list_potential_model_names = string(list_potential_models);
if ~exist('mdl_n','var') || isempty(mdl_n) || ~ismember(mdl_n, list_potential_models)
    mdl_n_name_idx = listdlg('PromptString','What is the model number?',...
        'ListString',list_potential_model_names,...
        'SelectionMode','single');
    mdl_n_nm = list_potential_model_names{mdl_n_name_idx};
    mdl_n = str2double(mdl_n_nm);
else % if mdl_n already defined, just transform the model number in a string for the output
    mdl_n_nm = num2str(mdl_n);
end

end % function