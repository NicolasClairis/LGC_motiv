function[mdlType, mdlN] = behavioral_model_selection(mdlType, mdlN)
% [mdlType, mdlN] = behavioral_model_selection(mdlType, mdlN)
% behavioral_model_selection will ask you which model type (bayesian or
% simple) and which model number you want to use for the current analysis
% in case they haven't been initialized already.
%
% INPUTS
% mdlType: which model type do you want to use
% 'bayesian'/'simple' (will be asked if left empty)
%
% mdlN: model number as a string (will be asked if left empty)
%
% OUTPUTS
% mdlType: which model type do you want to use
% 'bayesian'/'simple'
%
% mdlN: model number as a string

%% select model type if not already defined
if ~exist('mdlType','var') || isempty(mdlType)
    listPossibleModels = {'bayesian','simple'};
    mdlType_idx = listdlg('promptstring','Which model type?',...
        'ListString',listPossibleModels);
    mdlType = listPossibleModels{mdlType_idx};
end

%% which model number to use?
if ~exist('mdlN','var') || isempty(mdlN)
    switch mdlType
        case 'bayesian'
            listPossibleModelNumbers = {'3'};
        case 'simple'
            listPossibleModelNumbers = {'1','2','3','4'};
    end
    mdlN_idx = listdlg('promptstring','Which model number?',...
        'ListString',listPossibleModelNumbers);
    mdlN = listPossibleModelNumbers{mdlN_idx};
end
end % function