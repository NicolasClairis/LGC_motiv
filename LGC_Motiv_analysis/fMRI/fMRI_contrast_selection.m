function[selectedContrast, selectedCon_nm] = fMRI_contrast_selection(con_names)
% [selectedContrast, selectedCon_nm] = fMRI_contrast_selection(con_names)
% fMRI_contrast_selection asks you to select one contrast of interest.
%
% INPUTS
% con_names: list of contrast names
%
% OUTPUTS
% selectedContrast: index of the selected contrast
%
% selectedCon_nm: full name of the selected contrast

%% select the contrast of interest
selectedContrast = listdlg('PromptString','Which contrast?',...
    'SelectionMode','single',...
    'ListString',con_names);
selectedCon_nm = con_names{selectedContrast};

end % function