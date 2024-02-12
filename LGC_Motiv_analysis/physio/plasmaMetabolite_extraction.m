function[metabolite_allSubs, metabolite_nm] = plasmaMetabolite_extraction(subject_id)
% [metabolite_allSubs, metabolite_nm] = plasmaMetabolite_extraction(subject_id)
% plasmaMetabolite_extraction will load the plasma metabolites levels for 
% the subjects and metabolite selected and it will output it in the 
% variable metabolite_allSubs along with the name of the  selected metabolite.
%
% INPUTS
% study_nm:
% 'study1': study with dmPFC and anterior insula (aIns)
% 'study2': ventral striatum study
%
% subject_id: list of subjects to consider
%
% OUTPUTS
% metabolite_allSubs: vector with all the levels of plasma metabolites 
% across subjects
%
% metabolite_nm: name of the metabolite that you extracted

%% extract all metabolites for the subjects selected
[plasmaMetabolites, mb_names] = load_plasma_metabolites(subject_id);

%% select metabolite of interest
which_metab_idx = listdlg('PromptString','Which metabolite to focus on?',...
    'ListString',mb_names,'SelectionMode','single');
metabolite_nm = mb_names{which_metab_idx};

%% focus on metabolite and brain area selected
metabolite_allSubs = plasmaMetabolites.(metabolite_nm);

end % function