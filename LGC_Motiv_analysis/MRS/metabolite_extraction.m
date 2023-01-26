function[metabolite_allSubs, MRS_ROI_nm, metabolite_nm] = metabolite_extraction(study_nm, subject_id)
% [metabolite_allSubs, MRS_ROI_nm, metabolite_nm] = metabolite_extraction(study_nm, subject_id)
% metabolite_extraction will load the metabolites levels for the study and
% the subjects selected and output it in the variable metabolite_allSubs.
%
% INPUTS
% study_nm:
% 'study1': study with dmPFC and anterior insula (aIns)
% 'study2': ventral striatum study
%
% subject_id: list of subjects to consider
%
% OUTPUTS
% metabolite_allSubs: vector with all the levels of metabolites across
% subjects
%
% MRS_ROI_nm: name of the selected ROI for the metabolite extraction
%
% metabolite_nm: name of the metabolite that you extracted

%% extract all metabolites for the subjects selected
[metabolites] = metabolite_load(subject_id);

%% extract 
switch study_nm
    case 'study1'
        %% which ROI?
        ROIs = {'dmPFC','aIns'};
        ROI_idx = listdlg('PromptString','Metabolites in which brain area?',...
            'ListString',ROIs,'SelectionMode','single');
        MRS_ROI_nm = ROIs{ROI_idx};
        %% select metabolite of interest
        metabolite_names = fieldnames(metabolites.(MRS_ROI_nm));
    otherwise
        error(['not ready yet for ',study_nm]);
end
which_metab_idx = listdlg('PromptString','Which metabolite to focus on?',...
    'ListString',metabolite_names,'SelectionMode','single');
metabolite_nm = metabolite_names{which_metab_idx};

%% focus on metabolite and brain area selected
metabolite_allSubs = metabolites.(MRS_ROI_nm).(metabolite_nm);

% in some cases, replace name by generic name:
if strcmp(metabolite_nm,'Glu_Gln')
    metabolite_nm = 'Glx';
end

end % function