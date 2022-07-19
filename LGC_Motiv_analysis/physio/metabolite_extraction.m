function[metabolite_allSubs, ROI_nm, metabolite_nm] = metabolite_extraction(study_nm, subject_id)
% [metabolite_allSubs, ROI_nm, metabolite_nm] = metabolite_extraction(study_nm, subject_id)
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
% ROI_nm: name of the selected ROI for the metabolite extraction
%
% metabolite_nm: name of the metabolite that you extracted

%% extract 
switch study_nm
    case 'study1'
        %% which ROI?
        ROIs = {'dmPFC','aIns'};
        nROIs = length(ROIs);
        ROI_idx = spm_input('Metabolites in which brain area?',1,'m',...
            ROIs,1:nROIs,0);
        ROI_nm = ROIs{ROI_idx};
        %% select metabolite of interest
        metabolites = {'Mac','Ala','Asp','PCho','Cr','PCr','GABA',...
            'Gln','Glu','GSH','Gly','Ins','Lac','NAA','Scyllo','Tau',...
            'Asc','Glc','NAAG','GPC','PE','Ser',...
            'NAA_NAAG','Glu_Gln','GPC_PCho','Cr_PCr','Gly_Ins','Gln_div_Glu'};
    otherwise
        error(['not ready yet for ',study_nm]);
end
n_met = length(metabolites);
metabolite_idx = spm_input('Which metabolite to focus on?',1,'m',...
    metabolites,1:n_met,0);
metabolite_nm = metabolites{metabolite_idx};

%% extract all metabolites for the subjects selected
[metabolites] = metabolite_load(subject_id);
%% focus on metabolite and brain area selected
metabolite_allSubs = metabolites.(ROI_nm).(metabolite_nm);

% in some cases, replace name by generic name:
if strcmp(metabolite_nm,'Glu_Gln')
    metabolite_nm = 'Glx';
end

end % function