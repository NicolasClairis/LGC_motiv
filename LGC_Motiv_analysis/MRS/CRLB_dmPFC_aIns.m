function[m_CRLB, sem_CRLB, std_CRLB, MRS_ROI_nm, mb_nm, N_goodS] = CRLB_dmPFC_aIns()
% [m_CRLB, sem_CRLB, std_CRLB, MRS_ROI_nm, mb_nm, N_goodS] = CRLB_dmPFC_aIns()
% CRLB_dmPFC_aIns extracts the Cramér-Rao Lower Bound (CRLB) average, sem
% and SD across participants selected for the ROI and metabolite selected
%
% OUTPUTS
% m_CRLB, sem_CRLB, std_CRLB: mean, SEM and SD CRLB for the metabolite
% selected
%
% MRS_ROI_nm: name of the MRS ROI selected
%
% mb_nm: name of the metabolite selected
%
% N_goodS: number of subjects taken into account for the computation (ie
% after removing those with a too high CRLB)

%% subject selection
[study_nm, condition, gender, subject_id, NS] = sub_id;

%% load MRS data
[metabolites, CRLB] = metabolite_load(subject_id);

%% ROI selection
MRS_ROI_names = fieldnames(metabolites);
MRS_ROI_idx = listdlg('PromptString','Which ROI?',...
    'ListString',MRS_ROI_names,'SelectionMode','single');
MRS_ROI_nm = MRS_ROI_names{MRS_ROI_idx};

%% metabolite selection
metabolite_names = fieldnames(metabolites.(MRS_ROI_nm));
MRS_mb_idx = listdlg('PromptString','Which metabolite?',...
    'ListString',metabolite_names,'SelectionMode','single');
mb_nm = metabolite_names{MRS_mb_idx};

%% extract values
[m_CRLB, sem_CRLB, std_CRLB] = mean_sem_sd(CRLB.(MRS_ROI_nm).(mb_nm),2);
N_goodS = sum(~isnan(CRLB.(MRS_ROI_nm).(mb_nm)));

end % function