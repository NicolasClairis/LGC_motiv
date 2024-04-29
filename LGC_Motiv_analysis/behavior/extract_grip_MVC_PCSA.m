function[MVC, predictedF] = extract_grip_MVC_PCSA(study_nm, subject_id, NS)
% [MVC, predictedF] = extract_grip_MVC_PCSA(study_nm, subject_id, NS)
% extract_grip_MVC_PCSA will extract the maximum voluntary contraction
% (MVC) force and the Physiological Cross-sectional area (PCSA)
%
% INPUTS
% study_nm: study name 'study1' in principle
%
% subject_id: cell with list of subject ids
%
% NS: number of subjects included in subject_id
%
% OUTPUTS
% MVC: maximum voluntary contraction force (in Newtons)
%
% predictedF: theoretical Fmax based on the physiological cros-sectional 
% area (PCSA) measure of the forearm 

%% working directory
list_pcs = {'Lab','Home'};
which_pc_idx = listdlg('PromptString',{'Lab or home pc?'},...
    'SelectionMode','single','ListString',list_pcs);
curr_pc = list_pcs{which_pc_idx};
switch curr_pc
    case 'Lab'
        gitPath = [fullfile('C:','Users','clairis','Desktop','GitHub',...
            'LGC_motiv','LGC_Motiv_results'),filesep];
        dataRoot = ['E:',filesep];
    case 'Home'
        gitPath = [fullfile('C:','Users','Loco','Documents','GitHub',...
            'LGC_motiv','LGC_Motiv_results'),filesep];
        dataRoot = [fullfile('L:','human_data_private','raw_data_subject'),filesep];
end

%% load infos (including info about maximal theoretical force)
excelReadInfosFile = readtable([gitPath,study_nm,filesep,'summary_participants_infos.xlsx'],...
    'Sheet','all_subjects_no_pilots');

%% extract data
[MVC, predictedF] = deal(NaN(1,NS));
for iS = 1:NS
    sub_nm = subject_id{iS};
    
    % extract calibrated Fmax (in Newtons)
    subDataPath = [dataRoot,study_nm,filesep,'CID',sub_nm,filesep,'behavior',filesep];
    MVC_volts = getfield(load([subDataPath,'CID',sub_nm,'_physicalCalib.mat'],'MVC'),'MVC'); % MVC in volts
    MVC(iS) = grip_biopac_volts_to_newtons_conversion(MVC_volts); % convert voltage in Newtons
    
    % extract theoretical Fmax based on arm measurement
    sub_info_idx = strcmp(['CID',sub_nm],excelReadInfosFile.CID);
    pli_a = excelReadInfosFile.PliCutan_Ant_rieur_mm_(sub_info_idx);
    pli_p = excelReadInfosFile.PliCutan_Post_rieur_mm_(sub_info_idx);
    circ = excelReadInfosFile.Circonf_renceDeL_avant_bras_mm_(sub_info_idx);
    armLength = excelReadInfosFile.LongeurDeL_avant_bras_mm_(sub_info_idx);
    predictedF(iS) = Emax_morpho(pli_a, pli_p, circ, armLength);
end % subject loop

end % function