function[NMP] = extract_Nback_NMP(study_nm, subject_id, NS)
% [NMP] = extract_Nback_NMP(study_nm, subject_id, NS)
% extract_Nback_NMP will extract the number for maximal performance (NMP)
%
% INPUTS
% study_nm: study name 'study1' in principle
%
% subject_id: cell with list of subject ids
%
% NS: number of subjects included in subject_id
%
% OUTPUTS
% NMP: number for maximal performance (NMP) in the 2-back task

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

%% extract data
NMP = NaN(1,NS);
for iS = 1:NS
    sub_nm = subject_id{iS};
    
    % extract calibrated Fmax (in Newtons)
    subDataPath = [dataRoot,study_nm,filesep,'CID',sub_nm,filesep,'behavior',filesep];
    NMP(iS) = getfield(load([subDataPath,'CID',sub_nm,'_mentalCalib.mat'],'NMP'),'NMP');
end % subject loop

end % function