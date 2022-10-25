function[metabolites] = metabolite_load(subject_id)
% [metabolites] = metabolite_load(subject_id)
% metabolite_load will load the metabolite concentrations based on the
% excel file prepared by Arthur Barakat.
%
% INPUTS
% subject_id: list of subject names (if left empty, will load all subjects)
%
% OUTPUTS
% metabolites: structure with all the metabolites

%% working directory
root = pwd;
list_pcs = {'Lab','Home'};
% which_pc_idx = listdlg('PromptString',{'Lab or home pc?'},...
%     'SelectionMode','single','ListString',list_pcs);
which_pc_idx = 1;
switch list_pcs{which_pc_idx}
    case 'Lab'
        metaboliteFolder = fullfile('M:','human_data_private','analyzed_data','study1');
    case 'Home'
        metaboliteFolder = fullfile('L:','human_data_private','analyzed_data','study1');
end
cd(metaboliteFolder);
%% define subject list
if ~exist('subject_id','var') || isempty(subject_id)
    condition = subject_condition;
    [subject_id, NS] = LGCM_subject_selection('study1', condition);
else
    NS = length(subject_id);
end

%% list metabolites
all_metabolites = {'Asp','GABA',...
    'Gln','Glu','GSH','Gly','Ins','Lac','Scyllo','Tau',...
    'NAA','NAAG','PE',...
    'Glu_Gln','GPC_PCho','Cr_PCr'};
% metabolites with bad signal have been removed but you can find the full
% list of theoretically measurable metabolites below:
% all_metabolites = {'Mac','Ala','Asp','PCho','Cr','PCr','GABA',...
%     'Gln','Glu','GSH','Gly','Ins','Lac','NAA','Scyllo','Tau',...
%     'Asc','Glc','NAAG','GPC','PE','Ser',...
%     'NAA_NAAG','Glu_Gln','GPC_PCho','Cr_PCr','Gly_Ins'};
n_metabolites = length(all_metabolites);

%% Cramer-Rao Lower Bound threshold (in percentage)
CRLB_threshold = 0.5;

%% loop dmPFC and aIns
ROIs = {'dmPFC','aIns'};
nROIs = length(ROIs);
for iROI = 1:nROIs
    ROI_nm = ROIs{iROI};
    
    %% prepare metabolites
    for iMet = 1:n_metabolites
        metabolites.(ROI_nm).(all_metabolites{iMet}) = NaN(1,NS);
    end % metabolite loop
    metabolites.(ROI_nm).Gln_div_Glu = NaN(1,NS);
    
    %% load the data
    switch ROI_nm
        case 'dmPFC'
            sheet_nm = 'AdaptedDMPFCMetabolites';
            CRLB_sheet_nm = 'DMPFC_CRLB';
        case 'aIns'
            sheet_nm = 'AdaptedAIMetabolites';
            CRLB_sheet_nm = 'AI_CRLB';
    end
    excelRead_tmp = readtable([metaboliteFolder,filesep,'Metabolites.xlsx'],...
        'Sheet',sheet_nm);
    CRLB_excelRead_tmp = readtable([metaboliteFolder,filesep,'Metabolites.xlsx'],...
        'Sheet',CRLB_sheet_nm);
    
    %% extract the data for each subject
    for iMet = 1:n_metabolites
        met_nm = all_metabolites{iMet};
        % extract data for everybody for the current metabolite
        all_met_data = excelRead_tmp.(met_nm);
        % compute the mean and SD across participants
        mean_met_aSubs = mean(all_met_data, 2, 'omitnan');
        sd_met_aSubs = std(all_met_data, 0, 2, 'omitnan');
        % identify subjects with metabolite values too far from the SD
        bad_met_subs = (all_met_data > mean_met_aSubs + sd_met_aSubs.*3) |...
            (all_met_data < mean_met_aSubs + sd_met_aSubs.*3);
        for iS = 1:NS
            sub_nm = subject_id{iS};
            subj_line = strcmp(excelRead_tmp.CID, sub_nm);
            % check Cramer Lower Bound for the subject, ROI and metabolite
            % currently assessed
            subj_CRLB_line = strcmp(CRLB_excelRead_tmp.CID, sub_nm);
            CRLB_tmp = CRLB_excelRead_tmp.(met_nm)(subj_CRLB_line);
            % check if the subject is not an outlier
            is_sub_sd_outlier_tmp  = bad_met_subs(subj_line);
            if excelRead_tmp.(met_nm)(subj_line) > 0 &&...
                    CRLB_tmp <= CRLB_threshold &&...
                    is_sub_sd_outlier_tmp == 0
                % ignore subjects where data has not been extracted (=0 when not extracted)
                % and subjects for which the Cramer Lower Bound is too bad
                % (>50%) and subjects which are clearly outliers because
                % their values are too far from the others (>/< mean+3*SD)
                metabolites.(ROI_nm).(met_nm)(iS) = all_met_data(subj_line);
            end
        end % subject loop
    end % metabolites
    
    % perform also division Glu/Gln
    metabolites.(ROI_nm).Gln_div_Glu = metabolites.(ROI_nm).Gln./metabolites.(ROI_nm).Glu;
end % ROI loop

%% go back to root
cd(root);

end % function