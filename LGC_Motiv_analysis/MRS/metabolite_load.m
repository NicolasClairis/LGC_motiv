function[metabolites, CRLB] = metabolite_load(subject_id)
% [metabolites, CRLB] = metabolite_load(subject_id)
% metabolite_load will load the metabolite concentrations based on the
% excel file prepared by Arthur Barakat.
%
% INPUTS
% subject_id: list of subject names (if left empty, will load all subjects)
%
% OUTPUTS
% metabolites: structure with all the metabolites
%
% CRLB: structure with Cramér-Rao Lower Bound which gives a lower estimate
% for the variance of an unbiased estimator for each single metabolite (not the
% ratios)

%% working directory
% root = pwd;
list_pcs = {'Lab','Home'};
% which_pc_idx = listdlg('PromptString',{'Lab or home pc?'},...
%     'SelectionMode','single','ListString',list_pcs);
which_pc_idx = 1;
switch list_pcs{which_pc_idx}
    case 'Lab'
        metaboliteFolder = fullfile('M:','human_data_private','analyzed_data','study1');
    case 'Home'
        serverRoot = fullfile(filesep,filesep,'sv-nas1.rcp.epfl.ch',filesep,'Sandi-lab');
        metaboliteFolder = fullfile(serverRoot,'human_data_private','analyzed_data','study1');
end
% cd(metaboliteFolder);
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
        switch all_metabolites{iMet}
            case 'Glu_Gln'
                [metabolites.(ROI_nm).Glx,...
                    CRLB.(ROI_nm).Glx] = deal(NaN(1,NS));
            otherwise
                [metabolites.(ROI_nm).(all_metabolites{iMet}),...
                    CRLB.(ROI_nm).(all_metabolites{iMet})] = deal(NaN(1,NS));
        end
    end % metabolite loop
    [metabolites.(ROI_nm).Gln_div_Glu,...
        metabolites.(ROI_nm).z_antiox,...
        metabolites.(ROI_nm).antiox,...
        metabolites.(ROI_nm).Glu_div_GSH,...
        metabolites.(ROI_nm).Glu_div_Tau,...
        metabolites.(ROI_nm).Glu_div_antiox,...
        metabolites.(ROI_nm).Glu_div_z_antiox,...
        metabolites.(ROI_nm).Glu_div_GABA_div_GSH] = deal(NaN(1,NS));
    
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
        switch met_nm
            case 'Glu_Gln'
                met_nm_bis = 'Glx';
            otherwise
                met_nm_bis = met_nm;
        end
        %% extract data for everybody for the current metabolite
        all_met_data = excelRead_tmp.(met_nm);
        
        %% remove subjects too far (3*SD) from the median
        % compute the median and SD across participants
        median_met_aSubs = median(all_met_data, 1, 'omitnan');
        sd_met_aSubs = std(all_met_data, 0, 1, 'omitnan');
        % identify subjects with metabolite values too far from the SD
        bad_met_subs = (all_met_data > (median_met_aSubs + (sd_met_aSubs.*3))) |...
            (all_met_data < (median_met_aSubs - (sd_met_aSubs.*3)));
        
        %% extract data, removing bad CRLB and far from the mean
        for iS = 1:NS
            sub_nm = subject_id{iS};
            subj_line = strcmp(excelRead_tmp.CID, sub_nm);
            
            %% check Cramer Lower Bound for the subject, ROI and metabolite
            % currently assessed
            subj_CRLB_line = strcmp(CRLB_excelRead_tmp.CID, sub_nm);
            CRLB_tmp = CRLB_excelRead_tmp.(met_nm)(subj_CRLB_line);
            
            %% check if the subject is not an outlier compared to mean and SD of the group
            is_sub_sd_outlier_tmp  = bad_met_subs(subj_line);
            if excelRead_tmp.(met_nm)(subj_line) > 0 &&...
                    CRLB_tmp <= CRLB_threshold &&...
                    is_sub_sd_outlier_tmp == 0
                % ignore subjects where data has not been extracted (=0 when not extracted)
                % and subjects for which the Cramer Lower Bound is too bad
                % (>50%) and subjects which are clearly outliers because
                % their values are too far from the others (>/< mean+3*SD)
                metabolites.(ROI_nm).(met_nm_bis)(iS) = all_met_data(subj_line);
                CRLB.(ROI_nm).(met_nm_bis)(iS) = CRLB_tmp;
            end
        end % subject loop
    end % metabolites loop
    
    %% bonus metabolites: combination of other metabolites
    % perform also division Gln/Glu
    metabolites.(ROI_nm).Gln_div_Glu = metabolites.(ROI_nm).Gln./metabolites.(ROI_nm).Glu;
    
    % extract a pool of antioxidants
    metabolites.(ROI_nm).antiox = metabolites.(ROI_nm).GSH +...
        metabolites.(ROI_nm).Tau;
    metabolites.(ROI_nm).z_antiox = nanzscore(nanzscore(metabolites.(ROI_nm).GSH) +...
        nanzscore(metabolites.(ROI_nm).Tau));
    
    % perform also division Glu/GABA (known as Excitation/Inhibition (E-I) ratio)
    metabolites.(ROI_nm).Glu_div_GABA = metabolites.(ROI_nm).Glu./metabolites.(ROI_nm).GABA;
    % perform (Glu/GABA)/GSH ((E-I) ratio divided by GSH)
    metabolites.(ROI_nm).Glu_div_GABA_div_GSH = metabolites.(ROI_nm).Glu_div_GABA./metabolites.(ROI_nm).GSH;
    
    % Glu/antioxidants
    metabolites.(ROI_nm).Glu_div_GSH = metabolites.(ROI_nm).Glu./metabolites.(ROI_nm).GSH;
    metabolites.(ROI_nm).Glu_div_Tau = metabolites.(ROI_nm).Glu./metabolites.(ROI_nm).Tau;
    metabolites.(ROI_nm).Glu_div_antiox = metabolites.(ROI_nm).Glu./metabolites.(ROI_nm).antiox;
    metabolites.(ROI_nm).Glu_div_z_antiox = metabolites.(ROI_nm).Glu./metabolites.(ROI_nm).z_antiox;
    
    % combination Asp+Lac
    metabolites.(ROI_nm).Asp_plus_Lac = metabolites.(ROI_nm).Asp + metabolites.(ROI_nm).Lac;
    metabolites.(ROI_nm).zAsp_plus_zLac = nanzscore(metabolites.(ROI_nm).Asp) + nanzscore(metabolites.(ROI_nm).Lac);
    
    % Glu U-shape, GSH U-shape, (Glu U-shape)/GSH, Glu/(GSH U-shape), (Glu
    % U-shape)/(GSH U-shape) and (Glu/GSH) U-shape
    metabolites.(ROI_nm).Glu_Ushape = (metabolites.(ROI_nm).Glu - mean(metabolites.(ROI_nm).Glu,'omitnan')).^2;
    metabolites.(ROI_nm).GSH_Ushape = (metabolites.(ROI_nm).GSH - mean(metabolites.(ROI_nm).GSH,'omitnan')).^2;
    metabolites.(ROI_nm).Glu_Ushape_div_GSH = metabolites.(ROI_nm).Glu_Ushape./metabolites.(ROI_nm).GSH;
    metabolites.(ROI_nm).Glu_div_GSH_Ushape = metabolites.(ROI_nm).Glu./metabolites.(ROI_nm).GSH_Ushape;
    metabolites.(ROI_nm).Glu_Ushape_div_GSH_Ushape = metabolites.(ROI_nm).Glu_Ushape./metabolites.(ROI_nm).GSH_Ushape;
    metabolites.(ROI_nm).Glu_div_GSH_ratio_Ushape = (metabolites.(ROI_nm).Glu_div_GSH - mean(metabolites.(ROI_nm).Glu_div_GSH,'omitnan')).^2;
    
    % add Glu/Asp based on Arthur's results with ML on mental effort
    metabolites.(ROI_nm).Glu_div_Asp = metabolites.(ROI_nm).Glu./metabolites.(ROI_nm).Asp;
    % add Lac/GSH, Lac/(GSH+Tau), (Glu+Lac)/(GSH) and
    % (Glu+Lac)/(GSH+Tau) considering Glu+Lac = "toxic" and GSH+Tau prevent
    metabolites.(ROI_nm).Lac_div_GSH = metabolites.(ROI_nm).Lac./metabolites.(ROI_nm).GSH;
    metabolites.(ROI_nm).Lac_div_antiox = metabolites.(ROI_nm).Lac./(metabolites.(ROI_nm).GSH + metabolites.(ROI_nm).Tau);
    metabolites.(ROI_nm).Glu_plus_Lac_div_GSH = (metabolites.(ROI_nm).Glu+metabolites.(ROI_nm).Lac)./metabolites.(ROI_nm).GSH;
    metabolites.(ROI_nm).Glu_plus_Lac_div_antiox = (metabolites.(ROI_nm).Glu + metabolites.(ROI_nm).Lac)./(metabolites.(ROI_nm).GSH + metabolites.(ROI_nm).Tau);
end % ROI loop

% %% go back to root
% cd(root);

end % function