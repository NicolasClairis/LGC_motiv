function[plasmaM, mb_names, n_mb] = load_plasma_metabolites(subject_id)
% [plasmaM, mb_names, n_mb] = load_plasma_metabolites(subject_id)
% load_plasma_metabolites will load all the metabolite concentrations
% coming from the plasma (amino-acids + metabolic molecules).
%
% INPUTS
% subject_id: list of subjects to include (all subjects by default if left
% empty)
%
% OUTPUTS
% plasmaM: structure with metabolite data with the following subfields:
%   .CID: subject identification number (NS subjects in total)
%   .XXX: concentration of the corresponding metabolic molecule for each subject
%   .filtered.(XXX).CID: list of subjects after filtering mean+3*SD
%   .filtered.(XXX).(XXX): concentration of metabolites after filtering
%   subjects with mean+3*SD
%
% mb_names: cell with list of all the metabolites extracted
%
% n_mb: number of metabolites
%
% Note: all amino-acids are present but Tryptophane (not analyzed by Nestlé).
%
% N. Clairis - january 2024

%% do not allow empty subject_id
if ~exist('subject_id','var') || isempty(subject_id)
    [~, ~, ~, subject_id, NS] = sub_id;
else
    % extract number of subjects
    NS = length(subject_id);
end

%% working directories
plasma_path = 'P:\boulot\postdoc_CarmenSandi\results\plasma';

%% load metabolic data
metabo_excelReadTable = readtable([plasma_path,filesep,'plasma_blood_results.xlsx'],...
    'Sheet','simplified');
metabo_column_names = metabo_excelReadTable.Properties.VariableNames;
nMetaboColumns = length(metabo_column_names);

%% extract the metabolic data
for iMb = 1:nMetaboColumns
    mb_nm = metabo_column_names{iMb};
    
    % initialize variable of interest
    switch mb_nm
        case 'LacticAcid'
            plasmaM.Lac = NaN(1,NS);
        case 'CID'
            plasmaM.CID = cell(1,NS);
        otherwise
            plasmaM.(mb_nm) = NaN(1,NS);
    end
    
    % extract data
    for iS = 1:NS
        sub_nm = subject_id{iS};
        % find subject in the excel file
        plasmaM.CID{iS} = sub_nm;
        sub_idx = find(strcmp(metabo_excelReadTable.CID, ['CID',sub_nm]));
        if ~isempty(sub_idx) && size(sub_idx,2) == 1
            switch mb_nm
                case 'LacticAcid'
                    plasmaM.Lac(iS) = metabo_excelReadTable.LacticAcid(sub_idx);
                case 'CID' % avoid replacing by subject id from the excel file if subjects already entered in inputs
                otherwise
                    plasmaM.(mb_nm)(iS) = metabo_excelReadTable.(mb_nm)(sub_idx);
            end
        else
            error(['Problem with ',mb_nm,' extraction in subject ',sub_nm]);
        end
    end % subject loop
end % loop through metabolites

%% load amino-acids
AA_excelReadTable = readtable([plasma_path,filesep,'Blood_AA.xlsx'],...
    'Sheet','filtered_data_clean');
AA_column_names = AA_excelReadTable.Properties.VariableNames;
nAAColumns = length(AA_column_names);
%% extract amino-acid data and add it to the global structure
for iAA = 1:nAAColumns
    AA_nm = AA_column_names{iAA};
    
    % initialize variable of interest
    switch AA_nm
        case 'L_Alanine'
            plasmaM.Ala = NaN(1,NS);
        case 'L_Arginine'
            plasmaM.Arg = NaN(1,NS);
        case 'L_AsparticA'
            plasmaM.Asp = NaN(1,NS);
        case 'L_Asparagine'
            plasmaM.Asn = NaN(1,NS);
        case 'L_GlutamicA'
            plasmaM.Glu = NaN(1,NS);
        case 'L_Glutamine'
            plasmaM.Gln = NaN(1,NS);
        case 'L_Glycine'
            plasmaM.Gly = NaN(1,NS);
        case 'L_Histidine'
            plasmaM.His = NaN(1,NS);
        case 'L_Isoleucine'
            plasmaM.Ile = NaN(1,NS);
        case 'L_Lysine'
            plasmaM.Lys = NaN(1,NS);
        case 'L_Methionine'
            plasmaM.Met = NaN(1,NS);
        case 'L_Phenylalanine'
            plasmaM.Phe = NaN(1,NS);
        case 'L_Proline'
            plasmaM.Pro = NaN(1,NS);
        case 'L_Serine'
            plasmaM.Ser = NaN(1,NS);
        case 'L_Threonine'
            plasmaM.Thr = NaN(1,NS);
        case 'L_Tyrosine'
            plasmaM.Tyr = NaN(1,NS);
        case 'Valine'
            plasmaM.Val = NaN(1,NS);
        case 'Taurine'
            plasmaM.Tau = NaN(1,NS);
        case 'CID'
        otherwise
            plasmaM.(AA_nm) = NaN(1,NS);
    end
    
    % extract data
    for iS = 1:NS
        sub_nm = subject_id{iS};
        % find subject in the excel file
        sub_idx = find(strcmp(AA_excelReadTable.CID, ['CID',sub_nm]));
        if ~isempty(sub_idx) && size(sub_idx,2) == 1
            switch AA_nm
                case 'L_Alanine'
                    plasmaM.Ala(iS) = AA_excelReadTable.(AA_nm)(sub_idx);
                case 'L_Arginine'
                    plasmaM.Arg(iS) = AA_excelReadTable.(AA_nm)(sub_idx);
                case 'L_AsparticA'
                    plasmaM.Asp(iS) = AA_excelReadTable.(AA_nm)(sub_idx);
                case 'L_Asparagine'
                    plasmaM.Asn(iS) = AA_excelReadTable.(AA_nm)(sub_idx);
                case 'L_GlutamicA'
                    plasmaM.Glu(iS) = AA_excelReadTable.(AA_nm)(sub_idx);
                case 'L_Glutamine'
                    plasmaM.Gln(iS) = AA_excelReadTable.(AA_nm)(sub_idx);
                case 'L_Glycine'
                    plasmaM.Gly(iS) = AA_excelReadTable.(AA_nm)(sub_idx);
                case 'L_Histidine'
                    plasmaM.His(iS) = AA_excelReadTable.(AA_nm)(sub_idx);
                case 'L_Isoleucine'
                    plasmaM.Ile(iS) = AA_excelReadTable.(AA_nm)(sub_idx);
                case 'L_Lysine'
                    plasmaM.Lys(iS) = AA_excelReadTable.(AA_nm)(sub_idx);
                case 'L_Methionine'
                    plasmaM.Met(iS) = AA_excelReadTable.(AA_nm)(sub_idx);
                case 'L_Phenylalanine'
                    plasmaM.Phe(iS) = AA_excelReadTable.(AA_nm)(sub_idx);
                case 'L_Proline'
                    plasmaM.Pro(iS) = AA_excelReadTable.(AA_nm)(sub_idx);
                case 'L_Serine'
                    plasmaM.Ser(iS) = AA_excelReadTable.(AA_nm)(sub_idx);
                case 'L_Threonine'
                    plasmaM.Thr(iS) = AA_excelReadTable.(AA_nm)(sub_idx);
                case 'L_Tyrosine'
                    plasmaM.Tyr(iS) = AA_excelReadTable.(AA_nm)(sub_idx);
                case 'Valine'
                    plasmaM.Val(iS) = AA_excelReadTable.(AA_nm)(sub_idx);
                case 'Taurine'
                    plasmaM.Tau(iS) = AA_excelReadTable.(AA_nm)(sub_idx);
                case 'CID' % avoid replacing by subject id from the excel file if subjects already entered in inputs
                otherwise
                    plasmaM.(AA_nm)(iS) = AA_excelReadTable.(AA_nm)(sub_idx);
            end % amino-acid name
        else
            error(['Problem with ',AA_nm,' extraction in subject ',sub_nm]);
        end
    end % subject loop
end % loop through amino-acids

%% prepare list of metabolites
mb_names = fieldnames(plasmaM);
mb_names(strcmp(mb_names,'CID')) = [];
n_mb = length(mb_names);

%% add filter for outliers
for iM = 1:n_mb
    mb_nm = mb_names{iM};
    [~, ~,...
        cleaned_plasmaM_tmp,...
        idx_goodS] = rmv_outliers_3sd(plasmaM.(mb_nm));
    plasmaM.filtered.(mb_nm).CID = plasmaM.CID(idx_goodS);
    plasmaM.filtered.(mb_nm).(mb_nm) = cleaned_plasmaM_tmp;
end % metabolite loop



end % function