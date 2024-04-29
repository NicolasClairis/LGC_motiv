function[blood, sub_List] = load_blood_NAD(study_nm, subject_id)
% [blood, sub_List] = load_blood_NAD(study_nm, subject_id)
% load_blood_NAD will load the blood NAD values for the study defined in
% input 'study_nm'.
%
% INPUTS
% study_nm: string with study name 'study1'/'study2'
%
% subject_id: list of subjects to extract
%
% OUTPUTS
% blood: structure with blood info for all subjects
%
% sub_List: list of subjects corresponding to each raw in the blood
% variables

%% working directories
which_pc = 'Lab';
switch which_pc
    case 'Lab'
        gitPath = fullfile('C:','Users','clairis','Desktop');
    case 'home'
        gitPath = fullfile('C:','Users','Loco','Documents');
end 
dataResultsPath = fullfile(gitPath, 'GitHub',...
    'LGC_motiv','LGC_Motiv_results',study_nm,'blood','blood_NAD');
%% extract all blood data
switch study_nm
    case 'study1'
        bloodFilePath = fullfile(dataResultsPath,...
            '20221114_CS_NAD_whole_blood_data.xlsx');
    otherwise
        error(['study ',study_nm,' not ready']);
end
NAD_blood_table = readtable(bloodFilePath,...
    'Sheet','clean');
blood_metabolites = {'Nam','NMN','NR',...
    'NAD','NADH','NADP','NADPH',...
    'MeNam','MeXPY'};
nNADHmetab = length(blood_metabolites);

for iBloodM = 1:nNADHmetab
    blood_m_nm = blood_metabolites{iBloodM};
    blood.(blood_m_nm) = NAD_blood_table.(blood_m_nm)'; % re-shape vector from N*1 to 1*N for my further scripts which are more adapted with this shape
end % metabolite loop

%% add ratios
blood.NAD_div_NADH = blood.NAD./blood.NADH;
blood.NADP_div_NADPH = blood.NADP./blood.NADPH;

%% add pools of different metabolites
blood.total_NAD_precursors = blood.Nam + blood.NMN + blood.NR;
blood.total_NAD = blood.NAD + blood.NADH;
blood.total_NAD_with_precursors = blood.Nam + blood.NMN + blood.NR +...
    blood.NAD + blood.NADH +...
    blood.NADP + blood.NADPH;
blood.total_NAD_with_byproducts = blood.Nam + blood.NMN + blood.NR +...
    blood.NAD + blood.NADH +...
    blood.NADP + blood.NADPH +...
    blood.MeNam + blood.MeXPY;
blood.total_NAD_byproducts = blood.MeNam + blood.MeXPY;

%% extract list of subjects and blood-related data only for those subjects 
% entered in input if a list was provided
if ~exist('subject_id','var') || isempty(subject_id)
    sub_List = NAD_blood_table.ID;
else % filter subjects entered in input
    blood_metabolites_bis = fieldnames(blood);
    % extract list of subjects
    NS = length(subject_id);
    sub_List = cell(NS,1);
    for iS = 1:NS
        sub_nm = subject_id{iS};
        sub_List{iS} = ['CID',sub_nm];
    end
    % extract blood data
    for iB = 1:length(blood_metabolites_bis)
        blood_nm = blood_metabolites_bis{iB};
        
        % extract the data for the subjects entered in input
        blood_tmp = NaN(1,NS);
        for iS = 1:NS
            sub_nm = subject_id{iS};
            sub_idx = strcmp(NAD_blood_table.ID, ['CID',sub_nm]);
            if ~isempty(sub_idx) && size(sub_idx,2) == 1
                blood_tmp(iS) = blood.(blood_nm)(sub_idx);
            else
                error(['Problem with ',blood_nm,' extraction in subject ',sub_nm]);
            end
        end % subject loop
        
        % replace output with list filtered with the correct number of
        % subjects
        blood.(blood_nm) = blood_tmp;
    end % loop over blood metabolites

end % subject filter
end % function