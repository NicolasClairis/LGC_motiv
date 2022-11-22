function[blood, sub_List] = load_blood_NAD(study_nm)
% [blood, sub_List] = load_blood_NAD(study_nm)
% load_blood_NAD will load the blood NAD values for the study defined in
% input 'study_nm'.
%
% INPUTS
% study_nm: string with study name 'study1'/'study2'
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
    'LGC_motiv','LGC_Motiv_results',study_nm);
%% extract all blood data
switch study_nm
    case 'study1'
        bloodFilePath = fullfile(dataResultsPath,...
            'blood_NAD','20221114_CS_NAD_whole_blood_data.xlsx');
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
    blood.(blood_m_nm) = NAD_blood_table.(blood_m_nm);
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

%% extract list of subjects
sub_List = NAD_blood_table.ID;

end % function