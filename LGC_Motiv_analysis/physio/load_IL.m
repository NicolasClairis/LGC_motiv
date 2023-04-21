function[IL_data] = load_IL(study_nm)
% [IL_data] = load_IL(study_nm)
% load_IL will load the interleukins data for the study defined in input.
%
% INPUTS
% study_nm: study name
%
% OUTPUTS
% IL_data: structure with  interleukines data with the following subfields:
%   .CID: subject identification number
%   .IL1b: interleukin IL1-b
%   .IL6: interleukin IL6
%   .IL18: interleukin IL18
%
% N. Clairis - november 2022

%% working directories
which_pc = 'Lab';
switch which_pc
    case 'Lab'
        gitPath = 'C:\Users\clairis\Desktop\';
    case 'home'
        gitPath = 'C:\Users\Loco\Documents\';
end 
ILpath = fullfile(gitPath,'GitHub','LGC_motiv',...
    'LGC_Motiv_results',study_nm,'Interleukins');
%% load data
% excelReadTable = readtable([ILpath,filesep,'saliva_IL_samples - results elisa031122.xlsx'],...
%     'Sheet','sample C.1');
excelReadTable = readtable([ILpath,filesep,'saliva_IL_samples.xlsx'],...
    'Sheet','sample C.1');

%% extract the data
IL_data.CID = excelReadTable.CID;
IL_data.IL1b = excelReadTable.Ilb_pg_mL_;
IL_data.IL6 = excelReadTable.IL6_pg_mL_;
IL_data.IL18 = excelReadTable.IL18_pg_mL_;

end % function