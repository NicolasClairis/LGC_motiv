function[CORT_data] = load_CORT(study_nm)
% [CORT_data] = load_CORT(study_nm)
% load_CORT will load the cortisol data.
%
% INPUTS
% study_nm: string with study name
%
% OUTPUTS
% CORT_data: structure with cortisol data with the following subfields:
%   .CID: subject identification number (NS subjects in total)
%   .CORT: 4*NS matrix with each timepoint for each subject
%   .AUCg: 1xNS vector with area under the curve total cortisol during the
%   experiment
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
CORTpath = fullfile(gitPath,'GitHub','LGC_motiv',...
    'LGC_Motiv_results',study_nm,'CORT');
%% load data
excelReadTable = readtable([CORTpath,filesep,'saliva_CORT_samples_study1.xlsx'],...
    'Sheet','Sheet1');

%% extract the data
CORT_data.CID = excelReadTable.CID';
CORT_data.AUCg = excelReadTable.AUCg';
CORT_data.CORT = [excelReadTable.A1CORT';...
    excelReadTable.BCORT';...
    excelReadTable.DCORT';...
    excelReadTable.ECORT'];

end % function