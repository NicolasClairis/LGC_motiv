function[TESTO_data] = load_TESTO(study_nm)
% [TESTO_data] = load_TESTO(study_nm)
% load_TESTO will load the testosterone data.
%
% INPUTS
% study_nm: string with study name
%
% OUTPUTS
% TESTO_data: structure with testosterone data with the following subfields:
%   .CID: subject identification number (NS subjects in total)
%   .TESTO: 4*NS matrix with each timepoint for each subject
%   .AUCg: 1xNS vector with area under the curve total testosterone during the
%   experiment
%   .timings: 4*NS matrix with timing of each timepoint for each subject
%
% N. Clairis - april 2023

%% working directories
which_pc = 'home';
switch which_pc
    case 'Lab'
        gitPath = 'C:\Users\clairis\Desktop\';
    case 'home'
        gitPath = 'C:\Users\Loco\Documents\';
end 
TESTOpath = fullfile(gitPath,'GitHub','LGC_motiv',...
    'LGC_Motiv_results',study_nm,'TESTO');
%% load data
excelReadTable = readtable([TESTOpath,filesep,'results_testosterone_2023_04_03.xlsx'],...
    'Sheet','results_to_load');

%% extract the data
TESTO_data.CID = excelReadTable.CID';
TESTO_data.AUCg = excelReadTable.AUCg';
TESTO_A1 = excelReadTable.A_1';
TESTO_B = excelReadTable.B';
TESTO_D = excelReadTable.D';
TESTO_E = excelReadTable.E';
if iscell(TESTO_B)
    TESTO_B = str2double(TESTO_B);
end
if iscell(TESTO_D)
    TESTO_D = str2double(TESTO_D);
end
TESTO_data.TESTO = [TESTO_A1;...
    TESTO_B;...
    TESTO_D;...
    TESTO_E];
TESTO_data.timings = [excelReadTable.SalivaA1Cortisol';...
    excelReadTable.SalivaBCortisol_post_MRS_';...
    excelReadTable.SalivaDCortisol_pr__IRMf_';...
    excelReadTable.SalivaECortisol_post_IRMf_'].*24; % store timings in hours

end % function