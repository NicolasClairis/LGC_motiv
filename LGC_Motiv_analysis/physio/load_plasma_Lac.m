function[Lac] = load_plasma_Lac()
% [Lac] = load_plasma_Lac()
% load_plasma_Lac will load the Lactate data from the plasma.
%
% INPUTS
% study_nm: string with study name
%
% OUTPUTS
% Lac: structure with Lac data with the following subfields:
%   .CID: subject identification number (NS subjects in total)
%   .Lac: concentration of lactate in the plasma
%
% N. Clairis - november 2023

%% working directories
plasma_path = 'P:\boulot\postdoc_CarmenSandi\results\plasma';
%% load data
excelReadTable = readtable([plasma_path,filesep,'plasma_blood_results.xlsx'],...
    'Sheet','simplified');

%% extract the data
Lac.CID = excelReadTable.CID';
Lac.Lac = excelReadTable.LacticAcid';

end % function