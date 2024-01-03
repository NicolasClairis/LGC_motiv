function[Lac] = load_plasma_Lac(subject_id)
% [Lac] = load_plasma_Lac(subject_id)
% load_plasma_Lac will load the Lactate data from the plasma.
%
% INPUTS
% subject_id: list of subjects to include (all subjects by default if left
% empty)
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
if ~exist('subject_id','var') || isempty(subject_id)
    Lac.CID = excelReadTable.CID';
    Lac.Lac = excelReadTable.LacticAcid';
else
    NS = length(subject_id);
    Lac.CID = cell(1,NS);
    Lac.Lac = NaN(1,NS);
    for iS = 1:NS
        sub_nm = subject_id{iS};
        Lac.CID{iS} = sub_nm;
        sub_idx = find(strcmp(excelReadTable.CID, ['CID',sub_nm]));
        if ~isempty(sub_idx) && size(sub_idx,2) == 1
            Lac.Lac(iS) = excelReadTable.LacticAcid(sub_idx);
        else
            error(['Problem with Lactate extraction in subject ',sub_nm]);
        end
    end % subject loop
end

end % function