function[CORT_data] = load_CORT(study_nm, subject_id)
% [CORT_data] = load_CORT(study_nm, subject_id)
% load_CORT will load the cortisol data.
%
% INPUTS
% study_nm: string with study name
%
% subject_id: list of subjects to extract
%
% OUTPUTS
% CORT_data: structure with cortisol data with the following subfields:
%   .CID: subject identification number (NS subjects in total)
%   .CORT: 4*NS matrix with each timepoint for each subject
%   .AUCg: 1xNS vector with area under the curve total cortisol during the
%   experiment
%   .timings: 4*NS matrix with timing of each timepoint for each subject
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
CORT_data.timings = [excelReadTable.SalivaA1Cortisol';...
    excelReadTable.SalivaBCortisol_post_MRS_';...
    excelReadTable.SalivaDCortisol_pr__IRMf_';...
    excelReadTable.SalivaECortisol_post_IRMf_'].*24; % store timings in hours

%% filter based on subjects entered in input if subject_id is not empty
if exist('subject_id','var') && ~isempty(subject_id)
    NS = length(subject_id);
    sub_to_include = false(1,NS);
    for iS = 1:NS
        sub_nm = ['CID',subject_id{iS}];
        sub_idx = find(strcmp(sub_nm, CORT_data.CID));
        if ~isempty(sub_idx)
           sub_to_include(sub_idx) = true;
        end
    end % subject loop
    
    % apply filter on output
    CORT_data.CID = CORT_data.CID(sub_to_include);
    CORT_data.AUCg = CORT_data.AUCg(sub_to_include);
    CORT_data.CORT = CORT_data.CORT(:,sub_to_include);
    CORT_data.timings = CORT_data.timings(:,sub_to_include);
end % subject filter

end % function