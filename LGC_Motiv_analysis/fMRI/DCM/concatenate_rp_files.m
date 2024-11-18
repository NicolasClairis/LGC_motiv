function[mvmtFolder, rp_file_nm] = concatenate_rp_files(subj_scans_folder, subj_scan_folders_names, n_scans_PerRun)
%[mvmtFolder, rp_file_nm] = concatenate_rp_files(subj_scans_folder, subj_scan_folders_names, n_scans_PerRun)
% concatenate_rp_files will concatenate all rp files from all the sessions
% to consider in the current GLM together and it will create a new file
% inside the "mvmtFolder" folder with the "rp_file_nm" name.
%
% INPUTS
% subj_scans_folder: path where subject scans are located
%
% subj_scan_folders_names: cell with nScans rows containing the name of the
% folder corresponding to each of the nScans to concatenate
%
% n_scans_PerRun: 1*nScans vector with number of scans for each fMRI run
%
% OUTPUTS
% mvmtFolder: folder where new file will be stored
%
% rp_file_nm: name of the new rp file
%
% See also First_level_btach_concat_DCM.m

%% extract number of scans
nRuns = size(subj_scan_folders_names,1);
n_totalScans = sum(n_scans_PerRun);
n_mvmt_prm = 6;

%% initialize variable of interest
rp_matrix = NaN(n_totalScans, n_mvmt_prm);

%% loop over runs to load and concatenate the rp files together
for iRun = 1:nRuns
    scan_idx = (1:n_scans_PerRun(iRun)) + sum(n_scans_PerRun(1:(iRun-1)));
    % erase useless spaces from folder with run name
    subj_runFoldername_tmp = strrep(subj_scan_folders_names(iRun, :),' ','');
    % load rp file of current run
    mvmtFolder_tmp = [subj_scans_folder, filesep, subj_runFoldername_tmp, filesep];
    movement_file_tmp = ls([mvmtFolder_tmp, 'rp*']);
    movement_filePath_tmp = [mvmtFolder_tmp, movement_file_tmp];
    rp_data_tmp = readmatrix(movement_filePath_tmp);
    rp_matrix(scan_idx, :) = rp_data_tmp;
end % run loop

%% save resulting file
mvmtFolder = [subj_scans_folder,filesep,'mvmt_concatenated',filesep];
% create folder if does not exist yet
if ~exist(mvmtFolder,'dir')
    mkdir(mvmtFolder);
end
% create file with concatenated movement across all sessions
rp_file_nm = ['rp_file_pool_',num2str(nRuns),'runs.txt'];
writematrix(rp_matrix,[mvmtFolder, rp_file_nm]);

end % function