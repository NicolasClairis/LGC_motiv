function [computerRoot, spmFolderPath] = LGCM_root_paths()
%[computerRoot, spmFolderPath] = LGCM_root_paths()
% LGCM_root_paths define the paths for the computer currently in use (since
% different paths are being used in each of our computers.
%
% OUTPUTS
% computerRoot: path where the behavioral and fMRI data are stored
%
% spmFolderPath: path where SPM12 is

%% define possible paths
% Nicolas_homePath = [fullfile('C:','Users','Loco','Downloads'),filesep];
% Nicolas_labPath = [fullfile('C:','Users','clairis','Desktop'),filesep];
Nicolas_labPath = [fullfile('E:'),filesep];
human_serverPath = [filesep,filesep,fullfile('sv-nas1.rcp.epfl.ch','Sandi-lab','human_data_private','raw_data_subject'),filesep];
% human_serverPath = [fullfile('L:','human_data_private','raw_data_subject'),filesep];
Arthur_serverPath = [fullfile('sv-nas1.rcp.epfl.ch','Sandi-Lab','Arthur'),filesep];
% Nicolas_homePath = human_serverPath;
Nicolas_homePath = [fullfile('F:'),filesep]; % path to hard drive with data

% spm specific
Nicolas_homeSPMpath = fullfile('C:','Program Files','MATLAB','spm');
Nicolas_labSPMpath = fullfile('C:','Users','clairis','Desktop');
Arthur_SPMpath = fullfile('D:','Matlab extensions');

%% path where results are stored
rootPathList = {Nicolas_labPath,...
    Nicolas_homePath,...
    Arthur_serverPath,...
    human_serverPath};
rootListIdx = listdlg('PromptString','results path',...
    'SelectionMode','single',...
    'ListString',rootPathList);
computerRoot = rootPathList{rootListIdx};

%% path where SPM anatomical template is (for preprocessing)
switch computerRoot
    case Nicolas_homePath
        spmFolderPath = Nicolas_homeSPMpath;
    case Nicolas_labPath
        spmFolderPath = Nicolas_labSPMpath;
    case Arthur_serverPath
        spmFolderPath = Arthur_SPMpath;
    case human_serverPath
        % need to define where SPM is in this case
        spmPathList = {Nicolas_labSPMpath,...
            Nicolas_homeSPMpath,...
            Arthur_SPMpath};
        spmListIdx = listdlg('PromptString','SPM path',...
            'SelectionMode','single',...
            'ListString',spmPathList);
        spmFolderPath = spmPathList{spmListIdx};
end


end % function