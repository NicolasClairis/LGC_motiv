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
Nicolas_homePath = [fullfile('C:','Users','Loco','Downloads'),filesep];
% Nicolas_labPath = [fullfile('C:','Users','clairis','Desktop'),filesep];
Nicolas_labPath = [fullfile('E:'),filesep];
% human_serverPath = [fullfile('svfas5.epfl.ch','Sandi-Lab','human_data_private'),filesep];
human_serverPath = [fullfile('L:','human_data_private'),filesep];
Arthur_serverPath = [fullfile('svfas5.epfl.ch','Sandi-Lab','Arthur'),filesep];
Catherine_serverPath = [fullfile('to be defined'),filesep];
Catherine_persoPath = [fullfile('to be defined'),filesep];

% spm specific
Nicolas_homeSPMpath = fullfile('C:','Program Files','MATLAB','spm');
Nicolas_labSPMpath = fullfile('C:','Users','clairis','Desktop');
Arthur_SPMpath = fullfile('D:','Matlab extensions');

%% path where github folder scripts are stored
rootPathList = {Nicolas_labPath,...
    Nicolas_homePath,...
    Arthur_serverPath,...
    human_serverPath,...
    Catherine_serverPath,...
    Catherine_persoPath};
rootListIdx = listdlg('ListString',rootPathList);
computerRoot = rootPathList{rootListIdx};

%% path where SPM anatomical template is (for preprocessing)
switch computerRoot
    case Nicolas_homePath
        spmFolderPath = Nicolas_homeSPMpath;
    case Nicolas_labPath
        spmFolderPath = Nicolas_labSPMpath;
    case Arthur_serverPath
        spmFolderPath = Arthur_SPMpath;
    case Catherine_serverPath
        spmFolderPath = fullfile('to be defined');
    case Catherine_persoPath
        spmFolderPath = fullfile('to be defined');
    case human_serverPath
        % need to define where SPM is in this case
        spmPathList = {Nicolas_labSPMpath,...
            Nicolas_homeSPMpath,...
            Arthur_SPMpath,...
            fullfile('Catherine_ordi')};
        spmListIdx = listdlg('ListString',spmPathList);
        spmFolderPath = spmPathList{spmListIdx};
end


end % function