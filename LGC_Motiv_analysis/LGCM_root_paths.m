function [computerRoot, spmFolderPath] = LGCM_root_paths()
%[computerRoot, spmFolderPath] = LGCM_root_paths()
% LGCM_root_paths define the paths for the computer currently in use (since
% different paths are being used in each of our computers.
%
% OUTPUTS
% computerRoot: path where the behavioral and fMRI data are stored
%
% spmFolderPath: path where SPM12 is

%% path where to do the analysis
Nicolas_homePath = [fullfile('C:','Users','Loco','Downloads'),filesep];
Nicolas_labPath = [fullfile('C:','Users','clairis','Desktop'),filesep];
% human_serverPath = [fullfile('svfas5.epfl.ch','Sandi-Lab','human_data_private'),filesep];
human_serverPath = [fullfile('L:','human_data_private'),filesep];
Arthur_serverPath = [fullfile('svfas5.epfl.ch','Sandi-Lab','Arthur'),filesep];
Catherine_serverPath = [fullfile('to be defined'),filesep];
Catherine_persoPath = [fullfile('to be defined'),filesep];
rootPathList = {Nicolas_homePath,...
    Nicolas_labPath,...
    Arthur_serverPath,...
    human_serverPath,...
    Catherine_serverPath,...
    Catherine_persoPath};
rootListIdx = listdlg('ListString',rootPathList);
computerRoot = rootPathList{rootListIdx};

%% path where SPM anatomical template is (for preprocessing)
switch computerRoot
    case Nicolas_homePath
        spmFolderPath = fullfile('C:','Program Files','MATLAB','spm','spm12');
    case Nicolas_labPath
        spmFolderPath = fullfile('C:','Users','clairis','Desktop');
    case Arthur_serverPath
        spmFolderPath = fullfile('D:','Matlab extensions');
    case Catherine_serverPath
        spmFolderPath = fullfile('to be defined');
    case Catherine_persoPath
        spmFolderPath = fullfile('to be defined');
    case human_serverPath
        % need to define where SPM is in this case
        spmPathList = {fullfile('C:','Users','clairis','Desktop'),...
            fullfile('C:','Program Files','MATLAB','spm','spm12'),...
            fullfile('D:','Matlab extensions'),...
            fullfile('Catherine_ordi')};
        spmListIdx = listdlg('ListString',spmPathList);
        spmFolderPath = spmPathList{spmListIdx};
end


end % function