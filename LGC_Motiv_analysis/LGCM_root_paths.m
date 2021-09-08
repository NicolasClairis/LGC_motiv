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
rootPathList = {[fullfile('C:','Users','clairis','Desktop'),filesep],...
    [fullfile('C:','Users','Loco','Downloads'),filesep],...
    [fullfile('Arthur_ordi'),filesep],...
    [fullfile('Catherine_ordi'),filesep]};
rootListIdx = listdlg('ListString',rootPathList);
computerRoot = rootPathList{rootListIdx};

%% path where SPM anatomical template is (for preprocessing)
spmPathList = {fullfile('C:','Users','clairis','Desktop'),...
    fullfile('C:','Program Files','MATLAB','spm'),...
    fullfile('Arthur_ordi'),...
    fullfile('Catherine_ordi')};
spmListIdx = listdlg('ListString',spmPathList);
spmFolderPath = spmPathList{spmListIdx};

end % function