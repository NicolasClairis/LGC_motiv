function [] = check_currPath(folderName)
%check_currPath(folderName)
% check_currPath will check if you are currently in the folder containing 
% your functions of interest. If not, will display an error message
%
% INPUTS
% folderName: name of the current folder
%

currPath = pwd;
currPathSplit = split(currPath, filesep);
if ~strcmp(currPathSplit{end,1}, folderName)
    error('you should launch the script within the LGC_Motiv_task folder');
end

end % function