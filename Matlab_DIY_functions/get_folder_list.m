function [ outputList ] = get_folder_list(expr, path )
%[ outputList ] = get_folder_list(expr, path )
% get_folder_list outputs all the folder containing the expression given in
% 'expr' inside the 'path' path (or current path by default) inside
% outputList.
%
% INPUTS
% expr: string with expression to seek (you should add the stars where
% necessary)
%
% path: path (current by default if left empty)
%
% OUTPUTS
% outputList: string list of the folders containing the expr expression
%
% Created by N.Clairis - june 2019

if ~exist('path','var') || isempty(path)
    path = pwd;
end

%% get all the folders with expr in their name
A = dir([path, filesep, expr]);

n_folders = size(A,1);
outputList = cell(n_folders, 1);

for iFolder = 1:n_folders
    outputList{iFolder} = A(iFolder).name; 
end


end

