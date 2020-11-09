function [outputList] = ls_regexp(regexpPattern, dirPath)
%function [] = ls_regexp(regexpPattern, dirPath)
%
% This function returns a list of files or folders matching the
% regexpPattern in dirPath folder (or current folder if dirPath not
% specified)
%
% INPUTS
% regexpPattern: 
%
% dirPath: folder where to look (if left empty, checks in current folder)
%
% Written by E.Bioud - june 2018

    if nargin > 1
        currDir = pwd;
        cd(dirPath)
    end
    
    listContent = string(ls);
    
    startIdx = regexp(listContent, regexpPattern);
    containsIdx = ~cellfun(@isempty, startIdx);
    outputList = char(listContent(containsIdx));
    
    if nargin > 1
        cd(currDir)
    end
    
end