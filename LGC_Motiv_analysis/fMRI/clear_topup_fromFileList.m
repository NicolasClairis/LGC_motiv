function[subj_scan_folders_names] = clear_topup_fromFileList(subj_scan_folders_names)
% [subj_scan_folders_names] = clear_topup_fromFileList(subj_scan_folders_names)
% clear_topup_fromFileList will remove topup folder names from subj_scan_folders_names
%
% INPUTS
% subj_scan_folders_names: list of all folder names including topups
%
% OUTPUTS
% subj_scan_folders_names: same but after removing topup files
%

% remove AP/PA corrective runs
for iRunCorrect = size(subj_scan_folders_names,1):-1:1
    % delete references from the list (made for preprocessing with AP/PA correction of distorsions)
    if strcmp(subj_scan_folders_names(iRunCorrect,end-11:end-8),'_PA_') ||...
            strcmp(subj_scan_folders_names(iRunCorrect,end-13:end-10),'_PA_') ||...
                strcmp(subj_scan_folders_names(iRunCorrect,11:17),'_topup_')
        subj_scan_folders_names(iRunCorrect,:) = [];
    end
    
end

end % function