function[subj_scan_folders_names_fixed] = badRunsClearing(nRunToRemove, run_nm_toRemove, subj_scan_folders_names, sub_nm, reason)
% [subj_scan_folders_names] = badRunsClearing(nRunToRemove, run_nm_toRemove, subj_scan_folders_names, sub_nm, reason)
% badRunsClearing will remove the runs that have to be removed from the
% list
%
% INPUTS
% nRunsToRemove: 1*n vector with index of the runs to remove
%
% run_nm_toRemove: cell with list of full names of runs to remove
%
% subj_scan_folders_names: cell with list of subject scan folder names for
% each run
%
% sub_nm: string with subject identification number ('XXX')
%
% reason: reason why the runs have to excluded
%
% OUTPUT
% subj_scan_folders_names_fixed: cell with list of subject scan folder
% names for each run that has to be included after removing the bad ones


subj_scan_folders_names_fixed = subj_scan_folders_names;

n_badRuns = size(nRunToRemove,2); % if more than 1 bad run: procede in descending order
% check first if order is ok
if n_badRuns > 1
    isOrderOk = sum(nRunToRemove == sort(nRunToRemove,'descend')) == n_badRuns;
    if isOrderOk == false
        error(['Please order runs in descending order for subject ',sub_nm]);
    end
end
%% remove bad runs
for iRunToRmv = 1:n_badRuns
    runToRmvNb = nRunToRemove(iRunToRmv);
    runToRmvNb_nm = num2str(runToRmvNb);
    runToRmv_fullName = run_nm_toRemove{iRunToRmv};
    if strcmp(subj_scan_folders_names(runToRmvNb,:), runToRmv_fullName)
        subj_scan_folders_names_fixed(runToRmvNb,:) = [];
        disp(['run ',runToRmvNb_nm,' named ',runToRmv_fullName,' got removed for CID',sub_nm,...
            ' because of ',reason,'.']);
    else
        error(['file corresponding to run ',runToRmvNb_nm,' could not be identified for CID',sub_nm,...
            '. please fix it and remove it before going further in the analysis.']);
    end
end % loop on runs to remove
end % function