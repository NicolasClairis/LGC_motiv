function[subj_scan_folders_names_fixed] = badRunsClearing(run_nm_toRemove, subj_scan_folders_names, sub_nm, reason)
% [subj_scan_folders_names] = badRunsClearing(run_nm_toRemove, subj_scan_folders_names, sub_nm, reason)
% badRunsClearing will remove the runs that have to be removed from the
% list
%
% INPUTS
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

%% remove bad runs
subj_scan_folders_names_fixed = subj_scan_folders_names;
n_badRuns = size(run_nm_toRemove,2); % if more than 1 bad run: procede in descending order
j_runs_removed = 0;
for iRunToRmv = 1:n_badRuns
    runToRmv_fullName = run_nm_toRemove{iRunToRmv};
    nRuns = size(subj_scan_folders_names_fixed,1);

    for iRun = nRuns:(-1):1
        if strcmp(subj_scan_folders_names_fixed(iRun,:), runToRmv_fullName)
            subj_scan_folders_names_fixed(iRun,:) = [];
            disp(['run named ',runToRmv_fullName,' got removed for CID',sub_nm,...
                ' because of ',reason,'.']);
            j_runs_removed = j_runs_removed + 1;
        end % remove bad runs
    end % run loop on all runs to find the one to remove
end % loop on runs to remove

%% check if all went well
if j_runs_removed ~= n_badRuns
    error(['Some runs should have been removed but could not be found for CID',sub_nm,...
        '. Please check what happened.']); 
end
end % function