function[subj_scan_folders_names, jRun] = First_level_subRunFilter(study_nm, sub_nm, subj_scan_folders_names, iRun)
%[subj_scan_folders_names, jRun] = First_level_subRunFilter(study_nm, sub_nm, subj_scan_folders_names, iRun)
% First_level_subRunFilter will remove the files that cannot be used when
% appropriate and will adapt the run index accordingly.
% 
% INPUTS
% study_nm: string with study name ('fMRI_pilots'/'study1'/'study2')
%
% sub_nm: string with subject name
%
% sub_scan_folders_names: cell with list of fMRI runs
%
% iRun: index of the file checked now
%
% OUTPUTS
% subj_scan_folders_names: cell with list of fMRI runs after filtering useless runs
%
% jRun: index corrected depending on the runs that are ok

%% for run file names correction
if exist('subj_scan_folders_names','var') && ~isempty(subj_scan_folders_names)

    %% check if extraction worked
    if isempty(subj_scan_folders_names)
        error('The fMRI files could not be extracted. Check the name of the files maybe there is something wrong there.');
    end

    %% remove file names that are not useable
    switch study_nm
        case 'fMRI_pilots'
            switch sub_nm
                case {'pilot_s1','pilot_s2'} % remove top-up files (=AP/PA corrective runs) when they were performed
                    % (only 2 for the two first fMRI pilots)
                    [subj_scan_folders_names] = clear_topup_fromFileList(subj_scan_folders_names);
            end
        case 'study1'
            switch sub_nm
                case 'CID074' % remove run 1 (not enough trials)
                    if strcmp(subj_scan_folders_names(1,:),'3_007_run1_20211102')
                        subj_scan_folders_names(1,:) = [];
                        disp(['run1 named 3_007_run1_20211102 got removed for ',sub_nm]);
                    else
                        error(['file corresponding to run 1 could not be identified for ',sub_nm,...
                            '. please fix it and remove it before going further in the analysis.']);
                    end
            end
    end

else
    subj_scan_folders_names = [];

end % subj_scan_folders_names correction

%% fix index of the run
if exist('iRun','var') && ~isempty(iRun)
    switch study_nm
        case 'study1'
            switch sub_nm
                case 'CID074' % adapt index since run 1 has to be ignored
                    jRun = iRun + 1;
                otherwise % when everything worked, the index is the same
                    jRun = iRun;
            end
    end
else
    jRun = [];
end

end % function