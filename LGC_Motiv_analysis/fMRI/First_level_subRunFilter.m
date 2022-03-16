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

%% check if subject has been well taken into account
[~, n_runs] = runs_definition(study_nm, sub_nm, 'fMRI');
if n_runs < 4 &&...
        ( strcmp(study_nm,'study1') &&...
        ~ismember(sub_nm,{'017','040','074'}) ) ||...
        (strcmp(study_nm,'fMRI_pilots') &&...
        ~ismember(sub_nm,{'pilot_s1','pilot_s2'}) )
    error(['study ',study_nm,' subject ',sub_nm,': runs to include/exclude ',...
        'should be added in First_level_subRunFilter. ',...
        'Indeed there are only ',num2str(n_runs),' runs for this subject.']);
end

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
                case {'074'} % remove run 1 (not enough trials because fMRI crashed)
                    if strcmp(subj_scan_folders_names(1,:),'3_007_run1_20211102')
                        subj_scan_folders_names(1,:) = [];
                        disp(['run1 named ',subj_scan_folders_names(1,:),' got removed for CID',sub_nm]);
                    else
                        error(['file corresponding to run 1 could not be identified for CID',sub_nm,...
                            '. please fix it and remove it before going further in the analysis.']);
                    end
                case {'017'} % remove run 1 (not enough trials because fMRI crashed)
                    if strcmp(subj_scan_folders_names(1,:),'2_007_run1_20220104')
                        subj_scan_folders_names(1,:) = [];
                        disp(['run1 named ',subj_scan_folders_names(1,:),' got removed for CID',sub_nm]);
                    else
                        error(['file corresponding to run 1 could not be identified for CID',sub_nm,...
                            '. please fix it and remove it before going further in the analysis.']);
                    end
                case {'040'} % remove runs 3 and 4 (run3 crashed in the middle and run 4 never executed)
                    if strcmp(subj_scan_folders_names(1,:),'2_009_run3_20220304')
                        subj_scan_folders_names(1,:) = [];
                        disp(['run3 named ',subj_scan_folders_names(1,:),' got removed for CID',sub_nm]);
                    else
                        error(['file corresponding to run 3 could not be identified for CID',sub_nm,...
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
                case {'074','017'} % adapt index since run 1 has to be ignored
                    jRun = iRun + 1;
                otherwise % when everything worked, the index is the same
                    jRun = iRun;
            end
    end
else
    jRun = [];
end

end % function