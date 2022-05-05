function[subj_scan_folders_names, jRun] = First_level_subRunFilter(study_nm, sub_nm,....
    subj_scan_folders_names, iRun, condition)
%[subj_scan_folders_names, jRun] = First_level_subRunFilter(study_nm, sub_nm,...
%   subj_scan_folders_names, iRun, condition)
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
% condition:
% 'fMRI': all fMRI compatible data
% 'fMRI_no_move': remove runs with too much movement
%
% OUTPUTS
% subj_scan_folders_names: cell with list of fMRI runs after filtering useless runs
%
% jRun: index corrected depending on the runs that are ok

%% check if subject has been well taken into account
[runs, n_runs] = runs_definition(study_nm, sub_nm, condition);
% check for 'fMRI' condition (ie all fMRI data)
if n_runs < 4 &&...
        (strcmp(condition,'fMRI') || strcmp(condition,'fMRI_no_move_bis')) &&...
        ( strcmp(study_nm,'study1') &&...
        ~ismember(sub_nm,{'017','040','074'}) ) ||...
        (strcmp(study_nm,'fMRI_pilots') &&...
        ~ismember(sub_nm,{'pilot_s1','pilot_s2'}) )
    error(['study ',study_nm,' subject ',sub_nm,': runs to include/exclude ',...
        'should be added in First_level_subRunFilter. ',...
        'Indeed there are only ',num2str(n_runs),' runs for this subject.']);
end

% check for 'fMRI_no_move' condition
if n_runs < 4 &&...
        strcmp(condition,'fMRI_no_move') &&...
        ( strcmp(study_nm,'study1') &&...
        ~ismember(sub_nm,{'008','017','018','022','029','040','044',...
        '047','054','065','071','074','076','087','093'}) ) ||...
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
            % in any case, you should remove run 1 from subject 017 and 074
            % and run 3 from subject 040 because fMRI crashed and/or
            % subject crashed the fMRI
            switch sub_nm
                case {'017'} % remove run 1 (not enough trials because fMRI crashed)
                    sub017_run1_toRemove = '2_007_run1_20220104';
                    if strcmp(subj_scan_folders_names(1,:),sub017_run1_toRemove)
                        subj_scan_folders_names(1,:) = [];
                        disp(['run1 named ',sub017_run1_toRemove,' got removed for CID',sub_nm]);
                    else
                        error(['file corresponding to run 1 could not be identified for CID',sub_nm,...
                            '. please fix it and remove it before going further in the analysis.']);
                    end
                case {'040'} % remove runs 3 and 4 (run3 crashed in the middle and run 4 never executed)
                    sub040_run3_toRemove = '2_009_run3_20220304';
                    if strcmp(subj_scan_folders_names(3,:),sub040_run3_toRemove)
                        subj_scan_folders_names(3,:) = [];
                        disp(['run3 named ',sub040_run3_toRemove,' got removed for CID',sub_nm]);
                    else
                        error(['file corresponding to run 3 could not be identified for CID',sub_nm,...
                            '. please fix it and remove it before going further in the analysis.']);
                    end
                case {'074'} % remove run 1 (not enough trials because fMRI crashed)
                    sub074_run1_toRemove = '3_007_run1_20211102';
                    if strcmp(subj_scan_folders_names(1,:),sub074_run1_toRemove)
                        subj_scan_folders_names(1,:) = [];
                        disp(['run1 named ',sub074_run1_toRemove,' got removed for CID',sub_nm]);
                    else
                        error(['file corresponding to run 1 could not be identified for CID',sub_nm,...
                            '. please fix it and remove it before going further in the analysis.']);
                    end
            end % subject
            
            %% filter runs with too much movement
            if strcmp(condition,'fMRI_no_move')
                switch sub_nm
                    case '018'
                        sub018_run3_toRemove = '2_009_run3_20220426';
                        if strcmp(subj_scan_folders_names(3,:),sub018_run3_toRemove)
                            subj_scan_folders_names(3,:) = [];
                            disp(['run3 named ',sub018_run3_toRemove,' got removed for CID',sub_nm,...
                                ' because of too much movement.']);
                        else
                            error(['file corresponding to run 3 could not be identified for CID',sub_nm,...
                                '. please fix it and remove it before going further in the analysis.']);
                        end
                    case '029'
                        sub029_run4_toRemove = '2_010_run4_20220310';
                        if strcmp(subj_scan_folders_names(4,:),sub029_run4_toRemove)
                            subj_scan_folders_names(4,:) = [];
                            disp(['run4 named ',sub029_run4_toRemove,' got removed for CID',sub_nm,...
                                ' because of too much movement.']);
                        else
                            error(['file corresponding to run 4 could not be identified for CID',sub_nm,...
                                '. please fix it and remove it before going further in the analysis.']);
                        end
                    case '044'
                        sub044_run2_toRemove = '3_008_run2_20220408';
                        sub044_run4_toRemove = '3_010_run4_20220408';
                        % remove run 4 first
                        if strcmp(subj_scan_folders_names(4,:),sub044_run4_toRemove)
                            subj_scan_folders_names(4,:) = [];
                            disp(['run4 named ',sub044_run4_toRemove,' got removed for CID',sub_nm,...
                                ' because of too much movement.']);
                        else
                            error(['file corresponding to run 4 could not be identified for CID',sub_nm,...
                                '. please fix it and remove it before going further in the analysis.']);
                        end
                        % remove run 2 secondly
                        if strcmp(subj_scan_folders_names(2,:),sub044_run2_toRemove)
                            subj_scan_folders_names(2,:) = [];
                            disp(['run2 named ',sub044_run2_toRemove,' got removed for CID',sub_nm,...
                                ' because of too much movement.']);
                        else
                            error(['file corresponding to run 2 could not be identified for CID',sub_nm,...
                                '. please fix it and remove it before going further in the analysis.']);
                        end
                    case '047'
                        sub047_run3_toRemove = '2_009_run3_20220311';
                        if strcmp(subj_scan_folders_names(3,:),sub047_run3_toRemove)
                            subj_scan_folders_names(3,:) = [];
                            disp(['run3 named ',sub047_run3_toRemove,' got removed for CID',sub_nm,...
                                ' because of too much movement.']);
                        else
                            error(['file corresponding to run 3 could not be identified for CID',sub_nm,...
                                '. please fix it and remove it before going further in the analysis.']);
                        end
                    case '054'
                        sub054_run2_toRemove = '2_008_run2_20220208';
                        sub054_run4_toRemove = '2_010_run4_20220208';
                        % remove run 4 first
                        if strcmp(subj_scan_folders_names(4,:),sub054_run4_toRemove)
                            subj_scan_folders_names(4,:) = [];
                            disp(['run4 named ',sub054_run4_toRemove,' got removed for CID',sub_nm,...
                                ' because of too much movement.']);
                        else
                            error(['file corresponding to run 4 could not be identified for CID',sub_nm,...
                                '. please fix it and remove it before going further in the analysis.']);
                        end
                        % remove run 2 secondly
                        if strcmp(subj_scan_folders_names(2,:),sub054_run2_toRemove)
                            subj_scan_folders_names(2,:) = [];
                            disp(['run2 named ',sub054_run2_toRemove,' got removed for CID',sub_nm,...
                                ' because of too much movement.']);
                        else
                            error(['file corresponding to run 2 could not be identified for CID',sub_nm,...
                                '. please fix it and remove it before going further in the analysis.']);
                        end
                    case '065'
                        sub065_run3_toRemove = '2_009_run3_20220119';
                        if strcmp(subj_scan_folders_names(3,:),sub065_run3_toRemove)
                            subj_scan_folders_names(3,:) = [];
                            disp(['run3 named ',sub065_run3_toRemove,' got removed for CID',sub_nm,...
                                ' because of too much movement.']);
                        else
                            error(['file corresponding to run 3 could not be identified for CID',sub_nm,...
                                '. please fix it and remove it before going further in the analysis.']);
                        end
                    case '071'
                        sub071_run2_toRemove = '2_008_run2_20220405';
                        sub071_run4_toRemove = '2_010_run4_20220405';
                        % remove run 4 first
                        if strcmp(subj_scan_folders_names(4,:),sub071_run4_toRemove)
                            subj_scan_folders_names(4,:) = [];
                            disp(['run4 named ',sub071_run4_toRemove,' got removed for CID',sub_nm,...
                                ' because of too much movement.']);
                        else
                            error(['file corresponding to run 4 could not be identified for CID',sub_nm,...
                                '. please fix it and remove it before going further in the analysis.']);
                        end
                        % remove run 2 secondly
                        if strcmp(subj_scan_folders_names(2,:),sub071_run2_toRemove)
                            subj_scan_folders_names(2,:) = [];
                            disp(['run2 named ',sub071_run2_toRemove,' got removed for CID',sub_nm,...
                                ' because of too much movement.']);
                        else
                            error(['file corresponding to run 2 could not be identified for CID',sub_nm,...
                                '. please fix it and remove it before going further in the analysis.']);
                        end
                    case '076'
                        sub076_run1_toRemove = '2_007_run1_20220420';
                        if strcmp(subj_scan_folders_names(1,:),sub076_run1_toRemove)
                            subj_scan_folders_names(1,:) = [];
                            disp(['run1 named ',sub076_run1_toRemove,' got removed for CID',sub_nm,...
                                ' because of too much movement.']);
                        else
                            error(['file corresponding to run 1 could not be identified for CID',sub_nm,...
                                '. please fix it and remove it before going further in the analysis.']);
                        end
                    case '087'
                        sub087_run2_toRemove = '2_008_run2_20211215';
                        sub087_run3_toRemove = '2_009_run3_20211215';
                        sub087_run4_toRemove = '2_010_run4_20211215';
                        % remove run 4 first
                        if strcmp(subj_scan_folders_names(4,:),sub087_run4_toRemove)
                            subj_scan_folders_names(4,:) = [];
                            disp(['run4 named ',sub087_run4_toRemove,' got removed for CID',sub_nm,...
                                ' because of too much movement.']);
                        else
                            error(['file corresponding to run 4 could not be identified for CID',sub_nm,...
                                '. please fix it and remove it before going further in the analysis.']);
                        end
                        % remove run 3 first
                        if strcmp(subj_scan_folders_names(3,:),sub087_run3_toRemove)
                            subj_scan_folders_names(3,:) = [];
                            disp(['run3 named ',sub087_run3_toRemove,' got removed for CID',sub_nm,...
                                ' because of too much movement.']);
                        else
                            error(['file corresponding to run 3 could not be identified for CID',sub_nm,...
                                '. please fix it and remove it before going further in the analysis.']);
                        end
                        % remove run 2 secondly
                        if strcmp(subj_scan_folders_names(2,:),sub087_run2_toRemove)
                            subj_scan_folders_names(2,:) = [];
                            disp(['run2 named ',sub087_run2_toRemove,' got removed for CID',sub_nm,...
                                ' because of too much movement.']);
                        else
                            error(['file corresponding to run 2 could not be identified for CID',sub_nm,...
                                '. please fix it and remove it before going further in the analysis.']);
                        end
                    case '093'
                        sub093_run1_toRemove = '2_007_run1_20220309';
                        if strcmp(subj_scan_folders_names(1,:),sub093_run1_toRemove)
                            subj_scan_folders_names(1,:) = [];
                            disp(['run1 named ',sub093_run1_toRemove,' got removed for CID',sub_nm,...
                                ' because of too much movement.']);
                        else
                            error(['file corresponding to run 1 could not be identified for CID',sub_nm,...
                                '. please fix it and remove it before going further in the analysis.']);
                        end
                end % subject
            end % filter runs with too much movement
    end

else
    subj_scan_folders_names = [];

end % subj_scan_folders_names correction

%% fix index of the run
if exist('iRun','var') && ~isempty(iRun)
    jRun = runs.runsToKeep(iRun);
end

end % function