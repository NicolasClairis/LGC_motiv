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
        (strcmp(condition,'fMRI') || strcmp(condition,'fMRI_no_move_bis') || strcmp(condition,'fMRI_noSatTask')) &&...
        ( strcmp(study_nm,'study1') &&...
        ~ismember(sub_nm,{'017','040','043','074'}) ) ||...
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
            %% remove runs with acquisition problem
            % in any case, you should remove run 1 from subject 017 and 074
            % and run 3 from subject 040 because fMRI crashed and/or
            % subject crashed the fMRI
            switch sub_nm
                case {'017','040','043','074'}
                    removalReason = 'problem during the fMRI acquisition';
                    switch sub_nm
                        case {'017'} % remove run 1 (not enough trials because fMRI crashed)
                            nRunToRemove = 1;
                            run_nm_toRemove = {'2_007_run1_20220104'};
                        case {'040'} % remove runs 3 (and 4) (run3 crashed in the middle and run 4 never executed)
                            nRunToRemove = 3;
                            run_nm_toRemove = {'2_009_run3_20220304'};
                        case {'043'}
                            nRunToRemove = 1;
                            run_nm_toRemove = {'2_007_run1_20220531'};
                        case {'074'} % remove run 1 (not enough trials because fMRI crashed)
                            nRunToRemove = 1;
                            run_nm_toRemove = {'3_007_run1_20211102'};
                    end % subject
                    [subj_scan_folders_names] = badRunsClearing(nRunToRemove, run_nm_toRemove,...
                        subj_scan_folders_names, sub_nm, removalReason);
            end
            %% filter runs with behavioral saturation
            if strcmp(condition,'fMRI_noSatTask')
                removalReason = 'saturation';
                switch sub_nm
                    case {'002','005','012','027','032','047','048','052',...
                            '076','095','100'}
                        switch sub_nm
                            case '002'
                                nRunToRemove = 3;
                                run_nm_toRemove = {'2_009_run3_20220222'};
                            case '005'
                                nRunToRemove = 4;
                                run_nm_toRemove = {'2_010_run4_2022_0506'};
                            case '012'
                                nRunToRemove = 4;
                                run_nm_toRemove = {'2_010_run4_20221021'};
                            case '027'
                                nRunToRemove = [4,2];
                                run_nm_toRemove = {'2_010_run4_20220624',...
                                    '2_008_run2_20220624'};
                            case '032'
                                nRunToRemove = 3;
                                run_nm_toRemove = {'1_009_run3_20220428'};
                            case '047'
                                nRunToRemove = [4,2,1];
                                run_nm_toRemove = {'2_010_run4_20220311',...
                                    '2_008_run2_20220311',...
                                    '2_007_run1_20220311'};
                            case '048'
                                nRunToRemove = 2;
                                run_nm_toRemove = {'3_008_run2_20220518'};
                            case '052'
                                nRunToRemove = [3,1];
                                run_nm_toRemove = {'2_010_run3_20220216'};
                            case '076'
                                nRunToRemove = 4;
                                run_nm_toRemove = {'2_010_run4_20220420'};
                            case '095'
                                nRunToRemove = [4,2];
                                run_nm_toRemove = {'2_010_run4_20211119',...
                                    '2_008_run2_20211119'};
                            case '100'
                                nRunToRemove = [4,3];
                                run_nm_toRemove = {'2_010_run4_20220324',...
                                    '2_009_run3_20220324'};
                        end % subject filter
                        [subj_scan_folders_names] = badRunsClearing(nRunToRemove, run_nm_toRemove,...
                            subj_scan_folders_names, sub_nm, removalReason);
                end % subject filter
            end % filter for saturation
            
            %% filter runs with too much movement
            if strcmp(condition,'fMRI_no_move')
                removalReason = 'too much movement';
                switch sub_nm
                    case {'008','021','022','024','029',...
                            '044','047',...
                            '053','054','058','062',...
                            '071','076','078','080','083','087',...
                            '097','099'}
                        switch sub_nm
                            case {'008','022','024'} % too much movement ALL runs
                                error(['subject ',sub_nm,'should not be included at all']);
                            case '021'
                                nRunToRemove = [4,3,2];
                                run_nm_toRemove = {'2_010_run4_20220519',...
                                    '2_009_run3_20220519',...
                                    '2_008_run2_20220519'};
                            case '029'
                                nRunToRemove = 4;
                                run_nm_toRemove = {'2_010_run4_20220310'};
                            case '044'
                                nRunToRemove = [4,2];
                                run_nm_toRemove = {'3_010_run4_20220408',...
                                    '3_008_run2_20220408'};
                            case '047'
                                nRunToRemove = 3;
                                run_nm_toRemove = {'2_009_run3_20220311'};
                            case '053'
                                nRunToRemove = 3;
                                run_nm_toRemove = {'2_009_run3_20220727'};
                            case '054'
                                nRunToRemove = [4,2];
                                run_nm_toRemove = {'2_010_run4_20220208',...
                                    '2_008_run2_20220208'};
                            case '058'
                                nRunToRemove = [4,2];
                                run_nm_toRemove = {'3_010_run4_20220621',...
                                    '3_008_run2_20220621'};
                            case '062'
                                nRunToRemove = [4,3,2];
                                run_nm_toRemove = {'2_010_run4_20220722',...
                                    '2_009_run3_20220722',...
                                    '2_008_run2_20220722'};
                            case '071'
                                nRunToRemove = [4,2];
                                run_nm_toRemove = {'2_010_run4_20220405',...
                                    '2_008_run2_20220405'};
                            case '076'
                                nRunToRemove = 1;
                                run_nm_toRemove = {'2_007_run1_20220420'};
                            case '078'
                                nRunToRemove = [4,1];
                                run_nm_toRemove = {'2_010_run4_20220630',...
                                    '2_007_run1_20220630'};
                            case '080'
                                nRunToRemove = 3;
                                run_nm_toRemove = {'2_009_run3_20220629'};
                            case '083'
                                nRunToRemove = 3;
                                run_nm_toRemove = {'2_009_run3_20220520'};
                            case '087'
                                nRunToRemove = [4,3,2];
                                run_nm_toRemove = {'2_010_run4_20211215',...
                                    '2_009_run3_20211215',...
                                    '2_008_run2_20211215'};
                            case '097'
                                nRunToRemove = [4,3,2];
                                run_nm_toRemove = {'2_010_run4_20220721',...
                                    '2_009_run3_20220721',...
                                    '2_008_run2_20220721'};
                            case '099'
                                nRunToRemove = [4,2];
                                run_nm_toRemove = {'1_010_run4_20220705',...
                                    '1_008_run2_20220705'};
                        end % subject
                        [subj_scan_folders_names] = badRunsClearing(nRunToRemove, run_nm_toRemove,...
                            subj_scan_folders_names, sub_nm, removalReason);
                end
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