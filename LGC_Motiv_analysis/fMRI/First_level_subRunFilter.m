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
% 'fMRI_noMoveSub': remove subjects with too much movement in ALL runs
%
% OUTPUTS
% subj_scan_folders_names: cell with list of fMRI runs after filtering useless runs
%
% jRun: index corrected depending on the runs that are ok

%% check if subject has been well taken into account
[runs, n_runs] = runs_definition(study_nm, sub_nm, condition);
% check for 'fMRI' condition (ie all fMRI data)
if n_runs < 4 &&...
        (ismember(condition,{'fMRI',...
        'fMRI_noMoveSub','fMRI_noMoveSub_bis','fMRI_noMoveSub_ter',...
        'fMRI_noSatTaskSub','fMRI_noSatRunSub',...
        'fMRI_noSatTaskSub_noMove_bis_Sub'})) &&...
        ( strcmp(study_nm,'study1') && ~ismember(sub_nm,{'017','040','043','074'}) ) ||...
        (strcmp(study_nm,'fMRI_pilots') && ~ismember(sub_nm,{'pilot_s1','pilot_s2'}) )
    error(['study ',study_nm,' subject ',sub_nm,': runs to include/exclude ',...
        'should be added in First_level_subRunFilter. ',...
        'Indeed there are only ',num2str(n_runs),' runs for this subject.']);
end

%% for run file names correction
if exist('subj_scan_folders_names','var') && ~isempty(subj_scan_folders_names)

    %% remove file names that are not useable
    switch study_nm
        case 'fMRI_pilots'
            switch sub_nm
                case {'pilot_s1','pilot_s2'} % remove top-up files (=AP/PA corrective runs) when they were performed
                    % (only 2 for the two first fMRI pilots)
                    [subj_scan_folders_names] = clear_topup_fromFileList(subj_scan_folders_names);
            end
        case 'study1'
            %% remove runs with acquisition problem independent of the condition
            % in any case, you should remove run 1 from subject 017 and 074
            % and run 3 from subject 040 because fMRI crashed and/or
            % subject crashed the fMRI
            switch sub_nm
                case {'017','040','043','074'}
                    removalReason = 'problem during the fMRI acquisition';
                    switch sub_nm
                        case {'017'} % remove run 1 (not enough trials because fMRI crashed)
                            run_nm_toRemove = {'2_007_run1_20220104'};
                        case {'040'} % remove runs 3 (and 4) (run3 crashed in the middle and run 4 never executed)
                            run_nm_toRemove = {'2_009_run3_20220304'};
                        case {'043'}
                            run_nm_toRemove = {'2_007_run1_20220531'};
                        case {'074'} % remove run 1 (not enough trials because fMRI crashed)
                            run_nm_toRemove = {'3_007_run1_20211102'};
                        otherwise
                            error([sub_nm,' is missing.'])
                    end % subject
                    [subj_scan_folders_names] = badRunsClearing(run_nm_toRemove,...
                        subj_scan_folders_names, sub_nm, removalReason);
            end

            %% remove runs depending on the condition entered in input
            %% remove tasks fully saturated
            switch condition
                case 'fMRI_noSatTask'
                    removalReason = 'task saturation';
                    switch sub_nm
                        case {'027','047','052','095'}
                            switch sub_nm
                                case '027'
                                    run_nm_toRemove = {'2_010_run4_20220624',...
                                        '2_008_run2_20220624'};
                                case '047'
                                    error(['subject ',sub_nm,' should not be included (saturated ALL tasks)']);
                                case '052'
                                    run_nm_toRemove = {'2_010_run3_20220216',...
                                        '2_008_run1_20220216'};
                                case '069'
                                    run_nm_toRemove = {'2_010_run4_20220414',...
                                        '2_008_run2_20220414'};
                                case '076'
                                    run_nm_toRemove = {'2_010_run4_20220420',...
                                        '2_008_run2_20220420'};
                                case '095'
                                    run_nm_toRemove = {'2_010_run4_20211119',...
                                        '2_008_run2_20211119'};
                                otherwise
                                    error([sub_nm,' is missing.'])
                            end % subject individual
                            [subj_scan_folders_names] = badRunsClearing(run_nm_toRemove,...
                                subj_scan_folders_names, sub_nm, removalReason);
                    end % subject group
                    %% remove any run with saturation
                case 'fMRI_noSatRun'
                    removalReason = 'saturation';
                    switch sub_nm
                        case {'002','004','005','012','022','027',...
                                '032','038','044','047','048',...
                                '052','054','055','058',...
                                '061','062','069',...
                                '076','081','082','083','088',...
                                '095','097','099','100'}
                            switch sub_nm
                                case '002'
                                    run_nm_toRemove = {'2_009_run3_20220222'};
                                case '004'
                                    run_nm_toRemove = {'2_009_run3_20220902'};
                                case '005'
                                    run_nm_toRemove = {'2_010_run4_20220506'};
                                case '012'
                                    run_nm_toRemove = {'2_010_run4_20221021'};
                                case '022'
                                    run_nm_toRemove = {'2_009_run3_20220303'};
                                case '027'
                                    run_nm_toRemove = {'2_010_run4_20220624',...
                                        '2_009_run3_20220624',...
                                        '2_008_run2_20220624'};
                                case '032'
                                    run_nm_toRemove = {'1_009_run3_20220428'};
                                case '038'
                                    run_nm_toRemove = {'3_008_run2_20220701'};
                                case '044'
                                    run_nm_toRemove = {'3_009_run3_20220408'};
                                case '047'
                                    error(['subject ',sub_nm,' should not be included (saturated ALL tasks)']);
                                case '048'
                                    run_nm_toRemove = {'3_008_run2_20220518'};
                                case '052'
                                    run_nm_toRemove = {'2_010_run3_20220216',...
                                        '2_008_run1_20220216'};
                                case '054'
                                    run_nm_toRemove = {'2_007_run1_20220208'};
                                case '055'
                                    run_nm_toRemove = {'2_010_run4_20220202'};
                                case '058'
                                    run_nm_toRemove = {'3_010_run4_20220621'};
                                case '061'
                                    run_nm_toRemove = {'2_007_run1_20211126'};
                                case '062'
                                    run_nm_toRemove = {'2_009_run3_20220722'};
                                case '069'
                                    run_nm_toRemove = {'2_010_run4_20220414',...
                                        '2_008_run2_20220414'};
                                case '076'
                                    run_nm_toRemove = {'2_010_run4_20220420',...
                                        '2_008_run2_20220420'};
                                case '081'
                                    run_nm_toRemove = {'2_009_run3_20220107'};
                                case '082'
                                    run_nm_toRemove = {'2_010_run4_20220308',...
                                        '2_009_run3_20220308'};
                                case '083'
                                    run_nm_toRemove = {'2_010_run4_20220520'};
                                case '088'
                                    run_nm_toRemove = {'2_010_run4_20220512'};
                                case '095'
                                    run_nm_toRemove = {'2_010_run4_20211119',...
                                        '2_009_run3_20211119',...
                                        '2_008_run2_20211119'};
                                case '097'
                                    run_nm_toRemove = {'2_007_run1_20220721'};
                                case '099'
                                    run_nm_toRemove = {'1_009_run3_20220705'};
                                case '100'
                                    run_nm_toRemove = {'2_010_run4_20220324',...
                                        '2_009_run3_20220324'};
                                otherwise
                                    error([sub_nm,' is missing.'])
                            end % subject filter
                            [subj_scan_folders_names] = badRunsClearing(run_nm_toRemove,...
                                subj_scan_folders_names, sub_nm, removalReason);
                    end % subject group
                    
                %% remove any run with saturation and if not enough variability when splitting depending on choice
                case 'fMRI_noSatRun_choiceSplit_Elvl'
                    removalReason = 'saturation';
                    switch sub_nm
                        case {'002','003','004','005','009',...
                                '012','015','018',...
                                '020','027',...
                                '032','036','038',...
                                '042','044','046','047','048',...
                                '052','053','054','055','059',...
                                '062','065','069',...
                                '074','076','078',...
                                '080','081','082','083',...
                                '095','097','100'}
                            switch sub_nm
                                case '002'
                                    run_nm_toRemove = {'2_011_run4_20220222',...
                                        '2_009_run3_20220222'};
                                case '003'
                                    run_nm_toRemove = {'4_011_run4_20220406'};
                                case '004'
                                    run_nm_toRemove = {'2_010_run4_20220902',...
                                        '2_009_run3_20220902'};
                                case '005'
                                    run_nm_toRemove = {'2_010_run4_20220506'};
                                case '009'
                                    run_nm_toRemove = {'2_010_run4_20220217'};
                                case '012'
                                    run_nm_toRemove = {'2_010_run4_20221021',...
                                        '2_009_run3_20221021',...
                                        '2_007_run1_20221021'};
                                case '015'
                                    run_nm_toRemove = {'2_008_run2_20220412',...
                                        '2_007_run1_20220412'};
                                case '018'
                                    run_nm_toRemove = {'2_010_run4_20220426'};
                                case '020'
                                    run_nm_toRemove = {'2_009_run3_20220126'};
                                case '022'
                                    run_nm_toRemove = {'2_009_run3_20220303'};
                                case '027'
                                    run_nm_toRemove = {'2_010_run4_20220624',...
                                        '2_008_run2_20220624'};
                                case '032'
                                    run_nm_toRemove = {'1_009_run3_20220428'};
                                case '036'
                                    run_nm_toRemove = {'2_010_run4_20211110',...
                                        '2_008_run2_20211110'};
                                case '038'
                                    run_nm_toRemove = {'3_008_run2_20220701'};
                                case '042'
                                    run_nm_toRemove = {'2_010_run4_20220527'};
                                case '044'
                                    run_nm_toRemove = {'3_009_run3_20220408'};
                                case '046'
                                    run_nm_toRemove = {'2_008_run2_20220128'};
                                case '047'
                                    run_nm_toRemove = {'2_010_run4_20220311',...
                                        '2_008_run2_20220311',...
                                        '2_007_run1_20220311'};
                                case '048'
                                    run_nm_toRemove = {'3_008_run2_20220518'};
                                case '052'
                                    run_nm_toRemove = {'2_010_run3_20220216',...
                                        '2_008_run1_20220216'};
                                case '053'
                                    run_nm_toRemove = {'2_010_run4_20220727'};
                                case '054'
                                    run_nm_toRemove = {'2_007_run1_20220208'};
                                case '055'
                                    run_nm_toRemove = {'2_010_run4_20220202',...
                                        '2_009_run3_20220202'};
                                case '058'
                                    run_nm_toRemove = {'3_010_run4_20220621'};
                                case '059'
                                    run_nm_toRemove = {'3_009_run3_20220628'};
                                case '062'
                                    run_nm_toRemove = {'2_009_run3_20220722'};
                                case '065'
                                    run_nm_toRemove = {'2_010_run4_20220119'};
                                case '069'
                                    run_nm_toRemove = {'2_010_run4_20220414',...
                                        '2_008_run2_20220414'};
                                case '074'
                                    run_nm_toRemove = {'3_009_run3_20211102'};
                                case '076'
                                    run_nm_toRemove = {'2_010_run4_20220420',...
                                        '2_008_run2_20220420'};
                                case '078'
                                    run_nm_toRemove = {'2_009_run3_20220630'};
                                case '080'
                                    run_nm_toRemove = {'2_010_run4_20220629'};
                                case '081'
                                    run_nm_toRemove = {'2_009_run3_20220107'};
                                case '082'
                                    run_nm_toRemove = {'2_010_run4_20220308'};
                                case '083'
                                    run_nm_toRemove = {'2_010_run4_20220520'};
                                case '088'
                                    run_nm_toRemove = {'2_010_run4_20220512'};
                                case '095'
                                    run_nm_toRemove = {'2_010_run4_20211119',...
                                        '2_009_run3_20211119',...
                                        '2_008_run2_20211119'};
                                case '097'
                                    run_nm_toRemove = {'2_007_run1_20220721'};
                                case '100'
                                    run_nm_toRemove = {'2_010_run4_20220324',...
                                        '2_009_run3_20220324'};
                                otherwise
                                    error([sub_nm,' is missing.'])
                            end % subject filter
                            [subj_scan_folders_names] = badRunsClearing(run_nm_toRemove,...
                                subj_scan_folders_names, sub_nm, removalReason);
                    end % subject group

                    %% remove any run with saturation and if not enough variability when splitting depending on choice
                case 'fMRI_noSatRun_choiceSplit_Elvl_bis'
                    removalReason = 'saturation';
                    switch sub_nm
                        case {'012','032','039','047','055','073','095'}
                            error(['Subject ',sub_nm,' should not be included']);
                        case {'001','002','003','004','005','009',...
                                '013','015','017','018',...
                                '020','021','022','027',...
                                '035','036','038',...
                                '040','042','043','044','045','046','048',...
                                '052','053','054','056','058','059',...
                                '060','062','065','069',...
                                '072','074','076','078',...
                                '080','081','082','083','086','087','088',...
                                '093','094','097','099','100'}
                            switch sub_nm
                                case '001'
                                    run_nm_toRemove = {'2_010_run4_20220429',...
                                        '2_008_run2_20220429'};
                                case '002'
                                    run_nm_toRemove = {'2_011_run4_20220222',...
                                        '2_009_run3_20220222',...
                                        '2_008_run2_20220222'};
                                case '003'
                                    run_nm_toRemove = {'4_011_run4_20220406',...
                                        '4_009_run2_20220406'};
                                case '004'
                                    run_nm_toRemove = {'2_010_run4_20220902',...
                                        '2_009_run3_20220902',...
                                        '2_007_run1_20220902'};
                                case '005'
                                    run_nm_toRemove = {'2_010_run4_20220506',...
                                        '2_008_run2_20220506'};
                                case '009'
                                    run_nm_toRemove = {'2_010_run4_20220217'};
                                case '013'
                                    run_nm_toRemove = {'2_010_run4_20220707',...
                                        '2_008_run2_20220707'};
                                case '015'
                                    run_nm_toRemove = {'2_010_run4_20220412',...
                                        '2_008_run2_20220412',...
                                        '2_007_run1_20220412'};
                                case '017'
                                    run_nm_toRemove = {'2_010_run4_20220104'};
                                case '018'
                                    run_nm_toRemove = {'2_010_run4_20220426',...
                                        '2_008_run2_20220426'};
                                case '020'
                                    run_nm_toRemove = {'2_009_run3_20220126'};
                                case '021'
                                    run_nm_toRemove = {'2_010_run4_20220519',...
                                        '2_009_run3_20220519',...
                                        '2_008_run2_20220519'};
                                case '022'
                                    run_nm_toRemove = {'2_010_run4_20220303'};
                                case '027'
                                    run_nm_toRemove = {'2_010_run4_20220624',...
                                        '2_009_run3_20220624',...
                                        '2_008_run2_20220624'};
                                case '035'
                                    run_nm_toRemove = {'5_010_run4_20220413',...
                                        '5_009_run3_20220413'};
                                case '036'
                                    run_nm_toRemove = {'2_010_run4_20211110',...
                                        '2_008_run2_20211110',...
                                        '2_007_run1_20211110'};
                                case '038'
                                    run_nm_toRemove = {'3_010_run4_20220701',...
                                        '3_008_run2_20220701'};
                                case '040'
                                    run_nm_toRemove = {'2_007_run1_20220304'};
                                case '042'
                                    run_nm_toRemove = {'2_010_run4_20220527',...
                                        '2_009_run3_20220527',...
                                        '2_007_run1_20220527'};
                                case '043'
                                    run_nm_toRemove = {'2_010_run4_20220531',...
                                        '2_008_run2_20220531'};
                                case '044'
                                    run_nm_toRemove = {'3_010_run4_20220408',...
                                        '3_009_run3_20220408'};
                                case '045'
                                    run_nm_toRemove = {'3_010_run4_20220118',...
                                        '3_009_run3_20220118'};
                                case '046'
                                    run_nm_toRemove = {'2_010_run4_20220128',...
                                        '2_008_run2_20220128'};
                                case '047'
                                    run_nm_toRemove = {'2_010_run4_20220311',...
                                        '2_008_run2_20220311',...
                                        '2_007_run1_20220311'};
                                case '048'
                                    run_nm_toRemove = {'3_008_run2_20220518'};
                                case '052'
                                    run_nm_toRemove = {'2_010_run3_20220216',...
                                        '2_008_run1_20220216'};
                                case '053'
                                    run_nm_toRemove = {'2_010_run4_20220727',...
                                        '2_008_run2_20220727'};
                                case '054'
                                    run_nm_toRemove = {'2_007_run1_20220208'};
                                case '056'
                                    run_nm_toRemove = {'2_010_run4_20220218',...
                                        '2_009_run3_20220218',...
                                        '2_007_run1_20220218'};
                                case '058'
                                    run_nm_toRemove = {'3_010_run4_20220621',...
                                        '3_008_run2_20220621'};
                                case '059'
                                    run_nm_toRemove = {'3_009_run3_20220628',...
                                        '3_007_run1_20220628'};
                                case '060'
                                    run_nm_toRemove = {'3_009_run3_20220125'};
                                case '062'
                                    run_nm_toRemove = {'2_009_run3_20220722'};
                                case '065'
                                    run_nm_toRemove = {'2_010_run4_20220119',...
                                        '2_009_run3_20220119',...
                                        '2_008_run2_20220119'};
                                case '069'
                                    run_nm_toRemove = {'2_010_run4_20220414',...
                                        '2_009_run3_20220414',...
                                        '2_008_run2_20220414'};
                                case '072'
                                    run_nm_toRemove = {'2_011_run4_20220601',...
                                        '2_008_run2_20220601'};
                                case '074'
                                    run_nm_toRemove = {'3_009_run3_20211102'};
                                case '076'
                                    run_nm_toRemove = {'2_010_run4_20220420',...
                                        '2_008_run2_20220420'};
                                case '078'
                                    run_nm_toRemove = {'2_010_run4_20220630',...
                                        '2_009_run3_20220630'};
                                case '080'
                                    run_nm_toRemove = {'2_010_run4_20220629',...
                                        '2_008_run2_20220629'};
                                case '081'
                                    run_nm_toRemove = {'2_009_run3_20220107',...
                                        '2_008_run2_20220107',...
                                        '2_007_run1_20220107'};
                                case '082'
                                    run_nm_toRemove = {'2_010_run4_20220308',...
                                        '2_009_run3_20220308'};
                                case '083'
                                    run_nm_toRemove = {'2_010_run4_20220520',...
                                        '2_008_run2_20220520',...
                                        '2_007_run1_20220520'};
                                case '086'
                                    run_nm_toRemove = {'2_010_run4_20221006',...
                                        '2_009_run3_20221006',...
                                        '2_008_run2_20221006'};
                                case '087'
                                    run_nm_toRemove = {'2_010_run4_20211215',...
                                        '2_009_run3_20211215'};
                                case '088'
                                    run_nm_toRemove = {'2_010_run4_20220512'};
                                case '094'
                                    run_nm_toRemove = {'2_010_run4_20220914'};
                                case '097'
                                    run_nm_toRemove = {'2_009_run3_20220721',...
                                        '2_007_run1_20220721'};
                                case '099'
                                    run_nm_toRemove = {'1_010_run4_20220705',...
                                        '1_009_run3_20220705',...
                                        '1_007_run1_20220705'};
                                case '100'
                                    run_nm_toRemove = {'2_010_run4_20220324',...
                                        '2_009_run3_20220324',...
                                        '2_008_run2_20220324'};
                                otherwise
                                    error([sub_nm,' is missing.'])
                            end % subject filter
                            [subj_scan_folders_names] = badRunsClearing(run_nm_toRemove,...
                                subj_scan_folders_names, sub_nm, removalReason);
                    end % subject group
                    
                    %% filter runs with too much movement (lenient)
                case 'fMRI_noMove_bis'
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
                                    run_nm_toRemove = {'2_010_run4_20220519',...
                                        '2_009_run3_20220519',...
                                        '2_008_run2_20220519'};
                                case '029'
                                    run_nm_toRemove = {'2_010_run4_20220310'};
                                case '044'
                                    run_nm_toRemove = {'3_010_run4_20220408',...
                                        '3_008_run2_20220408'};
                                case '047'
                                    run_nm_toRemove = {'2_009_run3_20220311'};
                                case '053'
                                    run_nm_toRemove = {'2_009_run3_20220727'};
                                case '054'
                                    run_nm_toRemove = {'2_010_run4_20220208',...
                                        '2_008_run2_20220208'};
                                case '058'
                                    run_nm_toRemove = {'3_010_run4_20220621',...
                                        '3_008_run2_20220621'};
                                case '062'
                                    run_nm_toRemove = {'2_010_run4_20220722',...
                                        '2_009_run3_20220722',...
                                        '2_008_run2_20220722'};
                                case '071'
                                    run_nm_toRemove = {'2_010_run4_20220405',...
                                        '2_008_run2_20220405'};
                                case '076'
                                    run_nm_toRemove = {'2_007_run1_20220420'};
                                case '078'
                                    run_nm_toRemove = {'2_010_run4_20220630',...
                                        '2_007_run1_20220630'};
                                case '080'
                                    run_nm_toRemove = {'2_009_run3_20220629'};
                                case '083'
                                    run_nm_toRemove = {'2_009_run3_20220520'};
                                case '087'
                                    run_nm_toRemove = {'2_010_run4_20211215',...
                                        '2_009_run3_20211215',...
                                        '2_008_run2_20211215'};
                                case '097'
                                    run_nm_toRemove = {'2_010_run4_20220721',...
                                        '2_009_run3_20220721',...
                                        '2_008_run2_20220721'};
                                case '099'
                                    run_nm_toRemove = {'1_010_run4_20220705',...
                                        '1_008_run2_20220705'};
                                otherwise
                                    error([sub_nm,' is missing.'])
                            end % subject
                            [subj_scan_folders_names] = badRunsClearing(run_nm_toRemove,...
                                subj_scan_folders_names, sub_nm, removalReason);
                    end % subject group
                    %% filter too much movement (stringent)
                case 'fMRI_noMove_ter'
                    removalReason = 'too much movement';
                    switch sub_nm
                        case {'005','008','012','018',...
                                '021','022','024','029',...
                                '043','044','047',...
                                '050','052','053','054','056','058',...
                                '062','064','065',...
                                '071','076','078','079',...
                                '080','083','086','087',...
                                '090','093','094','095', '097','099'}
                            switch sub_nm
                                case {'008','022','024'} % too much movement ALL runs
                                    error(['subject ',sub_nm,'should not be included at all']);
                                case '005'
                                    run_nm_toRemove = {'2_010_run4_20220506',...
                                        '2_009_run3_20220506',...
                                        '2_007_run1_20220506'};
                                case '012'
                                    run_nm_toRemove = {'2_008_run2_20221021'};
                                case '018'
                                    run_nm_toRemove = {'2_009_run3_20220426'};
                                case '021'
                                    run_nm_toRemove = {'2_010_run4_20220519',...
                                        '2_009_run3_20220519',...
                                        '2_008_run2_20220519'};
                                case '029'
                                    run_nm_toRemove = {'2_010_run4_20220310'};
                                case '043'
                                    run_nm_toRemove = {'2_010_run4_20220531',...
                                        '2_008_run2_20220531',...
                                        '2_007_run1_20220531'};
                                case '044'
                                    run_nm_toRemove = {'3_010_run4_20220408',...
                                        '3_008_run2_20220408'};
                                case '047'
                                    run_nm_toRemove = {'2_009_run3_20220311'};
                                case '050'
                                    run_nm_toRemove = {'3_013_run4_20220401'};
                                case '052'
                                    run_nm_toRemove = {'2_011_run4_20220216'};
                                case '053'
                                    run_nm_toRemove = {'2_009_run3_20220727'};
                                case '054'
                                    run_nm_toRemove = {'2_010_run4_20220208',...
                                        '2_008_run2_20220208'};
                                case '056'
                                    run_nm_toRemove = {'2_009_run3_20220218',...
                                        '2_008_run2_20220218'};
                                case '058'
                                    run_nm_toRemove = {'3_010_run4_20220621',...
                                        '3_008_run2_20220621'};
                                case '062'
                                    run_nm_toRemove = {'2_010_run4_20220722',...
                                        '2_009_run3_20220722',...
                                        '2_008_run2_20220722'};
                                case '064'
                                    run_nm_toRemove = {'2_008_run2_20211123'};
                                case '065'
                                    run_nm_toRemove = {'2_009_run3_20220119'};
                                case '069'
                                    run_nm_toRemove = {'2_010_run4_20220414',...
                                        '2_009_run3_20220414'};
                                case '071'
                                    run_nm_toRemove = {'2_010_run4_20220405',...
                                        '2_008_run2_20220405'};
                                case '076'
                                    run_nm_toRemove = {'2_007_run1_20220420'};
                                case '078'
                                    run_nm_toRemove = {'2_010_run4_20220630',...
                                        '2_007_run1_20220630'};
                                case '079'
                                    run_nm_toRemove = {'2_009_run3_20211201'};
                                case '080'
                                    run_nm_toRemove = {'2_009_run3_20220629'};
                                case '083'
                                    run_nm_toRemove = {'2_009_run3_20220520'};
                                case '086'
                                    run_nm_toRemove = {'2_009_run3_20221006'};
                                case '087'
                                    run_nm_toRemove = {'2_010_run4_20211215',...
                                        '2_009_run3_20211215',...
                                        '2_008_run2_20211215'};
                                case '090'
                                    run_nm_toRemove = {'2_009_run3_20211130'};
                                case '093'
                                    run_nm_toRemove = {'2_009_run3_20220309',...
                                        '2_007_run1_20220309'};
                                case '094'
                                    run_nm_toRemove = {'2_009_run3_20220914'};
                                case '095'
                                    run_nm_toRemove = {'2_007_run1_20211119'};
                                case '097'
                                    run_nm_toRemove = {'2_010_run4_20220721',...
                                        '2_009_run3_20220721',...
                                        '2_008_run2_20220721'};
                                case '099'
                                    run_nm_toRemove = {'1_010_run4_20220705',...
                                        '1_008_run2_20220705'};
                                otherwise
                                    error([sub_nm,' is missing.'])
                            end % individual subject
                            [subj_scan_folders_names] = badRunsClearing(run_nm_toRemove,...
                                subj_scan_folders_names, sub_nm, removalReason);
                    end % subject group
                    %% remove saturated tasks AND movement
                case 'fMRI_noSatTask_noMove_bis'
                    removalReason = 'too much movement or saturated task';
                    switch sub_nm
                        case {'008','021','022','024','029',...
                                '044','047',...
                                '053','054','058','062','069',...
                                '071','076','078','080','083','087',...
                                '097','099'}
                            switch sub_nm
                                case {'008','022','024'} % too much movement ALL runs
                                    error(['subject ',sub_nm,'should not be included at all']);
                                case '021'
                                    run_nm_toRemove = {'2_010_run4_20220519',...
                                        '2_009_run3_20220519',...
                                        '2_008_run2_20220519'};
                                case '027'
                                    run_nm_toRemove = {'2_010_run4_20220624',...
                                        '2_008_run2_20220624'};
                                case '029'
                                    run_nm_toRemove = {'2_010_run4_20220310'};
                                case '044'
                                    run_nm_toRemove = {'3_010_run4_20220408',...
                                        '3_008_run2_20220408'};
                                case '047'
                                    error(['subject ',sub_nm,' should not be included (saturated ALL tasks)']);
                                case '052'
                                    run_nm_toRemove = {'2_010_run3_20220216',...
                                        '2_008_run1_20220216'};
                                case '053'
                                    run_nm_toRemove = {'2_009_run3_20220727'};
                                case '054'
                                    run_nm_toRemove = {'2_010_run4_20220208',...
                                        '2_008_run2_20220208'};
                                case '058'
                                    run_nm_toRemove = {'3_010_run4_20220621',...
                                        '3_008_run2_20220621'};
                                case '062'
                                    run_nm_toRemove = {'2_010_run4_20220722',...
                                        '2_009_run3_20220722',...
                                        '2_008_run2_20220722'};
                                case '069'
                                    run_nm_toRemove = {'2_010_run4_20220414',...
                                        '2_008_run2_20220414'};
                                case '071'
                                    run_nm_toRemove = {'2_010_run4_20220405',...
                                        '2_008_run2_20220405'};
                                case '076'
                                    run_nm_toRemove = {'2_010_run4_20220420',...
                                        '2_008_run2_20220420'};
                                case '078'
                                    run_nm_toRemove = {'2_010_run4_20220630',...
                                        '2_007_run1_20220630'};
                                case '080'
                                    run_nm_toRemove = {'2_009_run3_20220629'};
                                case '083'
                                    run_nm_toRemove = {'2_009_run3_20220520'};
                                case '087'
                                    run_nm_toRemove = {'2_010_run4_20211215',...
                                        '2_009_run3_20211215',...
                                        '2_008_run2_20211215'};
                                case '095'
                                    run_nm_toRemove = {'2_010_run4_20211119',...
                                        '2_008_run2_20211119'};
                                case '097'
                                    run_nm_toRemove = {'2_010_run4_20220721',...
                                        '2_009_run3_20220721',...
                                        '2_008_run2_20220721'};
                                case '099'
                                    run_nm_toRemove = {'1_010_run4_20220705',...
                                        '1_008_run2_20220705'};
                                otherwise
                                    error([sub_nm,' is missing.'])
                            end % subject
                            [subj_scan_folders_names] = badRunsClearing(run_nm_toRemove,...
                                subj_scan_folders_names, sub_nm, removalReason);
                    end % subject group
                    %% remove saturated runs AND movement
                case 'fMRI_noSatRun_noMove_bis'
                    removalReason = 'too much movement or saturated run';
                    switch sub_nm
                        case {'002','004','005','008',...
                                '012',...
                                '021','022','024','027','029',...
                                '032','038',...
                                '044','047','048',...
                                '052','053','054','055','058',...
                                '061','062','069',...
                                '071','076','078',...
                                '080','081','082','083','087','088',...
                                '095','097','099','100'}
                            switch sub_nm
                                case {'008','022','024','047','097'} % bad in ALL runs
                                    error(['subject ',sub_nm,'should not be included at all']);
                                case '002'
                                    run_nm_toRemove = {'2_009_run3_20220222'};
                                case '004'
                                    run_nm_toRemove = {'2_009_run3_20220902'};
                                case '005'
                                    run_nm_toRemove = {'2_010_run4_20220506'};
                                case '012'
                                    run_nm_toRemove = {'2_010_run4_20221021'};
                                case '021'
                                    run_nm_toRemove = {'2_010_run4_20220519',...
                                        '2_009_run3_20220519',...
                                        '2_008_run2_20220519'};
                                case '027'
                                    run_nm_toRemove = {'2_010_run4_20220624',...
                                        '2_009_run3_20220624',...
                                        '2_008_run2_20220624'};
                                case '029'
                                    run_nm_toRemove = {'2_010_run4_20220310'};
                                case '032'
                                    run_nm_toRemove = {'1_009_run3_20220428'};
                                case '038'
                                    run_nm_toRemove = {'3_002_run2_20220701'};
                                case '044'
                                    run_nm_toRemove = {'3_010_run4_20220408',...
                                        '3_009_run3_20220408',...
                                        '3_008_run2_20220408'};
                                case '048'
                                    run_nm_toRemove = {'3_008_run2_20220518'};
                                case '052'
                                    run_nm_toRemove = {'2_010_run3_20220216',...
                                        '2_008_run1_20220216'};
                                case '053'
                                    run_nm_toRemove = {'2_009_run3_20220727'};
                                case '054'
                                    run_nm_toRemove = {'2_010_run4_20220208',...
                                        '2_008_run2_20220208'
                                        '2_007_run1_20220208'};
                                case '055'
                                    run_nm_toRemove = {'2_010_run4_20220202'};
                                case '058'
                                    run_nm_toRemove = {'3_010_run4_20220621',...
                                        '3_008_run2_20220621'};
                                case '061'
                                    run_nm_toRemove = {'2_007_run1_20211126'};
                                case '062'
                                    run_nm_toRemove = {'2_010_run4_20220722',...
                                        '2_009_run3_20220722',...
                                        '2_008_run2_20220722'};
                                case '069'
                                    run_nm_toRemove = {'2_010_run4_20220414',...
                                        '2_008_run2_20220414'};
                                case '071'
                                    run_nm_toRemove = {'2_010_run4_20220405',...
                                        '2_008_run2_20220405'};
                                case '076'
                                    run_nm_toRemove = {'2_010_run4_20220420',...
                                        '2_008_run2_20220420',...
                                        '2_007_run1_20220420'};
                                case '078'
                                    run_nm_toRemove = {'2_010_run4_20220630',...
                                        '2_007_run1_20220630'};
                                case '080'
                                    run_nm_toRemove = {'2_009_run3_20220629'};
                                case '081'
                                    run_nm_toRemove = {'2_009_run3_20220107'};
                                case '082'
                                    run_nm_toRemove = {'2_010_run4_20220308',...
                                        '2_009_run3_20220308'};
                                case '083'
                                    run_nm_toRemove = {'2_010_run4_20220520',...
                                        '2_009_run3_20220520'};
                                case '087'
                                    run_nm_toRemove = {'2_010_run4_20211215',...
                                        '2_009_run3_20211215',...
                                        '2_008_run2_20211215'};
                                case '088'
                                    run_nm_toRemove = {'2_010_run4_20220512'};
                                case '095'
                                    run_nm_toRemove = {'2_010_run4_20211119',...
                                        '2_009_run3_20211119',...
                                        '2_008_run2_20211119'};
                                case '099'
                                    run_nm_toRemove = {'1_010_run4_20220705',...
                                        '1_009_run3_20220705',...
                                        '1_008_run2_20220705'};
                                case '100'
                                    run_nm_toRemove = {'2_010_run4_20220324',...
                                        '2_009_run3_20220324'};
                                    otherwise
                                    error([sub_nm,' is missing.'])
                            end % subject
                            [subj_scan_folders_names] = badRunsClearing(run_nm_toRemove,...
                                subj_scan_folders_names, sub_nm, removalReason);
                    end % subject group
            end % condition
    end % study
else
    subj_scan_folders_names = [];
end % subj_scan_folders_names correction

%% fix index of the run
if exist('iRun','var') && ~isempty(iRun)
    jRun = runs.runsToKeep(iRun);
end

end % function