function[subject_id, NS] = LGCM_subject_selection(study_nm, condition, genderFilter)
% [subject_id, NS] = LGCM_subject_selection(study_nm, condition, genderFilter)
% LGCM_subject_selection will select the subject names and number of
% subjects for the current study.
%
% INPUTS
% study_nm: definition of the study on which you want to analyze the data
% 'fMRI_pilots': pilots
% 'study1': first study (dmPFC + AI)
% 'study2': second study (clinical trial)
%
% condition:
% 'fullList': full list of subjects, including those who did not perform
% the behavioral task and fMRI
% 'behavior': all subjects (who did perform the behavioral task)
% 'behavior_noSatRunSub': all subjects but removing all subjects who saturated 
% completely any of the runs
% 'behavior_noSatTaskSub': all subjects but removing all subjects who saturated 
% completely any of the tasks (ie the two runs of a given task)
% 'behavior_noSatRun': like 'behavior' list but saturation runs will be
% removed from the subjects concerned
% 'behavior_noSatTask': like 'behavior' list but saturation tasks will be
% removed from the subjects concerned
% 'respiration_and_noSatRun': respiration ok and no saturation run
% 'fMRI': all subjects who performed the behavioral task in the fMRI
% (removing only the runs where fMRI crashed)
% 'fMRI_noMoveSub': remove subjects with too much movement in ALL runs
% 'fMRI_noMoveSub_bis': remove subjects with too much movement in ANY run
% 'fMRI_noMoveSub_ter': remove subjects with too much movement in ANY run,
% even including borderline movement
% 'fMRI_noSatTaskSub': remove subjects who saturated completely one of the
% tasks
% 'fMRI_noSatRunSub': remove subjects who saturated completely any of the
% runs
% 'fMRI_noSatTaskSub_noMove_bis_Sub': remove subjects who saturated 
% completely any of the runs or who moved too much in any run
% 'fMRI_noSatRun': like 'fMRI' list but saturation runs will be
% removed from the subjects concerned
% 'fMRI_noSatTask': like 'fMRI' list but saturation tasks will be
% removed from the subjects concerned
% 'fMRI_noMove_bis': like 'fMRI' but runs with too much movement will be
% removed
% 'fMRI_noMove_ter': like 'fMRI' but runs with too much movement will be
% removed (even including borderline movement)
% 'fMRI_noSatTask_noMove_bis': like 'fMRI' but any run with saturation or
% too much movement will be removed
% 'behavior_noSatRun_bayesianMdl': like 'behavior_noSatRun' but removing 
% subjects who saturated because bayesian model was not applied on them
% 'fMRI_noSatRun_bayesianMdl': like 'fMRI_noSatRun' but removing subjects
% who saturated because bayesian model was not applied on them
% 'behavior_noSatTask_bayesianMdl': like 'behavior_noSatTask' but removing 
% subjects who saturated because bayesian model was not applied on them
% 'fMRI_noSatTask_bayesianMdl': like 'fMRI_noSatTask' but removing subjects
% who saturated because bayesian model was not applied on them
% 'fMRI_noSatTask_noMove_bis_bayesianMdl': like 'fMRI_noSatTask_noMove_bis' 
% but removing subjects who saturated because bayesian model was not 
% applied on them
%
% genderFilter:
% 'all': by default, include all subjects
% 'males': remove females
% 'females': remove males
%
% OUTPUTS
% subject_id: list of subject names
%
% NS: number of subjects included in the final list

% by default include all subjects
if ~exist('genderFilter','var') || isempty(genderFilter)
    genderFilter = 'all';
end

% extract list of subjects
switch study_nm
    case 'fMRI_pilots'
        subject_id = {'pilot_s1','pilot_s2','pilot_s3'};
        %         subject_id = {'pilot_s1','pilot_s2','pilot_s3'};
    case 'study1'
        %% full list of all subjects included in the study
        fullSubList = {'001','002','003','004','005','008','009',...
            '011','012','013','015','017','018','019',...
            '020','021','022','024','027','029',...
            '030','032','034','035','036','038','039',...
            '040','042','043','044','045','046','047','048','049',...
            '050','052','053','054','055','056','058','059',...
            '060','061','062','064','065','068','069',...
            '071','072','073','074','075','076','078','079',...
            '080','081','082','083','085','086','087','088',...
            '090','091','093','094','095','097','099','100'};
        %% firstly remove subjects where behavior and fMRI could not be performed:
        switch condition
            case 'fullList'
                bad_subs = false(1,length(fullSubList));
            otherwise
                bad_subs1 = ismember(fullSubList,{'030','049'});
                fullSubList(bad_subs1) = [];
        end
        %% initialize the list of subjects to consider
        all_subs = fullSubList;
        
        %% remove some subjects depending on the condition entered as input
        switch condition
            case {'behavior','fMRI',...
                    'fMRI_noMove_bis','fMRI_noMove_ter',...
                    'fMRI_noSatRun_choiceSplit_Elvl'} % all subjects
                % (but removing the bad runs if the condition requires it)
                bad_subs = false(1,length(fullSubList));
            case {'behavior_noSatRun','behavior_noSatTask',...
                    'respiration_and_noSatRun',...
                    'fMRI_noSatRun','fMRI_noSatTask','fMRI_noSatTask_noMove_bis'}
                % removing subjects where no run survives after removing
                % runs with too much saturation
                bad_subs = ismember(fullSubList,{'047'});
            case {'behavior_noSatRun_bayesianMdl','fMRI_noSatRun_bayesianMdl',...
                    'behavior_noSatTask_bayesianMdl','fMRI_noSatTask_bayesianMdl',...
                    'fMRI_noSatTask_noMove_bis_bayesianMdl'}
                % removing subjects where Arthur did not apply the
                % computational model (027, 052, 095) because one task was
                % saturated
                bad_subs = ismember(fullSubList,{'027','047','052','095'});
            case 'fMRI_noSatRun_choiceSplit_Elvl_bis'
                % removing subjects where no run survives after removing
                % runs with too much saturation
                bad_subs = ismember(fullSubList,{'012','032','039','047',...
                    '055','073','095'});
            case 'fMRI_noSatRun_noMove_bis'
                % removing subjects where no run survives after removing
                % runs with too much movement or saturation
                bad_subs = ismember(fullSubList,{'047','097'});
            case {'behavior_noSatRunSub','fMRI_noSatRunSub'}
                % remove subjects who saturated the behavioral task in any
                % run
                bad_subs = ismember(fullSubList,{'002','004','005',...
                    '012',...
                    '022','027','032','038',...
                    '044','047','048',...
                    '052','054','055','058',...
                    '061','062','069',...
                    '076','081','082','083','088',...
                    '095','097','099','100'});
                % subjects with a full task saturated
                % 027: all ND for Em task (runs 2 and 4) and all D run 3 Ep task
                % 047: all ND for all runs
                % 052: all ND for Em task (runs 1 and 3)
                % 069: almost all ND for Em task (runs 2 and 4)
                % 076: almost all ND for Em task (runs 2 and 4)
                % 095: all ND for Ep task (runs 2 and 4) and run 3 (Em)
                %
                % subjects with a full run saturated
                % 002: run 3 ND for Ep task
                % 004: run 3 D for Ep task
                % 005: run 4 ND for Em task
                % 012: run 4 ND for Em task
                % 022: run 3 ND for Em task
                % 032: run 3 ND for Em task
                % 038: run 2 ND for Em task
                % 044: run 3 ND for Em task
                % 048: run 2 ND for Em task
                % 054: run 1 ND for Em task
                % 055: run 4 ND for Ep task
                % 058: run 4 D for Ep task
                % 061: run 1 ND for Em task
                % 062: run 3 ND for Em task
                % 081: run 3 ND for Em task
                % 082: run 3 (Em) and run 4 (Ep) ND
                % 083: run 4 ND for Em task
                % 097: run 1 ND for Ep task
                % 099: run 3 ND for Em task
                % 100: run 3 (Em) and run 4 (Ep) ND
            case {'behavior_noSatTaskSub','fMRI_noSatTaskSub'}
                % remove subjects for which either mental (Em) or physical
                % (Ep) task was fully saturated during choices and remove
                % runs that were saturating if only one saturated
                bad_subs = ismember(fullSubList,{'027','047',...
                    '052','069','076','095'});
                % 027: all ND for Em task (runs 2 and 4)
                % 047: all ND for all tasks
                % 052: all ND for Em task (runs 1 and 3)
                % 069: almost all ND for Em task (runs 2 and 4)
                % 076: almost all ND for Em task (runs 2 and 4)
                % 095: all ND for Ep task (runs 2 and 4)
            case 'fMRI_noMoveSub'
                % ignore subjects with too much movement in ALL runs (runs
                % with too much movement will be filtered for each subject)
                bad_subs = ismember(fullSubList,{'008','022','024'});
                % 008, 022 and 024 completely removed (because all runs are bad
                % in terms of movement)
                % for the other subjects, the bad runs should be removed by
                % First_level_subRunFilter.m
            case 'fMRI_noMoveSub_bis'
                % ignore subjects with too much movement in ANY RUN (ie no
                % filtering of bad runs, but only of bad subjects who had
                % only bad runs)
                bad_subs = ismember(fullSubList,{'008',...
                    '021','022','024','029',...
                    '044','047',...
                    '053','054','058',...
                    '062',...
                    '071','076','078',...
                    '080','083','087',...
                    '097','099'});
                % too much movement for 008 (all runs)
                % 021 (run 2, 3 and 4)
                % 022 (all runs)
                % 024 (all runs with at least some movement)
                % 029 (run 4 physical),
                % 044 (run 2 and 4 physical)
                % 047 (run 3 mostly, physical run),
                % 053 (run3 mostly, physical run),
                % 054 (run 2 and run 4 mostly, ie physical runs),
                % 058 (run 2 and run 4, ie physical runs)
                % 062 (run 2 + 3 + 4, but mostly 4 physical run)
                % 071 (run 2 and run 4 ie physical runs)
                % 076 (run 1)
                % 078 (run 1, run 4 and also a bit run 2)
                % 080 (run 3 physical and a bit run 4 mental)
                % 083 (run 3 physical)
                % 087(runs 2-4)
                % 097 (runs 2, 3 and 4 mostly)
                % 099 (runs 2 and 4, physical runs)
                
                % borderline movement:
                % 005(run 1, 3 and 4 a little bit of movement)
                % 012 (run 2)
                % 018 (run 3)
                % 040 (run 4 but runs 3 and 4 should be ignored anyway)
                % 043 (run 2 and run 4)
                % 050 (run 4)
                % 052 (run 2),
                % 056 (run 2 + 3)
                % 064 (run 2)
                % 065 (run 3)
                % 069 (run 3 and 4)
                % 079 (run 3)
                % 086 (run 3)
                % 090 (run 3)
                % 093 (run 1 + 3)
                % 094 (run 3)
                % 095 (run 1)
                
                % no movement:
                % 001, 002, 003, 004, 009, 011, 013, 015, 017, 019, 020, 
                % 032, 035, 036, 038, 039, 042, 045, 046, 048, 055, 059,
                % 060, 061, 062, 068, 072, 073, 074, 075, 081, 082, 100
            case 'fMRI_noMoveSub_ter'
                % ignore subjects with too much movement in ANY RUN even if
                % the movement was really scarce
                bad_subs = ismember(fullSubList,{'005','008',...
                    '012','018',...
                    '021','022','024','029',...
                    '040','043','044','047',...
                    '050','052','053','054','056','058',...
                    '062','064','065','069',...
                    '071','076','078','079',...
                    '080','083','086','087',...
                    '090','093','094','095','097','099'});
            case 'fMRI_noSatTaskSub_noMove_bis_Sub' % remove subjects who 
                % either saturated a full task or moved too much
                bad_subs = ismember(fullSubList,{'008','022','024',...
                    '027','047','052','069','076','095'});
                % '008','022','024': too much movement in all runs
                % '027','047','052','069','076','095': one or two tasks fully saturated
        end
        %% split subjects based on gender
        males = {'002','004',...
            '013','017',...
            '020','021','022',...
            '032','034','035','036',...
            '040','043','045','047','049',...
            '053','054','056','058','059',...
            '060','065','068','069',...
            '073','074','076',...
            '085','086','087',...
            '090','091','094','097'};
        females = {'001','003','005','008','009',...
            '011','012','015','018','019',...
            '024','027','029',...
            '030','038','039',...
            '042','044','046','048',...
            '050','052','055',...
            '061','062','064',...
            '071','072','075','078','079',...
            '080','081','082','083','088',...
            '093','095','099','100'};
        % pilots male: 028, 092
        % pilots female: 066
        switch genderFilter
            case 'males'
                bad_subs(ismember(fullSubList, females)) = true;
            case 'females'
                bad_subs(ismember(fullSubList, males)) = true;
        end
        %% remove subjects who did behavior but not fMRI
        if strcmp(condition(1:4),'fMRI')
            bad_subs(ismember(fullSubList,{'034','091'})) = true;
        end
        %% remove irrelevant subjects from the current analysis
        all_subs(bad_subs) = [];

        %% restrict to subjects of interest
        subject_id = all_subs;
        %% Notes:
        % '064': creates bug when focusing on chosen option because too few
        % trials during run 3
        %
        % '054': very weird behavior, could be considered as outlier
    case 'study2'
        %         subject_id = {}; % 'XXX'
        error('experiment hasn''t started yet...');
    otherwise
        error('error in study definition');
end
NS = length(subject_id);

end % function