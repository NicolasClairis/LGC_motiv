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
% 'behavior': behavioral files
% 'behavior_noSatRun': behavior but removing all subjects who saturated 
% completely any of the runs
% 'behavior_noSatTask': remove subjects who saturated completely one of the
% tasks
% 'fMRI': all fMRI compatible data
% 'fMRI_no_move': remove subjects with too much movement in ALL runs
% 'fMRI_no_move_bis': remove subjects with too much movement in ANY run
% 'fMRI_noSatTask': remove subjects who saturated completely one of the
% tasks
% 'fMRI_noSatRun': remove subjects who saturated completely any of the runs
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
        % full list of all subjects included in the study
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
        % firstly remove subjects where behavior and fMRI could not be performed:
        bad_subs1 = ismember(fullSubList,{'030','049'});
        fullSubList(bad_subs1) = [];
        % initialize the list of subjects to consider
        all_subs = fullSubList;
        
        % remove some subjects depending on the condition entered as input
        switch condition
            case {'behavior','fMRI'} % all subjects
                % for confidence, you should remove saturated subjects
                bad_subs = false(1,length(fullSubList));
                warning(['if you want to look at confidence in your GLM, ',...
                    'you should remove the subjects saturating behavior using ',...
                    'behavior_noSatRun/fMRI_noSatRun/behavior_noSatTask/fMRI_noSatTask',...
                    ' conditions to filter.']);
            case {'behavior_noSatRun','fMRI_noSatRun'}
                % remove subjects who saturated the behavioral task in any
                % run
                bad_subs = ismember(fullSubList,{'002','005','012',...
                    '027','032',...
                    '047','048','052',...
                    '076','095','100'});
                % subjects with a full task saturated
                % 047: all ND for Em task (runs 2 and 4) and for Ep run 1
                % 052: all ND for Em task (runs 1 and 3)
                % 095: all ND for Ep task (runs 2 and 4)
                %
                % subjects with a full run saturated
                % 002: run 3 ND for Ep task
                % 005: run 4 ND for Em task
                % 012: run 4 ND for Em task
                % 032: run 3 ND for Em task
                % 048: run 2 ND for Em task
                % 076: run 4 ND for Em task
                % 100: run 3 (Em) and run 4 (Ep) ND
            case {'behavior_noSatTask','fMRI_noSatTask'}
                % remove subjects for which either mental (Em) or physical
                % (Ep) task was fully saturated during choices
                bad_subs = ismember(fullSubList,{'027','047','052','095'});
                % 027: all ND for Em task (runs 2 and 4)
                % 047: all ND for Em task (runs 2 and 4) and for Ep run 1
                % 052: all ND for Em task (runs 1 and 3)
                % 095: all ND for Ep task (runs 2 and 4)
            case 'fMRI_no_move'
                % ignore subjects with too much movement in ALL runs (runs
                % with too much movement will be filtered for each subject)
                bad_subs = ismember(fullSubList,{'008','022','024'});
                % 008 and 022 completely removed (because all runs are bad
                % in terms of movement)
            case 'fMRI_no_move_bis'
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
                % 022 (all)
                % 040 (run 4)
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
                % 032, 035, 036, 038, 039, 042, 045, 046, 048, 055, 059, 060,
                % 061, 062, 068, 072, 073, 074, 075, 081, 082, 100
            case 'fMRI_GLM59'
                bad_subs = ismember(fullSubList,{'002','005',...
                    '013','017',...
                    '022','027',...
                    '032','038',...
                    '043','044','047','048',...
                    '052','054','055','058',...
                    '062','069',...
                    '072','074','076','078',...
                    '081','082','083','088',...
                    '095','097','100'});
        end
        % remove also based on gender
        males = {'002','004',...
            '013','017',...
            '020','021','022','028',...
            '032','034','035','036',...
            '040','043','045','047','049',...
            '053','054','056','058','059',...
            '060','065','068','069',...
            '073','074','076',...
            '085','086','087',...
            '090','091','092','094','097'};
        females = {'001','003','005','008','009',...
            '011','012','015','018','019',...
            '024','027','029',...
            '030','038','039',...
            '042','044','046','048',...
            '050','052','055',...
            '061','062','064','066',...
            '071','072','075','078','079',...
            '080','081','082','083','088',...
            '093','095','099','100'};
        switch genderFilter
            case 'males'
                bad_subs(ismember(fullSubList, females)) = true;
            case 'females'
                bad_subs(ismember(fullSubList, males)) = true;
        end
        % remove subjects who did behavior but not fMRI
        if strcmp(condition(1:4),'fMRI')
            bad_subs(ismember(fullSubList,{'034','091'})) = true;
        end
        % remove irrelevant subjects from the current analysis
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