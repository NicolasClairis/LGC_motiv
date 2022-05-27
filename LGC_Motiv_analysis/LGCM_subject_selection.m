function[subject_id, NS] = LGCM_subject_selection(study_nm, condition)
% [subject_id, NS] = LGCM_subject_selection(study_nm, condition)
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
% 'behavior_noSat': behavior but removing all runs (or subjects with all
% runs) with saturation
% 'fMRI': all fMRI compatible data
% 'fMRI_no_move': remove runs with too much movement
%
% OUTPUTS
% subject_id: list of subject names
%
% NS: number of subjects

switch study_nm
    case 'fMRI_pilots'
        subject_id = {'pilot_s1','pilot_s2','pilot_s3'};
        %         subject_id = {'pilot_s1','pilot_s2','pilot_s3'};
    case 'study1'
        % full list of all subjects included in the study
        fullSubList = {'001','002','003','005','008','009',...
            '015','017','018',...
            '020','021','022','029',...
            '032','035','036','039',...
            '040','042','044','045','046','047','048',...
            '050','052','054','055','056',...
            '060','061','064','065','068','069',...
            '071','074','075','076','079',...
            '081','082','083','087','088',...
            '090','093','095','100'};
        all_subs = fullSubList;
        % remove some subjects depending on the condition entered as input
        switch condition
            case {'behavior','fMRI'} % all subjects
%                 % for confidence remove saturated subjects:
% also remove 042 because data not added yet.
                bad_subs = ismember(fullSubList,{'042','047','095'});
                
            case {'behavior_noSatRun','fMRI_noSatRun'}
                % remove subjects who saturated the behavioral task in any
                % run
                bad_subs = ismember(fullSubList,{'002','005','032',...
                    '047','048','052','076','095','100'});
                % subjects with a full task saturated
                % 047: all ND for Em task (runs 2 and 4) and for Ep run 1
                % 052: all ND for Em task (runs 1 and 3)
                % 095: all ND for Ep task (runs 2 and 4)
                %
                % subjects with a full run saturated
                % 002: run 3 ND for Ep task
                % 005: run 4 ND for Em task
                % 032: run 3 ND for Em task
                % 048: run 2 ND for Em task
                % 076: run 4 ND for Em task
                % 100: run 3 (Em) and run 4 (Ep) ND
                 warning('check last subs');
            case {'behavior_noSatTask','fMRI_noSatTask'}
                % remove subjects for which either mental (Em) or physical
                % (Ep) task was fully saturated during choices
                bad_subs = ismember(fullSubList,{'047','052','095'});
                % 047: all ND for Em task (runs 2 and 4) and for Ep run 1
                % 052: all ND for Em task (runs 1 and 3)
                % 095: all ND for Ep task (runs 2 and 4)
                warning('check last subs');
            case 'fMRI_no_move'
                % ignore subjects with too much movement in ALL runs (runs
                % with too much movement will be filtered for each subject)
                bad_subs = ismember(fullSubList,{'008','022'});
                % 008 and 022 completely removed (because all runs are bad
                % in terms of movement)
            case 'fMRI_no_move_bis'
                % ignore subjects with too much movement in ANY RUN (ie no
                % filtering of bad runs, but only of bad subjects who had
                % only bad runs)
                bad_subs = ismember(fullSubList,{'008',...
                    '021','022','029',...
                    '044','047','054',...
                    '071','076','083','087'});
                % too much movement for 008 (all runs)
                % 021 (run 2, 3 and 4)
                % 029 (run 4 physical),
                % 044 (run2 and 4 physical), 047 (run 3 mostly, physical run),
                % 054 (run2 and run4 mostly, ie physical runs),
                % 071 (run2 and run 4 ie physical runs), 076 (run 1)
                % 083 (run3 physical), 087(runs 2-4)
                %
                % borderline movement: 005(run1,3 and 4 a little bit of movement), 022 (all),
                % 040 (run 4), 050 (run4), 052 (run 2),
                % 056 (run 2 + 3), 064 (run 2), 065 (run3), 069 (run3 and 4),
                % 079 (run 3), 090 (run 3), 093 (run1 + 3), 095 (run 1)
        end
        % remove irrelevant subjects from the current analysis
        all_subs(bad_subs) = [];
        % reminder to check the last acquired subjects:
        warning('check 042,... and other last subjects for movement');
        
        %% restrict to subjects of interest
        subject_id = all_subs;
        %% Notes:
        % '064': creates bug when focusing on chosen option because too few
        % trials during run 3
        %
        % '054': very weird behavior, could be considered as outlier
    case 'study2'
        %         subject_id = {}; % 'XXX'
        error('experiment hasn''started yet...');
    otherwise
        error('error in study definition');
end
NS = length(subject_id);

end % function