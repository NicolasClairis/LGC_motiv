function[subject_id, NS] = LGCM_subject_selection(study_nm)
% [subject_id, NS] = LGCM_subject_selection(study_nm)
% LGCM_subject_selection will select the subject names and number of
% subjects for the current study.
%
% INPUTS
% study_nm: definition of the study on which you want to analyze the data
% 'fMRI_pilots': pilots
% 'study1': first study (dmPFC + AI)
% 'study2': second study (clinical trial)
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
        %% all subjects
        all_subs = {'002','003','008','009','015','017','020','022','029',...
            '035','036','039','040','044','045','046','047','050','052',...
            '054','055',...
            '056','060','061','064','065','068','069','071','074','075',...
            '076','079',...
            '081','082','087','090','093','095','100'};
        all_subs_no_mvmt = {'002','003','009','017','020','022','035','036',...
            '039','040','045','046','050','052','055','056',...
            '060','061','064','065','068','069','074','079',...
            '081','090','093','095','100'};
        % too much movement for 008, 029 (run 4 physical), 044 (run2 and 4 physical),
        % 047 (run 3 mostly, physical run), 054 (run2 and run4 mostly, ie physical runs),
        % 071 (run2 and run 4 ie physical runs), 076 (run 1), 087
        %
        % borderline movement: 022, 040 (run 4), 050 (run4), 052 (run 2),
        % 056 (run 2 + 3), 064 (run 2), 065 (run3), 069 (run3 and 4),
        % 079 (run 3), 090 (run 3), 093 (run1 + 3), 095 (run 1)
        warning('check 015, 035, 044, 068, 069, 076 and... for movement');


        %% restrict to subjects of interest
        subject_id = all_subs;
        %% Notes:
        % bad because of too much movement:
        % '008', '054', '087'
        %
        % '095': always chose non-default in physical task => hard to
        % analyze
        %
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