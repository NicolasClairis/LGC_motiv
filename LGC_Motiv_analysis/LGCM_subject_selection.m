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
        all_subs = {'002','008','009','017','020','022','029',...
            '036','039','040','045','046','047','050','052','054','055',...
            '056','060','061','064','065','074','079',...
            '081','082','087','090','093','095','100'};
        all_subs_no_mvmt = {'002','009','017','020','022','036',...
            '039','040','045','046','052','055','056',...
            '060','061','064','065','074','079',...
            '081','090','093','095'};
        % too much movement for 008, 029 (run 4), 047 (run 3 mostly),
        % 054, 087
        % borderline movement: 022, 040 (run 4), 052 (run 2), 056 (run 2 + 3),
        % 064 (run 2), 079 (run 3), 090 (run 3), 093 (run1 + 3), 095 (run 1)
        warning('check 100 and... for movement');


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