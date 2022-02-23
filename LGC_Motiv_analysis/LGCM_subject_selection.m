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
        all_subs = {'002','008','009','017','020','036',...
            '039','045','046','052','054','055','056',...
            '060','061','064','065','074','079',...
            '081','087','090','095'};
        all_subs_no_mvmt = {'002','009','017','020','036',...
            '039','045','046','052','055','056',...
            '060','061','064','065','074','079',...
            '081','090','095'};


        %% restrict to subjects of interest
        subject_id = all_subs;
% subject_id={'002'};
% subject_id = {'074','036','064','061',...
%             '090','079','087','017','039','081','045','065',...
%             '060','020','046','055','008','009'};
% subject_id= {'074','036','061',...
%             '090','079','087','017','039','081','045','065',...
%             '060','020','046','055','008','009'};
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