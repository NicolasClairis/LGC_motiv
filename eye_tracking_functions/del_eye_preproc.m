function [  ] = del_eye_preproc( which_study )
% del_eye_preproc( which_study ) deletes all the preprocessed data for the
% selected study. Allows to relaunch the preprocessing of the eye data
% while being sure that everything was cleaned up before launching the
% script.
%
% INPUTS
% which_study
% (0) Motiscan1 multiseq april/may 2016 study
% (1) Motiscan1 MBB june 2016 study
% (2) Motiscan2 march 2017 study
%
% See also preproc_eyeGazePupil_MotiScan.m & preproc_eyeGazePupil.m

%% select subject list depending on selected study
switch which_study
    case {0,1}
        [subject_id, NS] = MS1_subject_id_selection(which_study,'behavior');
        if which_study == 0
            root = ['B:' filesep 'resultats' filesep 'multiseq_april_may2016' filesep];
        elseif which_study == 1
            root = ['B:' filesep 'resultats' filesep 'MBB_june2016' filesep];
        end
    case 2
        [subject_id, NS] = MS2_subject_id_selection('pupil');
        root = ['F:' filesep 'MBB_MotiScan2_march2017' filesep];
end

%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    cd([root, sub_nm, filesep,'eye_analysis',filesep,'preprocessed_files',filesep])
    
    %% delete all eye data
    eyedata_filenm = ls('*.mat');
    for iFile = size(eyedata_filenm,1):-1:1
        delete(eyedata_filenm(iFile,:));
    end
    
    pic_filenm = ls('*.png');
    for iFile = size(pic_filenm,1):-1:1
        delete(pic_filenm(iFile,:));
    end
    
    %% keep track of where you are
    disp(['Subject ',sub_nm, ' data deleted']);
end

end

