function[subject_id, NS] = preproc_already_made(computerRoot, study_nm, subject_id, smKernel)
% [subject_id, NS] = preproc_already_made(computerRoot, study_nm, subject_id, smKernel)
% preproc_already_made will check if the preprocessing has already been
% done or not for each subject for the preprocessing kernel entered as
% input.
%
% INPUTS
% computerRoot: root path where subject data is stored
%
% study_nm: study name
% 'study1': dmPFC/aINs MRS study
% 'clinical_study': clinical study
% 
% subject_id: list of subjects
%
% smKernel: smoothing kernel used (in mm)
%
% OUTPUTS
% subject_id: list of subjects after filtering those already preprocessed
%
% NS: number of subjects after filtering those already preprocessed

switch study_nm
    case {'study1','study2','fMRI_pilots'}
        studyPath = [computerRoot, study_nm, filesep];
    case 'study2_pilots'
        studyPath = [fullfile(computerRoot, 'study2','pilots','fMRI_pilots'),filesep];
end

%% check for each subject if data already exists
NS = length(subject_id);
subNotAlreadyPreprocessed = true(1,NS);
for iS = 1:NS
    switch study_nm
        case {'study1','study2'}
            sub_nm = ['CID',subject_id{iS}];
        case {'fMRI_pilots','study2_pilots'}
            sub_nm = subject_id{iS};
    end
    subj_scans_folder = [studyPath, sub_nm, filesep, 'fMRI_scans' filesep];
    cd(subj_scans_folder);
    subj_scan_folders_names = ls('*run*'); % takes all functional runs folders
    n_runs = size(subj_scan_folders_names,1);
    dataPreprocessed_tmp = false;
    for iRun = 1:n_runs
        runPath = [subj_scans_folder, subj_scan_folders_names(iRun,:),filesep,'preproc_sm_',num2str(smKernel),'mm'];
        if exist(runPath,'dir')
            disp([runPath,' already existed and was ignored from preprocessing']);
            % if any of the runs has already been preprocessed => ignore
            % the whole subject from preprocessing
            if dataPreprocessed_tmp == false
                dataPreprocessed_tmp = true; % if the subject has already 
                % been preprocessed for any run, the script will ignore the
                % subject as a whole => turn temporary variable
                % dataPreprocessed_tmp to true even if only one run is
                % concerned
                subNotAlreadyPreprocessed(iS) = false;
            end
        end % preprocessing folder already exists
    end % run loop
end % subject loop

%% back to main folder
cd(studyPath);


%% extract new list of subjects
subject_id = subject_id(subNotAlreadyPreprocessed);
NS = length(subject_id);

end % function