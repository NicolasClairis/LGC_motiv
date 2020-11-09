function [] = preproc_eyeGazePupil_MotiScan(which_study)
%preproc_eyeGazePupil function to make the preprocessing of the eye gaze
%and pupil data before using analytical scripts to check any correlation
%with the variables of interest.
%
% See also read_eyelink_1_asc_kk.m


%% select subjects to preprocess
if exist('which_study','var') == 0 || isempty(which_study)
    which_study = input('Sequence comparison (0) or MBB june2016(1) or both (2) or MS-2 (3)?');
end

if ismember(which_study,[0,1,2])
    [subject_id, ~, NS_mseq, NS_MBB] = MS1_subject_id_selection(which_study, 'pupil_rtg');
elseif which_study == 3
    [subject_id, NS_MS2] = MS2_subject_id_selection('pupil');
end

if which_study == 0
    min_st = 0;
    max_st = 0;
elseif which_study == 1
    min_st = 1;
    max_st = 1;
elseif which_study == 2
    min_st = 0;
    max_st = 1;
elseif which_study == 3
    min_st = 3;
    max_st = 3;
end

%% define main acquisition parameters (eye-tracker sampling rate, screen
% size, task names)
if which_study <= 3
    % sampling rate in Hz for CENIR fMRI eye-tracker
    samplingRate = 1000;
    
    % define screen coordinates
    x_screen = 1024;
    y_screen = 768;

    if which_study < 3
        taskName = {'rtgs','c1d','cRE'};
    elseif which_study == 3
        taskName = {'RL','grip','stroop'};
    end
    nTasks = length(taskName);
end

disp_gr = 1; % display resume of pupil graph for each run of each subject (1) or not (0)?

% define folders containing analysis scripts
scripts_folder = ['C:',filesep,'Users',filesep,'nicolas.clairis',filesep,...
    'Desktop',filesep,'resultats',filesep,'analysis_scripts',filesep,'eye_tracking_functions',filesep];
preproc_Git_folder = ['C:',filesep,'Users',filesep,'nicolas.clairis',filesep,...
    'Documents',filesep,'GitHub',filesep,'proud-pupil',filesep];

%% list of data to skip
list_skipped_files_mseq_aprilMay2016_0 = {'MSs_10r3.asc','MSs_10r7.asc','MSs_10r8.asc','MSs_10r9.asc'};
% s10r3 mseq: eye not captured at all during this run
list_skipped_files_MBBjune2016_0 = {'MS_s02r3.asc',...
    'MS_s07r7.asc','MS_s07r8.asc','MS_s07r9.asc'};% s02r3 MBBjune: corrupted file: TTLs missing;

%% loop across studies/subjects/tasks
for iStudy = min_st:max_st
    if ismember(which_study,[0,1,2])
        NS_study  = NS_mseq*(iStudy == 0) + NS_MBB*(iStudy == 1);
    elseif which_study == 3
        NS_study  = NS_MS2;
    end
    
    %% path definition
    switch iStudy
        case 0
            path = 'multiseq_april_may2016';
        case 1
            path = 'MBB_june2016';
        case 3
            path = 'MS2_march2017';
    end
    switch which_study
        case {0,1,2}
            root = ['B:',filesep,'resultats',filesep, path, filesep];
        case 3
            root = ['F:',filesep,'MBB_MotiScan2_march2017',filesep];
    end
    
    switch iStudy
        case 0
            list_skipped_files = list_skipped_files_mseq_aprilMay2016_0;
        case 1
            list_skipped_files = list_skipped_files_MBBjune2016_0;
        case 3
            list_skipped_files = {};
    end
    
    if which_study == 2 && iStudy == 0
        list_skipped_files2 = [list_skipped_files_mseq_aprilMay2016_0, list_skipped_files_MBBjune2016_0];
    end
    
    for iSubject = 1:NS_study
        if which_study == 2 && iStudy == 1
            sub_idx = NS_mseq + iSubject;
        else
            sub_idx = iSubject;
        end
        subName = subject_id{sub_idx};
        % take particular subject number
        if strcmp(path,'MBB_june2016')
            subid = subName(2:3);
        elseif strcmp(path,'multiseq_april_may2016') || strcmp(path,'MS2_march2017')
            if strcmp(subName(3),'_')
                subid = subName(2);
            elseif strcmp(subName(4),'_')
                subid = subName(2:3);
            end
        end
        
        % store data into one folder for each subject (and create it if it
        % does not exist already)
        eyeAnalysis_folder = [root, subName, filesep, 'eye_analysis', filesep];
        if exist(eyeAnalysis_folder,'dir') == 0
            mkdir(eyeAnalysis_folder);
        end
        % same for preprocessing folder
        eyePreproc_folder = [eyeAnalysis_folder, 'preprocessed_files', filesep];
        if exist(eyePreproc_folder,'dir') == 0
            mkdir(eyePreproc_folder);
        end
        
        %% loop through tasks
        for iTask = 1:nTasks
            if which_study < 3
                runs = (1:3)*(iTask == 1) + (4:6)*(iTask == 2) + (7:9)*(iTask == 3); % 3 first runs = task 1 = ratings, 4 to 6 = task 2 = 1D-choices, 7 to 9 = task 3 = RE-choices
            else
                runs = MS2_task_runs_extraction(taskName{iTask}, subName);
            end
            for iRun = runs
                runName = num2str(iRun);
                
                % extract raw eye signal
                eyeDataFolder = [root, subName, filesep, 'eye_data', filesep];
                
                % extract file name
                if iStudy == 0
                    if length(subid) < 2
                        fileName = ['MSs_s',subid,'r',runName,'.asc'];
                    elseif length(subid) == 2
                        fileName = ['MSs_',subid,'r',runName,'.asc'];
                    end
                elseif iStudy == 1
                    fileName = ['MS_s',subid,'r',runName,'.asc'];
                elseif iStudy == 3
                    fileName = ['MS2s',subid,'r',runName,'.asc'];
                end
                full_fileName = [eyeDataFolder, fileName];
                
                if ~ismember(fileName, list_skipped_files)
                    
                    % extract raw signal of eyes
                    asc = read_eyelink_1_asc_kk(full_fileName);
                    
                    %% check TTL and filter first wrong TTLs in corresponding files
                    % if no TTL or problem with TTL too big to be fixed, the run/subject should
                    % be abandoned and thus there is no use in spending time with this run/subject
                    
                    % extract TTL from asc file
                    TTL = cell2mat(struct2cell(asc.button)');
                    TTL = TTL(TTL(:,3) == 1,1);
                    
                    % fix TTL for each run and each subject
                    [TTL] = TTL_check_for_eye(fileName, TTL, iStudy, subName, runName);
                    
                    
                    %% use preprocessing script
                    [eyeTime, X_eyeCoord, Y_eyeCoord, pupil_diam, z_pupil_diam, nanEpisodes] = preproc_eyeGazePupil(asc, TTL(1), x_screen, y_screen, samplingRate, disp_gr);
                    
%                     %% modified for Antonius part
%                     cd('C:\Users\nicolas.clairis\Desktop\resultats\analysis_scripts\eye_tracking_functions\Antonius');
%                     timestamps = asc.dat(1,:).*1000';
%                     pupil_size = asc.dat(4,:)';
%                     valid_recording = asc.dat(5,:)'; % 0 when eye lost
%                     [pupil_size, timestamps, sampling_factor] = pupil_preprocessing_canonical_v01(timestamps, pupil_size, valid_recording);
%                     saveNameFiles = [eyePreproc_folder, 'eye_preproc_FOR_ANTONIUS_', subName,'_',taskName{iTask},'_run',runName,'.mat'];
%                     save(saveNameFiles, 'timestamps','pupil_size','sampling_factor');
                    
                    %% save the resulting files
                    % save resulting files (if data usable)
                    if ~isempty(eyeTime) && ~isempty(X_eyeCoord) && ~isempty(Y_eyeCoord) && ~isempty(pupil_diam) && ~isempty(z_pupil_diam)
                        saveNameFiles = [eyePreproc_folder, 'eye_preproc_', subName,'_',taskName{iTask},'_run',runName,'.mat'];
                        save(saveNameFiles, 'eyeTime','X_eyeCoord', 'Y_eyeCoord', 'pupil_diam', 'z_pupil_diam', 'nanEpisodes');
                        disp([subName,' ',taskName{iTask},' run',runName,' done.']);
                    else
                        warning([subName,' ',taskName{iTask},' run',runName,' skipped because data not usable.']);
                        list_skipped_files = [list_skipped_files, fileName];
                        list_skipped_files2 = [list_skipped_files2, fileName];
                    end
                    % save resulting resume figure
                    if disp_gr == 1
                        saveNameFig = [eyePreproc_folder, 'eye_preproc_', subName,'_',taskName{iTask},'_run',runName,'.png'];
                        saveas(gcf,saveNameFig)
                        close all;
                    end
                    
                end % run filter
            end % run loop
        end % task loop
    end % subject loop
    
    %% save files ignored because signal was too bad
    save([root, filesep, 'behavior_summary', filesep 'eye' filesep 'list_skipped_runs.mat'], 'list_skipped_files');
    
    if which_study == 2 && iStudy == 1
        list_skipped_files = list_skipped_files2;
        concat_path = ['B:',filesep,'resultats',filesep, 'concat_MBB_multiseq_2016', filesep, 'behavior_summary', filesep 'eye' filesep];
        save([concat_path,  'list_skipped_runs.mat'], 'list_skipped_files');
    end

end % study loop

end % function end

