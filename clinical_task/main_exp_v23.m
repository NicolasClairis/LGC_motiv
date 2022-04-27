% main script for visits 2 and 3 of the clinical trial
% 1) re-do the training of the task to remind them how it works
% 2) perform the reward/effort tradeoff incentivized task
% 3) train for the "social" task
% 4) perform the "social" task


%% clean workspace before starting
sca;
clearvars;
close all;
instrreset; % Disconnect and delete all instrument objects
clc;

%% working directories
% launch within the folder where scripts are stored or will not work
cd ..
main_folder                 = [pwd filesep]; % you have to be sure that you are in the correct path when you launch the script
main_task_folder            = [main_folder, 'LGC_Motiv_task' filesep];
results_folder              = [main_folder, 'LGC_Motiv_results' filesep];
% BioPac_folder               = [main_folder, 'BioPac_functions' filesep];
% pics_folder                 = [main_task_folder, 'Coin_PNG', filesep];
Matlab_DIY_functions_folder = [main_folder, 'Matlab_DIY_functions', filesep];

% add personal functions (needed for PTB opening at least)
addpath(genpath(main_task_folder));
% addpath(BioPac_folder);
addpath(Matlab_DIY_functions_folder);

% create results folder if no -subject has been acquired yet
if ~exist(results_folder,'dir')
    mkdir(results_folder);
end
% go back to folder with scripts
cd(main_task_folder);

%% Define subject ID

% Insert the initials, the number of the participants
iSubject = [];
langue = 'f';
IRMdisp = 0; % defines the screen parameters (0 for training screen, 1 for fMRI screen)
IRMbuttons = 1; % defines the buttons to use (1 = same as in fMRI)
testing_script = 0; % use all computer resources (particularly for mental calibration)
while isempty(iSubject) || length(iSubject) ~= 3
    % repeat until all questions are answered
    info = inputdlg({'Subject CID (XXX)','p/m','visit?'});
    [iSubject, p_or_m, visit] = info{[1,2,3]};
    %     warning('when real experiment starts, remember to block in french and IRM = 1');
end
% Create subjectCodeName which is used as a file saving name
subjectCodeName = strcat('CID',iSubject);
subResultFolder = [results_folder, subjectCodeName, filesep,'behavior',filesep];
if ~exist(subResultFolder,'dir')
    mkdir(subResultFolder);
end


if strcmp(p_or_m,'m') == 0 && strcmp(p_or_m,'p') == 0
    error('this letter has no definition');
end

n_visit = str2double(visit);
if ~ismember(n_visit,[1,2,3])
    error(['visit number other than 1, 2 or 3 not possible and you defined visit number as ',visit]);
end
visit_nm = ['_v',visit];

file_nm_training_Em = ['training_data_Em_CID',num2str(iSubject),visit_nm];
file_nm_training_Ep = ['training_data_Ep_CID',num2str(iSubject),visit_nm];
file_nm = ['training_data_CID',num2str(iSubject),visit_nm];
fullFileNm = [subResultFolder, file_nm,'.mat'];
if exist(fullFileNm,'file')
    error([fullFileNm,' file already exists. If the training script was launched and crashed, ',...
        ' please consider renaming the files already saved to avoid losing data.']);
end
Ep_visitCalib_filenm = [subResultFolder,subjectCodeName,'_physicalCalib',visit_nm,'.mat'];
Em_visitCalib_filenm = [subResultFolder,subjectCodeName,'_mentalCalib',visit_nm,'.mat'];
% baseline max
Ep_calib_filenm = [subResultFolder,subjectCodeName,'_physicalCalib_v1.mat'];
Em_calib_filenm = [subResultFolder,subjectCodeName,'_mentalCalib_v1.mat'];
% IP file name
file_nm_IP = ['delta_IP_CID',num2str(iSubject),visit_nm];
% convert subject CID into number (only if used to perform actual task)
if ischar(iSubject)
    iSubject = str2double(iSubject);
end