% function [  ] = eye_bins()
%eye_ratings_bins separates pupil diameter into bins of 20% to see whether
%there is a link with the difficulty/ratings

clear all; close all; clc;
% enter subjects identification
which_study = input('Sequence comparison (0) or MBB june2016(1) ?');
% which_study = 0;
if which_study == 0
    subject_id = {'s1_250416','s2_260416','s3_070516','s5_070516',...
        's6_100516','s7_110516','s8_200516','s9_220516','s10_220516',...
        's11_220516','s12_220516','s13_220516','s14_100616','s15_170616'};
    path = 'multiseq_april_may2016';
elseif which_study == 1
    subject_id = {'s00_030616','s01_050616','s02_050616','s03_050616',...
        's04_050616','s05_050616','s06_090616','s07_090616',...
        's08_090616','s09_090616','s10_090616','s11_110616',...
        's12_110616','s13_110616','s14_110616','s15_110616',...
        's16_120616','s17_120616','s19_190616','s20_190616',...
        's21_190616','s23_190616','s24_190616'};
    % s18 and s22 excluded because signal was too bad
    path = 'MBB_june2016';
else
    disp('error in study number, retry');
    return;
end
NS = length(subject_id);
root = ['B:\resultats\',path,'\'];
% ratings or choice task?
whichtask = input('Ratings (0) or choice (1) task?');
if whichtask == 0
    task = 'ratings';
elseif whichtask == 1
    task = 'choice_1D';
else
    disp('error in task number');
    return;
end
varname = [num2str(whichtask+1),'_',path,'_eye_',task,'_resume_',num2str(NS),'subs.mat'];

scripts_folder = [fullfile('C:','Users','nicolas.clairis','Desktop','resultats','analysis_scripts','eye_tracking_functions'), filesep];

%% load pupil diameter data for all the subjects with onset and RT lock
if whichtask == 0
    % load data saved from eye_ratings.m
    loadStruct = load([root,filesep,'behavior_summary',filesep,varname],'pupil_lock_onset','pupil_lock_RTdec','pupil_lock_RTval');
    pupil_E_lock_onset = loadStruct.pupil_lock_onset.E;
    pupil_E_lock_RTdec = loadStruct.pupil_lockRTdec.E;
    pupil_E_lock_RTval = loadStruct.pupil_lockRTval.E;
    pupil_R_lock_onset = loadStruct.pupil_lock_onset.R;
    pupil_R_lock_RTdec = loadStruct.pupil_lockRTdec.R;
    pupil_R_lock_RTval = loadStruct.pupil_lockRTval.R;
    pupil_Rim_lock_onset = loadStruct.pupil_lock_onset.Rim;
    pupil_Rim_lock_RTdec = loadStruct.pupil_lockRTdec.Rim;
    pupil_Rim_lock_RTval = loadStruct.pupil_lockRTval.Rim;
elseif whichtask == 1
    % % 1st possibility: relaunch the script:
%     [~,...
%         ~, ~,...
%         pupil_lock_onset, pupil_lock_RT] = eye_choice_one_dim(0, 2, 0);
    % faster: load data saved from this script
    loadStruct = load([root, filesep, 'behavior_summary', filesep, varname],'pupil_lock_onset','pupil_lock_RT');
    pupil_EE_lock_onset     = loadStruct.pupil_lock_onset.EE;
    pupil_EE_lock_RT        = loadStruct.pupil_lock_RT.EE;
    pupil_RR_lock_onset     = loadStruct.pupil_lock_onset.RR;
    pupil_RR_lock_RT        = loadStruct.pupil_lock_RT.RR;
    pupil_RimRim_lock_onset = loadStruct.pupil_lock_onset.RimRim;
    pupil_RimRim_lock_RT    = loadStruct.pupil_lock_RT.RimRim;
end
% variables
nb_trials = 72;
nb_trials_per_run = nb_trials/3;
% rating variables
binvar_E = zeros(nb_trials,NS); % variable used for binning (ratings/difficulty)
binvar_R = zeros(nb_trials,NS);
binvar_Rim = zeros(nb_trials, NS);
% bins for different rating values
if whichtask == 0
    timelength = size(pupil_E_lock_onset,2);
elseif whichtask == 1
    timelength = size(pupil_EE_lock_onset,2);
end
[pupil_E_lock_onset_0_20, pupil_E_lock_onset_20_40, pupil_E_lock_onset_40_60, pupil_E_lock_onset_60_80, pupil_E_lock_onset_80_100,...
    pupil_R_lock_onset_0_20, pupil_R_lock_onset_20_40, pupil_R_lock_onset_40_60, pupil_R_lock_onset_60_80, pupil_R_lock_onset_80_100,...
    pupil_Rim_lock_onset_0_20, pupil_Rim_lock_onset_20_40, pupil_Rim_lock_onset_40_60, pupil_Rim_lock_onset_60_80, pupil_Rim_lock_onset_80_100] = deal( NaN(nb_trials,timelength,NS) );
if whichtask == 0
    [pupil_E_lock_RT_0_20, pupil_E_lock_RT_20_40, pupil_E_lock_RT_40_60, pupil_E_lock_RT_60_80, pupil_E_lock_RT_80_100,...
        pupil_R_lock_RT_0_20, pupil_R_lock_RT_20_40, pupil_R_lock_RT_40_60, pupil_R_lock_RT_60_80, pupil_R_lock_RT_80_100,...
        pupil_Rim_lock_RT_0_20, pupil_Rim_lock_RT_20_40, pupil_Rim_lock_RT_40_60, pupil_Rim_lock_RT_60_80, pupil_Rim_lock_RT_80_100,...
        pupil_E_lock_RTval_0_20, pupil_E_lock_RTval_20_40, pupil_E_lock_RTval_40_60, pupil_E_lock_RTval_60_80, pupil_E_lock_RTval_80_100,...
        pupil_R_lock_RTval_0_20, pupil_R_lock_RTval_20_40, pupil_R_lock_RTval_40_60, pupil_R_lock_RTval_60_80, pupil_R_lock_RTval_80_100,...
        pupil_Rim_lock_RTval_0_20, pupil_Rim_lock_RTval_20_40, pupil_Rim_lock_RTval_40_60, pupil_Rim_lock_RTval_60_80, pupil_Rim_lock_RTval_80_100] = deal( NaN(nb_trials,timelength,NS) );
elseif whichtask == 1
    [pupil_E_lock_RT_0_20, pupil_E_lock_RT_20_40, pupil_E_lock_RT_40_60, pupil_E_lock_RT_60_80, pupil_E_lock_RT_80_100,...
        pupil_R_lock_RT_0_20, pupil_R_lock_RT_20_40, pupil_R_lock_RT_40_60, pupil_R_lock_RT_60_80, pupil_R_lock_RT_80_100,...
        pupil_Rim_lock_RT_0_20, pupil_Rim_lock_RT_20_40, pupil_Rim_lock_RT_40_60, pupil_Rim_lock_RT_60_80, pupil_Rim_lock_RT_80_100] = deal( NaN(nb_trials,timelength,NS) );
end


for subject = 1:NS
    % take particular subject number
    if which_study == 0 % multiseq
        if strcmp(subject_id{subject}(3),'_') % subjects 1 to 9
            subid = subject_id{subject}(2);
        elseif strcmp(subject_id{subject}(4),'_') % subjects 10 to 15
            subid = subject_id{subject}(2:3);
        end
    elseif which_study == 1 % MBB june 2016
        if strcmp(subject_id{subject}(2:3), '00') || strcmp(subject_id{subject}(2), '0') == 0 % subjects 00 and 10 to 24
            subid = subject_id{subject}(2:3);
        elseif strcmp(subject_id{subject}(2), '0') && strcmp(subject_id{subject}(3), '0') == 0 % subjects 01 to 09
            subid = subject_id{subject}(3);
        end
    end
    sub_behavior_dir = [root subject_id{subject} filesep 'behavior', filesep];
    
    
    if whichtask == 0
        for run = 1:3 % ratings runs
            % values (=ratings) for each trial
            runname = num2str(run);
            first_idx = 1 + nb_trials_per_run*(run - 1);
            last_idx = nb_trials_per_run*run;
            % efforts
            x_e = dir([sub_behavior_dir, 'MBB_battery_ratingE_sub',subid,'_sess',runname,'_*','*min.mat']);
            load([sub_behavior_dir, x_e.name],'rating')
            binvar_E(first_idx:last_idx, subject) = rating;
            
            % reward text
            x_r = dir([sub_behavior_dir, 'MBB_battery_ratingR_sub',subid,'_sess',runname,'_*','*min.mat']);
            load([sub_behavior_dir, x_r.name],'rating')
            binvar_R(first_idx:last_idx, subject) = rating;
            
            % reward images
            x_rim = dir([sub_behavior_dir, 'MBB_battery_ratingR_im_sub',subid,'_sess',runname,'_*','*min.mat']);
            load([sub_behavior_dir, x_rim.name],'rating')
            binvar_Rim(first_idx:last_idx, subject) = rating;
        end
    elseif whichtask == 1
        for run = 4:6 % choice 1D runs
            % values (= difficulty) for each trial
            runname = num2str(run);
            first_idx = 1 + nb_trials_per_run*(run - 4);
            last_idx = nb_trials_per_run*(run - 3);
            % efforts
            x_ee = dir([sub_behavior_dir, 'MBB_battery_choiceE_sub',subid,'_sess',runname,'_*','*min.mat']);
            load([sub_behavior_dir, x_ee.name],'ratingleft_E','ratingright_E')
            binvar_E(first_idx:last_idx, subject) = abs(ratingleft_E - ratingright_E);
            
            % reward text
            x_rr = dir([sub_behavior_dir, 'MBB_battery_choiceR_sub',subid,'_sess',runname,'_*','*min.mat']);
            load([sub_behavior_dir, x_rr.name],'ratingleft_R','ratingright_R')
            binvar_R(first_idx:last_idx, subject) = abs(ratingleft_R - ratingright_R);
            
            % reward images
            x_rimrim = dir([sub_behavior_dir, 'MBB_battery_choiceR_im_sub',subid,'_sess',runname,'_*','*min.mat']);
            load([sub_behavior_dir, x_rimrim.name],'ratingleft_Rim','ratingright_Rim')
            binvar_Rim(first_idx:last_idx, subject) = abs(ratingleft_Rim - ratingright_Rim);
        end
    end
    
    
    %% use ratings to separate pupil values into bins of 20%
    for trial = 1:nb_trials
        if whichtask == 0
            % effort
            if binvar_E(trial, subject) > 0 && binvar_E(trial, subject) <= 20 % rating starts at 1
                pupil_E_lock_onset_0_20(trial,:,subject) = pupil_E_lock_onset(trial,:,subject);
                pupil_E_lock_RT_0_20(trial,:,subject) = pupil_E_lock_RTdec(trial,:,subject);
                pupil_E_lock_RTval_0_20(trial,:,subject) = pupil_E_lock_RTval(trial,:,subject);
            elseif binvar_E(trial, subject) > 20 && binvar_E(trial, subject) <= 40
                pupil_E_lock_onset_20_40(trial,:,subject) = pupil_E_lock_onset(trial,:,subject);
                pupil_E_lock_RT_20_40(trial,:,subject) = pupil_E_lock_RTdec(trial,:,subject);
                pupil_E_lock_RTval_20_40(trial,:,subject) = pupil_E_lock_RTval(trial,:,subject);
            elseif binvar_E(trial, subject) > 40 && binvar_E(trial, subject) <= 60
                pupil_E_lock_onset_40_60(trial,:,subject) = pupil_E_lock_onset(trial,:,subject);
                pupil_E_lock_RT_40_60(trial,:,subject) = pupil_E_lock_RTdec(trial,:,subject);
                pupil_E_lock_RTval_40_60(trial,:,subject) = pupil_E_lock_RTval(trial,:,subject);
            elseif binvar_E(trial, subject) > 60 && binvar_E(trial, subject) <= 80
                pupil_E_lock_onset_60_80(trial,:,subject) = pupil_E_lock_onset(trial,:,subject);
                pupil_E_lock_RT_60_80(trial,:,subject) = pupil_E_lock_RTdec(trial,:,subject);
                pupil_E_lock_RTval_60_80(trial,:,subject) = pupil_E_lock_RTval(trial,:,subject);
            elseif binvar_E(trial, subject) > 80 && binvar_E(trial, subject) <= 100
                pupil_E_lock_onset_80_100(trial,:,subject) = pupil_E_lock_onset(trial,:,subject);
                pupil_E_lock_RT_80_100(trial,:,subject) = pupil_E_lock_RTdec(trial,:,subject);
                pupil_E_lock_RTval_80_100(trial,:,subject) = pupil_E_lock_RTval(trial,:,subject);
            end
            % reward text
            if binvar_R(trial, subject) > 0 && binvar_R(trial, subject) <= 20 % rating starts at 1
                pupil_R_lock_onset_0_20(trial,:,subject) = pupil_R_lock_onset(trial,:,subject);
                pupil_R_lock_RT_0_20(trial,:,subject) = pupil_R_lock_RTdec(trial,:,subject);
                pupil_R_lock_RTval_0_20(trial,:,subject) = pupil_R_lock_RTval(trial,:,subject);
            elseif binvar_R(trial, subject) > 20 && binvar_R(trial, subject) <= 40
                pupil_R_lock_onset_20_40(trial,:,subject) = pupil_R_lock_onset(trial,:,subject);
                pupil_R_lock_RT_20_40(trial,:,subject) = pupil_R_lock_RTdec(trial,:,subject);
                pupil_R_lock_RTval_20_40(trial,:,subject) = pupil_R_lock_RTval(trial,:,subject);
            elseif binvar_R(trial, subject) > 40 && binvar_R(trial, subject) <= 60
                pupil_R_lock_onset_40_60(trial,:,subject) = pupil_R_lock_onset(trial,:,subject);
                pupil_R_lock_RT_40_60(trial,:,subject) = pupil_R_lock_RTdec(trial,:,subject);
                pupil_R_lock_RTval_40_60(trial,:,subject) = pupil_R_lock_RTval(trial,:,subject);
            elseif binvar_R(trial, subject) > 60 && binvar_R(trial, subject) <= 80
                pupil_R_lock_onset_60_80(trial,:,subject) = pupil_R_lock_onset(trial,:,subject);
                pupil_R_lock_RT_60_80(trial,:,subject) = pupil_R_lock_RTdec(trial,:,subject);
                pupil_R_lock_RTval_60_80(trial,:,subject) = pupil_R_lock_RTval(trial,:,subject);
            elseif binvar_R(trial, subject) > 80 && binvar_R(trial, subject) <= 100
                pupil_R_lock_onset_80_100(trial,:,subject) = pupil_R_lock_onset(trial,:,subject);
                pupil_R_lock_RT_80_100(trial,:,subject) = pupil_R_lock_RTdec(trial,:,subject);
                pupil_R_lock_RTval_80_100(trial,:,subject) = pupil_R_lock_RTval(trial,:,subject);
            end
            % reward images
            if binvar_Rim(trial, subject) > 0 && binvar_Rim(trial, subject) <= 20 % rating starts at 1
                pupil_Rim_lock_onset_0_20(trial,:,subject) = pupil_Rim_lock_onset(trial,:,subject);
                pupil_Rim_lock_RT_0_20(trial,:,subject) = pupil_Rim_lock_RTdec(trial,:,subject);
                pupil_Rim_lock_RTval_0_20(trial,:,subject) = pupil_Rim_lock_RTval(trial,:,subject);
            elseif binvar_Rim(trial, subject) > 20 && binvar_Rim(trial, subject) <= 40
                pupil_Rim_lock_onset_20_40(trial,:,subject) = pupil_Rim_lock_onset(trial,:,subject);
                pupil_Rim_lock_RT_20_40(trial,:,subject) = pupil_Rim_lock_RTdec(trial,:,subject);
                pupil_Rim_lock_RTval_20_40(trial,:,subject) = pupil_Rim_lock_RTval(trial,:,subject);
            elseif binvar_Rim(trial, subject) > 40 && binvar_Rim(trial, subject) <= 60
                pupil_Rim_lock_onset_40_60(trial,:,subject) = pupil_Rim_lock_onset(trial,:,subject);
                pupil_Rim_lock_RT_40_60(trial,:,subject) = pupil_Rim_lock_RTdec(trial,:,subject);
                pupil_Rim_lock_RTval_40_60(trial,:,subject) = pupil_Rim_lock_RTval(trial,:,subject);
            elseif binvar_Rim(trial, subject) > 60 && binvar_Rim(trial, subject) <= 80
                pupil_Rim_lock_onset_60_80(trial,:,subject) = pupil_Rim_lock_onset(trial,:,subject);
                pupil_Rim_lock_RT_60_80(trial,:,subject) = pupil_Rim_lock_RTdec(trial,:,subject);
                pupil_Rim_lock_RTval_60_80(trial,:,subject) = pupil_Rim_lock_RTval(trial,:,subject);
            elseif binvar_Rim(trial, subject) > 80 && binvar_Rim(trial, subject) <= 100
                pupil_Rim_lock_onset_80_100(trial,:,subject) = pupil_Rim_lock_onset(trial,:,subject);
                pupil_Rim_lock_RT_80_100(trial,:,subject) = pupil_Rim_lock_RTdec(trial,:,subject);
                pupil_Rim_lock_RTval_80_100(trial,:,subject) = pupil_Rim_lock_RTval(trial,:,subject);
            end
        elseif whichtask == 1
            % effort
            if binvar_E(trial, subject) > 0 && binvar_E(trial, subject) <= 20 % rating starts at 1
                pupil_E_lock_onset_0_20(trial,:,subject) = pupil_EE_lock_onset(trial,:,subject);
                pupil_E_lock_RT_0_20(trial,:,subject) = pupil_EE_lock_RT(trial,:,subject);
            elseif binvar_E(trial, subject) > 20 && binvar_E(trial, subject) <= 40
                pupil_E_lock_onset_20_40(trial,:,subject) = pupil_EE_lock_onset(trial,:,subject);
                pupil_E_lock_RT_20_40(trial,:,subject) = pupil_EE_lock_RT(trial,:,subject);
            elseif binvar_E(trial, subject) > 40 && binvar_E(trial, subject) <= 60
                pupil_E_lock_onset_40_60(trial,:,subject) = pupil_EE_lock_onset(trial,:,subject);
                pupil_E_lock_RT_40_60(trial,:,subject) = pupil_EE_lock_RT(trial,:,subject);
            elseif binvar_E(trial, subject) > 60 && binvar_E(trial, subject) <= 80
                pupil_E_lock_onset_60_80(trial,:,subject) = pupil_EE_lock_onset(trial,:,subject);
                pupil_E_lock_RT_60_80(trial,:,subject) = pupil_EE_lock_RT(trial,:,subject);
            elseif binvar_E(trial, subject) > 80 && binvar_E(trial, subject) <= 100
                pupil_E_lock_onset_80_100(trial,:,subject) = pupil_EE_lock_onset(trial,:,subject);
                pupil_E_lock_RT_80_100(trial,:,subject) = pupil_EE_lock_RT(trial,:,subject);
            end
            % reward text
            if binvar_R(trial, subject) > 0 && binvar_R(trial, subject) <= 20 % rating starts at 1
                pupil_R_lock_onset_0_20(trial,:,subject) = pupil_RR_lock_onset(trial,:,subject);
                pupil_R_lock_RT_0_20(trial,:,subject) = pupil_RR_lock_RT(trial,:,subject);
            elseif binvar_R(trial, subject) > 20 && binvar_R(trial, subject) <= 40
                pupil_R_lock_onset_20_40(trial,:,subject) = pupil_RR_lock_onset(trial,:,subject);
                pupil_R_lock_RT_20_40(trial,:,subject) = pupil_RR_lock_RT(trial,:,subject);
            elseif binvar_R(trial, subject) > 40 && binvar_R(trial, subject) <= 60
                pupil_R_lock_onset_40_60(trial,:,subject) = pupil_RR_lock_onset(trial,:,subject);
                pupil_R_lock_RT_40_60(trial,:,subject) = pupil_RR_lock_RT(trial,:,subject);
            elseif binvar_R(trial, subject) > 60 && binvar_R(trial, subject) <= 80
                pupil_R_lock_onset_60_80(trial,:,subject) = pupil_RR_lock_onset(trial,:,subject);
                pupil_R_lock_RT_60_80(trial,:,subject) = pupil_RR_lock_RT(trial,:,subject);
            elseif binvar_R(trial, subject) > 80 && binvar_R(trial, subject) <= 100
                pupil_R_lock_onset_80_100(trial,:,subject) = pupil_RR_lock_onset(trial,:,subject);
                pupil_R_lock_RT_80_100(trial,:,subject) = pupil_RR_lock_RT(trial,:,subject);
            end
            % reward images
            if binvar_Rim(trial, subject) > 0 && binvar_Rim(trial, subject) <= 20 % rating starts at 1
                pupil_Rim_lock_onset_0_20(trial,:,subject) = pupil_RimRim_lock_onset(trial,:,subject);
                pupil_Rim_lock_RT_0_20(trial,:,subject) = pupil_RimRim_lock_RT(trial,:,subject);
            elseif binvar_Rim(trial, subject) > 20 && binvar_Rim(trial, subject) <= 40
                pupil_Rim_lock_onset_20_40(trial,:,subject) = pupil_RimRim_lock_onset(trial,:,subject);
                pupil_Rim_lock_RT_20_40(trial,:,subject) = pupil_RimRim_lock_RT(trial,:,subject);
            elseif binvar_Rim(trial, subject) > 40 && binvar_Rim(trial, subject) <= 60
                pupil_Rim_lock_onset_40_60(trial,:,subject) = pupil_RimRim_lock_onset(trial,:,subject);
                pupil_Rim_lock_RT_40_60(trial,:,subject) = pupil_RimRim_lock_RT(trial,:,subject);
            elseif binvar_Rim(trial, subject) > 60 && binvar_Rim(trial, subject) <= 80
                pupil_Rim_lock_onset_60_80(trial,:,subject) = pupil_RimRim_lock_onset(trial,:,subject);
                pupil_Rim_lock_RT_60_80(trial,:,subject) = pupil_RimRim_lock_RT(trial,:,subject);
            elseif binvar_Rim(trial, subject) > 80 && binvar_Rim(trial, subject) <= 100
                pupil_Rim_lock_onset_80_100(trial,:,subject) = pupil_RimRim_lock_onset(trial,:,subject);
                pupil_Rim_lock_RT_80_100(trial,:,subject) = pupil_RimRim_lock_RT(trial,:,subject);
            end
        end
    end
end

% mean for pupil arranged by bins
% pupil mean all trials for each subject
% onset locked
mean_pupil_E_lock_onset_0_20 = nanmean(pupil_E_lock_onset_0_20,1);
mean_pupil_E_lock_onset_20_40 = nanmean(pupil_E_lock_onset_20_40,1);
mean_pupil_E_lock_onset_40_60 = nanmean(pupil_E_lock_onset_40_60,1);
mean_pupil_E_lock_onset_60_80 = nanmean(pupil_E_lock_onset_60_80,1);
mean_pupil_E_lock_onset_80_100 = nanmean(pupil_E_lock_onset_80_100,1);
mean_pupil_R_lock_onset_0_20 = nanmean(pupil_R_lock_onset_0_20,1);
mean_pupil_R_lock_onset_20_40 = nanmean(pupil_R_lock_onset_20_40,1);
mean_pupil_R_lock_onset_40_60 = nanmean(pupil_R_lock_onset_40_60,1);
mean_pupil_R_lock_onset_60_80 = nanmean(pupil_R_lock_onset_60_80,1);
mean_pupil_R_lock_onset_80_100 = nanmean(pupil_R_lock_onset_80_100,1);
mean_pupil_Rim_lock_onset_0_20 = nanmean(pupil_Rim_lock_onset_0_20,1);
mean_pupil_Rim_lock_onset_20_40 = nanmean(pupil_Rim_lock_onset_20_40,1);
mean_pupil_Rim_lock_onset_40_60 = nanmean(pupil_Rim_lock_onset_40_60,1);
mean_pupil_Rim_lock_onset_60_80 = nanmean(pupil_Rim_lock_onset_60_80,1);
mean_pupil_Rim_lock_onset_80_100 = nanmean(pupil_Rim_lock_onset_80_100,1);
% RT 1st press (whichtask=0)/RT choice locked (whichtask=1)
mean_pupil_E_lock_RT_0_20 = nanmean(pupil_E_lock_RT_0_20,1);
mean_pupil_E_lock_RT_20_40 = nanmean(pupil_E_lock_RT_20_40,1);
mean_pupil_E_lock_RT_40_60 = nanmean(pupil_E_lock_RT_40_60,1);
mean_pupil_E_lock_RT_60_80 = nanmean(pupil_E_lock_RT_60_80,1);
mean_pupil_E_lock_RT_80_100 = nanmean(pupil_E_lock_RT_80_100,1);
mean_pupil_R_lock_RT_0_20 = nanmean(pupil_R_lock_RT_0_20,1);
mean_pupil_R_lock_RT_20_40 = nanmean(pupil_R_lock_RT_20_40,1);
mean_pupil_R_lock_RT_40_60 = nanmean(pupil_R_lock_RT_40_60,1);
mean_pupil_R_lock_RT_60_80 = nanmean(pupil_R_lock_RT_60_80,1);
mean_pupil_R_lock_RT_80_100 = nanmean(pupil_R_lock_RT_80_100,1);
mean_pupil_Rim_lock_RT_0_20 = nanmean(pupil_Rim_lock_RT_0_20,1);
mean_pupil_Rim_lock_RT_20_40 = nanmean(pupil_Rim_lock_RT_20_40,1);
mean_pupil_Rim_lock_RT_40_60 = nanmean(pupil_Rim_lock_RT_40_60,1);
mean_pupil_Rim_lock_RT_60_80 = nanmean(pupil_Rim_lock_RT_60_80,1);
mean_pupil_Rim_lock_RT_80_100 = nanmean(pupil_Rim_lock_RT_80_100,1);

if whichtask == 0
    % RT validation
    mean_pupil_E_lock_RTval_0_20 = nanmean(pupil_E_lock_RTval_0_20,1);
    mean_pupil_E_lock_RTval_20_40 = nanmean(pupil_E_lock_RTval_20_40,1);
    mean_pupil_E_lock_RTval_40_60 = nanmean(pupil_E_lock_RTval_40_60,1);
    mean_pupil_E_lock_RTval_60_80 = nanmean(pupil_E_lock_RTval_60_80,1);
    mean_pupil_E_lock_RTval_80_100 = nanmean(pupil_E_lock_RTval_80_100,1);
    mean_pupil_R_lock_RTval_0_20 = nanmean(pupil_R_lock_RTval_0_20,1);
    mean_pupil_R_lock_RTval_20_40 = nanmean(pupil_R_lock_RTval_20_40,1);
    mean_pupil_R_lock_RTval_40_60 = nanmean(pupil_R_lock_RTval_40_60,1);
    mean_pupil_R_lock_RTval_60_80 = nanmean(pupil_R_lock_RTval_60_80,1);
    mean_pupil_R_lock_RTval_80_100 = nanmean(pupil_R_lock_RTval_80_100,1);
    mean_pupil_Rim_lock_RTval_0_20 = nanmean(pupil_Rim_lock_RTval_0_20,1);
    mean_pupil_Rim_lock_RTval_20_40 = nanmean(pupil_Rim_lock_RTval_20_40,1);
    mean_pupil_Rim_lock_RTval_40_60 = nanmean(pupil_Rim_lock_RTval_40_60,1);
    mean_pupil_Rim_lock_RTval_60_80 = nanmean(pupil_Rim_lock_RTval_60_80,1);
    mean_pupil_Rim_lock_RTval_80_100 = nanmean(pupil_Rim_lock_RTval_80_100,1);
end

%% pupil
% mean the subjects
% onset locked
allsubs_mean_pupil_E_lock_onset_0_20 = nanmean(mean_pupil_E_lock_onset_0_20,3);
allsubs_mean_pupil_E_lock_onset_20_40 = nanmean(mean_pupil_E_lock_onset_20_40,3);
allsubs_mean_pupil_E_lock_onset_40_60 = nanmean(mean_pupil_E_lock_onset_40_60,3);
allsubs_mean_pupil_E_lock_onset_60_80 = nanmean(mean_pupil_E_lock_onset_60_80,3);
allsubs_mean_pupil_E_lock_onset_80_100 = nanmean(mean_pupil_E_lock_onset_80_100,3);
allsubs_mean_pupil_R_lock_onset_0_20 = nanmean(mean_pupil_R_lock_onset_0_20,3);
allsubs_mean_pupil_R_lock_onset_20_40 = nanmean(mean_pupil_R_lock_onset_20_40,3);
allsubs_mean_pupil_R_lock_onset_40_60 = nanmean(mean_pupil_R_lock_onset_40_60,3);
allsubs_mean_pupil_R_lock_onset_60_80 = nanmean(mean_pupil_R_lock_onset_60_80,3);
allsubs_mean_pupil_R_lock_onset_80_100 = nanmean(mean_pupil_R_lock_onset_80_100,3);
allsubs_mean_pupil_Rim_lock_onset_0_20 = nanmean(mean_pupil_Rim_lock_onset_0_20,3);
allsubs_mean_pupil_Rim_lock_onset_20_40 = nanmean(mean_pupil_Rim_lock_onset_20_40,3);
allsubs_mean_pupil_Rim_lock_onset_40_60 = nanmean(mean_pupil_Rim_lock_onset_40_60,3);
allsubs_mean_pupil_Rim_lock_onset_60_80 = nanmean(mean_pupil_Rim_lock_onset_60_80,3);
allsubs_mean_pupil_Rim_lock_onset_80_100 = nanmean(mean_pupil_Rim_lock_onset_80_100,3);
% RT 1st press (whichtask=0)/RT choice locked (whichtask=1)
allsubs_mean_pupil_E_lock_RT_0_20 = nanmean(mean_pupil_E_lock_RT_0_20,3);
allsubs_mean_pupil_E_lock_RT_20_40 = nanmean(mean_pupil_E_lock_RT_20_40,3);
allsubs_mean_pupil_E_lock_RT_40_60 = nanmean(mean_pupil_E_lock_RT_40_60,3);
allsubs_mean_pupil_E_lock_RT_60_80 = nanmean(mean_pupil_E_lock_RT_60_80,3);
allsubs_mean_pupil_E_lock_RT_80_100 = nanmean(mean_pupil_E_lock_RT_80_100,3);
allsubs_mean_pupil_R_lock_RT_0_20 = nanmean(mean_pupil_R_lock_RT_0_20,3);
allsubs_mean_pupil_R_lock_RT_20_40 = nanmean(mean_pupil_R_lock_RT_20_40,3);
allsubs_mean_pupil_R_lock_RT_40_60 = nanmean(mean_pupil_R_lock_RT_40_60,3);
allsubs_mean_pupil_R_lock_RT_60_80 = nanmean(mean_pupil_R_lock_RT_60_80,3);
allsubs_mean_pupil_R_lock_RT_80_100 = nanmean(mean_pupil_R_lock_RT_80_100,3);
allsubs_mean_pupil_Rim_lock_RT_0_20 = nanmean(mean_pupil_Rim_lock_RT_0_20,3);
allsubs_mean_pupil_Rim_lock_RT_20_40 = nanmean(mean_pupil_Rim_lock_RT_20_40,3);
allsubs_mean_pupil_Rim_lock_RT_40_60 = nanmean(mean_pupil_Rim_lock_RT_40_60,3);
allsubs_mean_pupil_Rim_lock_RT_60_80 = nanmean(mean_pupil_Rim_lock_RT_60_80,3);
allsubs_mean_pupil_Rim_lock_RT_80_100 = nanmean(mean_pupil_Rim_lock_RT_80_100,3);
% SEM pupil per subjects
% onset locked
allsubs_sem_pupil_E_lock_onset_0_20 = zeros(1,timelength);
allsubs_sem_pupil_E_lock_onset_20_40 = zeros(1,timelength);
allsubs_sem_pupil_E_lock_onset_40_60 = zeros(1,timelength);
allsubs_sem_pupil_E_lock_onset_60_80 = zeros(1,timelength);
allsubs_sem_pupil_E_lock_onset_80_100 = zeros(1,timelength);
allsubs_sem_pupil_R_lock_onset_0_20 = zeros(1,timelength);
allsubs_sem_pupil_R_lock_onset_20_40 = zeros(1,timelength);
allsubs_sem_pupil_R_lock_onset_40_60 = zeros(1,timelength);
allsubs_sem_pupil_R_lock_onset_60_80 = zeros(1,timelength);
allsubs_sem_pupil_R_lock_onset_80_100 = zeros(1,timelength);
allsubs_sem_pupil_Rim_lock_onset_0_20 = zeros(1,timelength);
allsubs_sem_pupil_Rim_lock_onset_20_40 = zeros(1,timelength);
allsubs_sem_pupil_Rim_lock_onset_40_60 = zeros(1,timelength);
allsubs_sem_pupil_Rim_lock_onset_60_80 = zeros(1,timelength);
allsubs_sem_pupil_Rim_lock_onset_80_100 = zeros(1,timelength);
% RT 1st press (whichtask=0)/RT choice locked (whichtask=1)
allsubs_sem_pupil_E_lock_RT_0_20 = zeros(1,timelength);
allsubs_sem_pupil_E_lock_RT_20_40 = zeros(1,timelength);
allsubs_sem_pupil_E_lock_RT_40_60 = zeros(1,timelength);
allsubs_sem_pupil_E_lock_RT_60_80 = zeros(1,timelength);
allsubs_sem_pupil_E_lock_RT_80_100 = zeros(1,timelength);
allsubs_sem_pupil_R_lock_RT_0_20 = zeros(1,timelength);
allsubs_sem_pupil_R_lock_RT_20_40 = zeros(1,timelength);
allsubs_sem_pupil_R_lock_RT_40_60 = zeros(1,timelength);
allsubs_sem_pupil_R_lock_RT_60_80 = zeros(1,timelength);
allsubs_sem_pupil_R_lock_RT_80_100 = zeros(1,timelength);
allsubs_sem_pupil_Rim_lock_RT_0_20 = zeros(1,timelength);
allsubs_sem_pupil_Rim_lock_RT_20_40 = zeros(1,timelength);
allsubs_sem_pupil_Rim_lock_RT_40_60 = zeros(1,timelength);
allsubs_sem_pupil_Rim_lock_RT_60_80 = zeros(1,timelength);
allsubs_sem_pupil_Rim_lock_RT_80_100 = zeros(1,timelength);
% create new var because grpstats cannot process 3D variables
% onset locked
mp_E_onset_0_20(:,:) = mean_pupil_E_lock_onset_0_20(1,:,:);
mp_E_onset_20_40(:,:) = mean_pupil_E_lock_onset_20_40(1,:,:);
mp_E_onset_40_60(:,:) = mean_pupil_E_lock_onset_40_60(1,:,:);
mp_E_onset_60_80(:,:) = mean_pupil_E_lock_onset_60_80(1,:,:);
mp_E_onset_80_100(:,:) = mean_pupil_E_lock_onset_80_100(1,:,:);
mp_R_onset_0_20(:,:) = mean_pupil_R_lock_onset_0_20(1,:,:);
mp_R_onset_20_40(:,:) = mean_pupil_R_lock_onset_20_40(1,:,:);
mp_R_onset_40_60(:,:) = mean_pupil_R_lock_onset_40_60(1,:,:);
mp_R_onset_60_80(:,:) = mean_pupil_R_lock_onset_60_80(1,:,:);
mp_R_onset_80_100(:,:) = mean_pupil_R_lock_onset_80_100(1,:,:);
mp_Rim_onset_0_20(:,:) = mean_pupil_Rim_lock_onset_0_20(1,:,:);
mp_Rim_onset_20_40(:,:) = mean_pupil_Rim_lock_onset_20_40(1,:,:);
mp_Rim_onset_40_60(:,:) = mean_pupil_Rim_lock_onset_40_60(1,:,:);
mp_Rim_onset_60_80(:,:) = mean_pupil_Rim_lock_onset_60_80(1,:,:);
mp_Rim_onset_80_100(:,:) = mean_pupil_Rim_lock_onset_80_100(1,:,:);
% RT 1st press (whichtask=0)/RT choice locked (whichtask=1)
mp_E_RT_0_20(:,:) = mean_pupil_E_lock_RT_0_20(1,:,:);
mp_E_RT_20_40(:,:) = mean_pupil_E_lock_RT_20_40(1,:,:);
mp_E_RT_40_60(:,:) = mean_pupil_E_lock_RT_40_60(1,:,:);
mp_E_RT_60_80(:,:) = mean_pupil_E_lock_RT_60_80(1,:,:);
mp_E_RT_80_100(:,:) = mean_pupil_E_lock_RT_80_100(1,:,:);
mp_R_RT_0_20(:,:) = mean_pupil_R_lock_RT_0_20(1,:,:);
mp_R_RT_20_40(:,:) = mean_pupil_R_lock_RT_20_40(1,:,:);
mp_R_RT_40_60(:,:) = mean_pupil_R_lock_RT_40_60(1,:,:);
mp_R_RT_60_80(:,:) = mean_pupil_R_lock_RT_60_80(1,:,:);
mp_R_RT_80_100(:,:) = mean_pupil_R_lock_RT_80_100(1,:,:);
mp_Rim_RT_0_20(:,:) = mean_pupil_Rim_lock_RT_0_20(1,:,:);
mp_Rim_RT_20_40(:,:) = mean_pupil_Rim_lock_RT_20_40(1,:,:);
mp_Rim_RT_40_60(:,:) = mean_pupil_Rim_lock_RT_40_60(1,:,:);
mp_Rim_RT_60_80(:,:) = mean_pupil_Rim_lock_RT_60_80(1,:,:);
mp_Rim_RT_80_100(:,:) = mean_pupil_Rim_lock_RT_80_100(1,:,:);

if whichtask == 0
    % pupil mean the subjects
    % RTval
    allsubs_mean_pupil_E_lock_RTval_0_20 = nanmean(mean_pupil_E_lock_RTval_0_20,3);
    allsubs_mean_pupil_E_lock_RTval_20_40 = nanmean(mean_pupil_E_lock_RTval_20_40,3);
    allsubs_mean_pupil_E_lock_RTval_40_60 = nanmean(mean_pupil_E_lock_RTval_40_60,3);
    allsubs_mean_pupil_E_lock_RTval_60_80 = nanmean(mean_pupil_E_lock_RTval_60_80,3);
    allsubs_mean_pupil_E_lock_RTval_80_100 = nanmean(mean_pupil_E_lock_RTval_80_100,3);
    allsubs_mean_pupil_R_lock_RTval_0_20 = nanmean(mean_pupil_R_lock_RTval_0_20,3);
    allsubs_mean_pupil_R_lock_RTval_20_40 = nanmean(mean_pupil_R_lock_RTval_20_40,3);
    allsubs_mean_pupil_R_lock_RTval_40_60 = nanmean(mean_pupil_R_lock_RTval_40_60,3);
    allsubs_mean_pupil_R_lock_RTval_60_80 = nanmean(mean_pupil_R_lock_RTval_60_80,3);
    allsubs_mean_pupil_R_lock_RTval_80_100 = nanmean(mean_pupil_R_lock_RTval_80_100,3);
    allsubs_mean_pupil_Rim_lock_RTval_0_20 = nanmean(mean_pupil_Rim_lock_RTval_0_20,3);
    allsubs_mean_pupil_Rim_lock_RTval_20_40 = nanmean(mean_pupil_Rim_lock_RTval_20_40,3);
    allsubs_mean_pupil_Rim_lock_RTval_40_60 = nanmean(mean_pupil_Rim_lock_RTval_40_60,3);
    allsubs_mean_pupil_Rim_lock_RTval_60_80 = nanmean(mean_pupil_Rim_lock_RTval_60_80,3);
    allsubs_mean_pupil_Rim_lock_RTval_80_100 = nanmean(mean_pupil_Rim_lock_RTval_80_100,3);
    % SEM pupil per subjects
    % RTval
    allsubs_sem_pupil_E_lock_RTval_0_20 = zeros(1,timelength);
    allsubs_sem_pupil_E_lock_RTval_20_40 = zeros(1,timelength);
    allsubs_sem_pupil_E_lock_RTval_40_60 = zeros(1,timelength);
    allsubs_sem_pupil_E_lock_RTval_60_80 = zeros(1,timelength);
    allsubs_sem_pupil_E_lock_RTval_80_100 = zeros(1,timelength);
    allsubs_sem_pupil_R_lock_RTval_0_20 = zeros(1,timelength);
    allsubs_sem_pupil_R_lock_RTval_20_40 = zeros(1,timelength);
    allsubs_sem_pupil_R_lock_RTval_40_60 = zeros(1,timelength);
    allsubs_sem_pupil_R_lock_RTval_60_80 = zeros(1,timelength);
    allsubs_sem_pupil_R_lock_RTval_80_100 = zeros(1,timelength);
    allsubs_sem_pupil_Rim_lock_RTval_0_20 = zeros(1,timelength);
    allsubs_sem_pupil_Rim_lock_RTval_20_40 = zeros(1,timelength);
    allsubs_sem_pupil_Rim_lock_RTval_40_60 = zeros(1,timelength);
    allsubs_sem_pupil_Rim_lock_RTval_60_80 = zeros(1,timelength);
    allsubs_sem_pupil_Rim_lock_RTval_80_100 = zeros(1,timelength);
    % create new var because grpstats cannot process 3D variables
    % RTval
    mp_E_RTval_0_20(:,:) = mean_pupil_E_lock_RTval_0_20(1,:,:);
    mp_E_RTval_20_40(:,:) = mean_pupil_E_lock_RTval_20_40(1,:,:);
    mp_E_RTval_40_60(:,:) = mean_pupil_E_lock_RTval_40_60(1,:,:);
    mp_E_RTval_60_80(:,:) = mean_pupil_E_lock_RTval_60_80(1,:,:);
    mp_E_RTval_80_100(:,:) = mean_pupil_E_lock_RTval_80_100(1,:,:);
    mp_R_RTval_0_20(:,:) = mean_pupil_R_lock_RTval_0_20(1,:,:);
    mp_R_RTval_20_40(:,:) = mean_pupil_R_lock_RTval_20_40(1,:,:);
    mp_R_RTval_40_60(:,:) = mean_pupil_R_lock_RTval_40_60(1,:,:);
    mp_R_RTval_60_80(:,:) = mean_pupil_R_lock_RTval_60_80(1,:,:);
    mp_R_RTval_80_100(:,:) = mean_pupil_R_lock_RTval_80_100(1,:,:);
    mp_Rim_RTval_0_20(:,:) = mean_pupil_Rim_lock_RTval_0_20(1,:,:);
    mp_Rim_RTval_20_40(:,:) = mean_pupil_Rim_lock_RTval_20_40(1,:,:);
    mp_Rim_RTval_40_60(:,:) = mean_pupil_Rim_lock_RTval_40_60(1,:,:);
    mp_Rim_RTval_60_80(:,:) = mean_pupil_Rim_lock_RTval_60_80(1,:,:);
    mp_Rim_RTval_80_100(:,:) = mean_pupil_Rim_lock_RTval_80_100(1,:,:);
end

for time = 1:timelength
    allsubs_sem_pupil_E_lock_onset_0_20(time) = grpstats(mp_E_onset_0_20(time,:), [], 'sem');
    allsubs_sem_pupil_E_lock_onset_20_40(time) = grpstats(mp_E_onset_20_40(time,:), [], 'sem');
    allsubs_sem_pupil_E_lock_onset_40_60(time) = grpstats(mp_E_onset_40_60(time,:), [], 'sem');
    allsubs_sem_pupil_E_lock_onset_60_80(time) = grpstats(mp_E_onset_60_80(time,:), [], 'sem');
    allsubs_sem_pupil_E_lock_onset_80_100(time) = grpstats(mp_E_onset_80_100(time,:), [], 'sem');
    allsubs_sem_pupil_R_lock_onset_0_20(time) = grpstats(mp_R_onset_0_20(time,:), [], 'sem');
    allsubs_sem_pupil_R_lock_onset_20_40(time) = grpstats(mp_R_onset_20_40(time,:), [], 'sem');
    allsubs_sem_pupil_R_lock_onset_40_60(time) = grpstats(mp_R_onset_40_60(time,:), [], 'sem');
    allsubs_sem_pupil_R_lock_onset_60_80(time) = grpstats(mp_R_onset_60_80(time,:), [], 'sem');
    allsubs_sem_pupil_R_lock_onset_80_100(time) = grpstats(mp_R_onset_80_100(time,:), [], 'sem');
    allsubs_sem_pupil_Rim_lock_onset_0_20(time) = grpstats(mp_Rim_onset_0_20(time,:), [], 'sem');
    allsubs_sem_pupil_Rim_lock_onset_20_40(time) = grpstats(mp_Rim_onset_20_40(time,:), [], 'sem');
    allsubs_sem_pupil_Rim_lock_onset_40_60(time) = grpstats(mp_Rim_onset_40_60(time,:), [], 'sem');
    allsubs_sem_pupil_Rim_lock_onset_60_80(time) = grpstats(mp_Rim_onset_60_80(time,:), [], 'sem');
    allsubs_sem_pupil_Rim_lock_onset_80_100(time) = grpstats(mp_Rim_onset_80_100(time,:), [], 'sem');
    % RT
    allsubs_sem_pupil_E_lock_RT_0_20(time) = grpstats(mp_E_RT_0_20(time,:), [], 'sem');
    allsubs_sem_pupil_E_lock_RT_20_40(time) = grpstats(mp_E_RT_20_40(time,:), [], 'sem');
    allsubs_sem_pupil_E_lock_RT_40_60(time) = grpstats(mp_E_RT_40_60(time,:), [], 'sem');
    allsubs_sem_pupil_E_lock_RT_60_80(time) = grpstats(mp_E_RT_60_80(time,:), [], 'sem');
    allsubs_sem_pupil_E_lock_RT_80_100(time) = grpstats(mp_E_RT_80_100(time,:), [], 'sem');
    allsubs_sem_pupil_R_lock_RT_0_20(time) = grpstats(mp_R_RT_0_20(time,:), [], 'sem');
    allsubs_sem_pupil_R_lock_RT_20_40(time) = grpstats(mp_R_RT_20_40(time,:), [], 'sem');
    allsubs_sem_pupil_R_lock_RT_40_60(time) = grpstats(mp_R_RT_40_60(time,:), [], 'sem');
    allsubs_sem_pupil_R_lock_RT_60_80(time) = grpstats(mp_R_RT_60_80(time,:), [], 'sem');
    allsubs_sem_pupil_R_lock_RT_80_100(time) = grpstats(mp_R_RT_80_100(time,:), [], 'sem');
    allsubs_sem_pupil_Rim_lock_RT_0_20(time) = grpstats(mp_Rim_RT_0_20(time,:), [], 'sem');
    allsubs_sem_pupil_Rim_lock_RT_20_40(time) = grpstats(mp_Rim_RT_20_40(time,:), [], 'sem');
    allsubs_sem_pupil_Rim_lock_RT_40_60(time) = grpstats(mp_Rim_RT_40_60(time,:), [], 'sem');
    allsubs_sem_pupil_Rim_lock_RT_60_80(time) = grpstats(mp_Rim_RT_60_80(time,:), [], 'sem');
    allsubs_sem_pupil_Rim_lock_RT_80_100(time) = grpstats(mp_Rim_RT_80_100(time,:), [], 'sem');
    if whichtask == 0
        % RTval
        allsubs_sem_pupil_E_lock_RTval_0_20(time) = grpstats(mp_E_RTval_0_20(time,:), [], 'sem');
        allsubs_sem_pupil_E_lock_RTval_20_40(time) = grpstats(mp_E_RTval_20_40(time,:), [], 'sem');
        allsubs_sem_pupil_E_lock_RTval_40_60(time) = grpstats(mp_E_RTval_40_60(time,:), [], 'sem');
        allsubs_sem_pupil_E_lock_RTval_60_80(time) = grpstats(mp_E_RTval_60_80(time,:), [], 'sem');
        allsubs_sem_pupil_E_lock_RTval_80_100(time) = grpstats(mp_E_RTval_80_100(time,:), [], 'sem');
        allsubs_sem_pupil_R_lock_RTval_0_20(time) = grpstats(mp_R_RTval_0_20(time,:), [], 'sem');
        allsubs_sem_pupil_R_lock_RTval_20_40(time) = grpstats(mp_R_RTval_20_40(time,:), [], 'sem');
        allsubs_sem_pupil_R_lock_RTval_40_60(time) = grpstats(mp_R_RTval_40_60(time,:), [], 'sem');
        allsubs_sem_pupil_R_lock_RTval_60_80(time) = grpstats(mp_R_RTval_60_80(time,:), [], 'sem');
        allsubs_sem_pupil_R_lock_RTval_80_100(time) = grpstats(mp_R_RTval_80_100(time,:), [], 'sem');
        allsubs_sem_pupil_Rim_lock_RTval_0_20(time) = grpstats(mp_Rim_RTval_0_20(time,:), [], 'sem');
        allsubs_sem_pupil_Rim_lock_RTval_20_40(time) = grpstats(mp_Rim_RTval_20_40(time,:), [], 'sem');
        allsubs_sem_pupil_Rim_lock_RTval_40_60(time) = grpstats(mp_Rim_RTval_40_60(time,:), [], 'sem');
        allsubs_sem_pupil_Rim_lock_RTval_60_80(time) = grpstats(mp_Rim_RTval_60_80(time,:), [], 'sem');
        allsubs_sem_pupil_Rim_lock_RTval_80_100(time) = grpstats(mp_Rim_RTval_80_100(time,:), [], 'sem');
    end
end

%% figures
min_pupil = 3000;
max_pupil = 7000;

%% onset locked
min_onset = 0;
max_onset = 2500;
% display mean +- SEM (very noisy)
figure
drawnow;
set(get(handle(gcf),'JavaFrame'),'Maximized',1); % maximize window size
% effort
subplot(1,3,1);
errorbar(allsubs_mean_pupil_E_lock_onset_0_20, allsubs_sem_pupil_E_lock_onset_0_20,'k-');
hold on
errorbar(allsubs_mean_pupil_E_lock_onset_20_40, allsubs_sem_pupil_E_lock_onset_20_40,'g-');
errorbar(allsubs_mean_pupil_E_lock_onset_40_60, allsubs_sem_pupil_E_lock_onset_40_60,'r-');
errorbar(allsubs_mean_pupil_E_lock_onset_60_80, allsubs_sem_pupil_E_lock_onset_60_80,'b-');
errorbar(allsubs_mean_pupil_E_lock_onset_80_100, allsubs_sem_pupil_E_lock_onset_80_100,'m-');
xlim([min_onset max_onset]);
ylim([min_pupil max_pupil]);
xlabel(sprintf(['Time after onset of stimulus (ms) \n',task,' Efforts']));
ylabel('Pupil diameter (micrometers) +/- SEM');
legend('0-20','20-40','40-60','60-80','80-100','Location','northeast');
legend('boxoff')
set(findobj(gca,'type','line'),'linew',2);
set(gca,'fontsize',20,'FontWeight','bold');
set(findall(gcf,'type','text'),'FontSize',20,'fontWeight','bold');
% reward
subplot(1,3,2);
errorbar(allsubs_mean_pupil_R_lock_onset_0_20, allsubs_sem_pupil_R_lock_onset_0_20,'k-');
hold on
errorbar(allsubs_mean_pupil_R_lock_onset_20_40, allsubs_sem_pupil_R_lock_onset_20_40,'g-');
errorbar(allsubs_mean_pupil_R_lock_onset_40_60, allsubs_sem_pupil_R_lock_onset_40_60,'r-');
errorbar(allsubs_mean_pupil_R_lock_onset_60_80, allsubs_sem_pupil_R_lock_onset_60_80,'b-');
errorbar(allsubs_mean_pupil_R_lock_onset_80_100, allsubs_sem_pupil_R_lock_onset_80_100,'m-');
xlim([min_onset max_onset]);
ylim([min_pupil max_pupil]);
xlabel(sprintf(['Time after onset of stimulus (ms) \n ',task,' Reward (text)']));
ylabel('Pupil diameter (micrometers) +/- SEM');
legend('0-20','20-40','40-60','60-80','80-100','Location','northeast');
legend('boxoff')
set(findobj(gca,'type','line'),'linew',2);
set(gca,'fontsize',20,'FontWeight','bold');
set(findall(gcf,'type','text'),'FontSize',20,'fontWeight','bold');
% reward images
subplot(1,3,3);
errorbar(allsubs_mean_pupil_Rim_lock_onset_0_20, allsubs_sem_pupil_Rim_lock_onset_0_20,'k-');
hold on
errorbar(allsubs_mean_pupil_Rim_lock_onset_20_40, allsubs_sem_pupil_Rim_lock_onset_20_40,'g-');
errorbar(allsubs_mean_pupil_Rim_lock_onset_40_60, allsubs_sem_pupil_Rim_lock_onset_40_60,'r-');
errorbar(allsubs_mean_pupil_Rim_lock_onset_60_80, allsubs_sem_pupil_Rim_lock_onset_60_80,'b-');
errorbar(allsubs_mean_pupil_Rim_lock_onset_80_100, allsubs_sem_pupil_Rim_lock_onset_80_100,'m-');
xlim([min_onset max_onset]);
ylim([min_pupil max_pupil]);
xlabel(sprintf(['Time after onset of stimulus (ms) \n ',task,' Reward (images)']));
ylabel('Pupil diameter (micrometers) +/- SEM');
legend('0-20','20-40','40-60','60-80','80-100','Location','northeast');
legend('boxoff')
set(findobj(gca,'type','line'),'linew',2);
set(gca,'fontsize',20,'FontWeight','bold');
set(findall(gcf,'type','text'),'FontSize',20,'fontWeight','bold');
% save image
cd([root,'\behavior_summary']);
img_name = ['1a_eye_pupil_diameter_',task,'_bins_',num2str(NS),'subs_onset_lock.png'];
if exist(img_name,'file') ~= 2
    set(gcf,'PaperPosition',[0 0 1 1]);
    set(gcf,'PaperPositionMode','auto');
    saveas(gcf,img_name);
end

% display only the means
figure
drawnow;
set(get(handle(gcf),'JavaFrame'),'Maximized',1); % maximize window size
% effort
subplot(1,3,1);
plot(allsubs_mean_pupil_E_lock_onset_0_20,'k-');
hold on
plot(allsubs_mean_pupil_E_lock_onset_20_40,'g-');
plot(allsubs_mean_pupil_E_lock_onset_40_60,'r-');
plot(allsubs_mean_pupil_E_lock_onset_60_80,'b-');
plot(allsubs_mean_pupil_E_lock_onset_80_100,'m-');
xlim([min_onset max_onset]);
ylim([min_pupil max_pupil]);
xlabel(sprintf(['Time after onset of stimulus (ms) \n ',task,' Efforts']));
ylabel('Mean Pupil diameter (micrometers)');
legend('0-20','20-40','40-60','60-80','80-100','Location','northeast');
legend('boxoff')
set(findobj(gca,'type','line'),'linew',2);
set(gca,'fontsize',20,'FontWeight','bold');
set(findall(gcf,'type','text'),'FontSize',20,'fontWeight','bold');
% reward
subplot(1,3,2);
plot(allsubs_mean_pupil_R_lock_onset_0_20,'k-');
hold on
plot(allsubs_mean_pupil_R_lock_onset_20_40,'g-');
plot(allsubs_mean_pupil_R_lock_onset_40_60,'r-');
plot(allsubs_mean_pupil_R_lock_onset_60_80,'b-');
plot(allsubs_mean_pupil_R_lock_onset_80_100,'m-');
xlim([min_onset max_onset]);
ylim([min_pupil max_pupil]);
xlabel(sprintf(['Time after onset of stimulus (ms) \n ',task,' Reward (text)']));
ylabel('Mean Pupil diameter (micrometers)');
legend('0-20','20-40','40-60','60-80','80-100','Location','northeast');
legend('boxoff')
set(findobj(gca,'type','line'),'linew',2);
set(gca,'fontsize',20,'FontWeight','bold');
set(findall(gcf,'type','text'),'FontSize',20,'fontWeight','bold');
% reward images
subplot(1,3,3);
plot(allsubs_mean_pupil_Rim_lock_onset_0_20,'k-');
hold on
plot(allsubs_mean_pupil_Rim_lock_onset_20_40,'g-');
plot(allsubs_mean_pupil_Rim_lock_onset_40_60,'r-');
plot(allsubs_mean_pupil_Rim_lock_onset_60_80,'b-');
plot(allsubs_mean_pupil_Rim_lock_onset_80_100,'m-');
xlim([min_onset max_onset]);
ylim([min_pupil max_pupil]);
xlabel(sprintf(['Time after onset of stimulus (ms) \n ',task,' Reward (images)']));
ylabel('Mean Pupil diameter (micrometers)');
legend('0-20','20-40','40-60','60-80','80-100','Location','northeast');
legend('boxoff')
set(findobj(gca,'type','line'),'linew',2);
set(gca,'fontsize',20,'FontWeight','bold');
set(findall(gcf,'type','text'),'FontSize',20,'fontWeight','bold');
% save image
cd([root,'\behavior_summary']);
img_name = ['1a_eye_pupil_diameter_',task,'_bins_mean_only_',num2str(NS),'subs_onset_lock.png'];
if exist(img_name,'file') ~= 2
    set(gcf,'PaperPosition',[0 0 1 1]);
    set(gcf,'PaperPositionMode','auto');
    saveas(gcf,img_name);
end

%% RT (1st press if rating/choice if choice 1D) locked
if whichtask == 0
    xlab_rt = 'Time before First Button Press (ms)';
    dec = 'dec';
elseif whichtask ==1
    xlab_rt = 'Time before Reaction Time (ms)';
    dec = '';
end
min_RT = timelength - 2500;
max_RT = timelength;
% display mean +- SEM (very noisy)
figure
drawnow;
set(get(handle(gcf),'JavaFrame'),'Maximized',1); % maximize window size
% effort
subplot(1,3,1);
errorbar(allsubs_mean_pupil_E_lock_RT_0_20, allsubs_sem_pupil_E_lock_RT_0_20,'k-');
hold on
errorbar(allsubs_mean_pupil_E_lock_RT_20_40, allsubs_sem_pupil_E_lock_RT_20_40,'g-');
errorbar(allsubs_mean_pupil_E_lock_RT_40_60, allsubs_sem_pupil_E_lock_RT_40_60,'r-');
errorbar(allsubs_mean_pupil_E_lock_RT_60_80, allsubs_sem_pupil_E_lock_RT_60_80,'b-');
errorbar(allsubs_mean_pupil_E_lock_RT_80_100, allsubs_sem_pupil_E_lock_RT_80_100,'m-');
xlim([min_RT max_RT]);
ylim([min_pupil max_pupil]);
xlabel(sprintf([xlab_rt, ' \n ',task,' Efforts']));
ylabel('Pupil diameter (micrometers) +/- SEM');
legend('0-20','20-40','40-60','60-80','80-100','Location','northwest');
legend('boxoff')
set(findobj(gca,'type','line'),'linew',2);
set(gca,'fontsize',20,'FontWeight','bold');
set(findall(gcf,'type','text'),'FontSize',20,'fontWeight','bold');
% reward
subplot(1,3,2);
errorbar(allsubs_mean_pupil_R_lock_RT_0_20, allsubs_sem_pupil_R_lock_RT_0_20,'k-');
hold on
errorbar(allsubs_mean_pupil_R_lock_RT_20_40, allsubs_sem_pupil_R_lock_RT_20_40,'g-');
errorbar(allsubs_mean_pupil_R_lock_RT_40_60, allsubs_sem_pupil_R_lock_RT_40_60,'r-');
errorbar(allsubs_mean_pupil_R_lock_RT_60_80, allsubs_sem_pupil_R_lock_RT_60_80,'b-');
errorbar(allsubs_mean_pupil_R_lock_RT_80_100, allsubs_sem_pupil_R_lock_RT_80_100,'m-');
xlim([min_RT max_RT]);
ylim([min_pupil max_pupil]);
xlabel(sprintf([xlab_rt, ' \n ',task,' Reward (text)']));
ylabel('Pupil diameter (micrometers) +/- SEM');
legend('0-20','20-40','40-60','60-80','80-100','Location','northwest');
legend('boxoff')
set(findobj(gca,'type','line'),'linew',2);
set(gca,'fontsize',20,'FontWeight','bold');
set(findall(gcf,'type','text'),'FontSize',20,'fontWeight','bold');
% reward images
subplot(1,3,3);
errorbar(allsubs_mean_pupil_Rim_lock_RT_0_20, allsubs_sem_pupil_Rim_lock_RT_0_20,'k-');
hold on
errorbar(allsubs_mean_pupil_Rim_lock_RT_20_40, allsubs_sem_pupil_Rim_lock_RT_20_40,'g-');
errorbar(allsubs_mean_pupil_Rim_lock_RT_40_60, allsubs_sem_pupil_Rim_lock_RT_40_60,'r-');
errorbar(allsubs_mean_pupil_Rim_lock_RT_60_80, allsubs_sem_pupil_Rim_lock_RT_60_80,'b-');
errorbar(allsubs_mean_pupil_Rim_lock_RT_80_100, allsubs_sem_pupil_Rim_lock_RT_80_100,'m-');
xlim([min_RT max_RT]);
ylim([min_pupil max_pupil]);
xlabel(sprintf([xlab_rt, ' \n ',task,' Reward (images)']));
ylabel('Pupil diameter (micrometers) +/- SEM');
legend('0-20','20-40','40-60','60-80','80-100','Location','northwest');
legend('boxoff')
set(findobj(gca,'type','line'),'linew',2);
set(gca,'fontsize',20,'FontWeight','bold');
set(findall(gcf,'type','text'),'FontSize',20,'fontWeight','bold');
% save image
cd([root,'\behavior_summary']);
img_name = ['1a_eye_pupil_diameter_',task,'_bins_',num2str(NS),'subs_RT',dec,'_lock.png'];
if exist(img_name,'file') ~= 2
    set(gcf,'PaperPosition',[0 0 1 1]);
    set(gcf,'PaperPositionMode','auto');
    saveas(gcf,img_name);
end

% display only the means
figure
drawnow;
set(get(handle(gcf),'JavaFrame'),'Maximized',1); % maximize window size
% effort
subplot(1,3,1);
plot(allsubs_mean_pupil_E_lock_RT_0_20,'k-');
hold on
plot(allsubs_mean_pupil_E_lock_RT_20_40,'g-');
plot(allsubs_mean_pupil_E_lock_RT_40_60,'r-');
plot(allsubs_mean_pupil_E_lock_RT_60_80,'b-');
plot(allsubs_mean_pupil_E_lock_RT_80_100,'m-');
xlim([min_RT max_RT]);
ylim([min_pupil max_pupil]);
xlabel(sprintf([xlab_rt, ' \n ',task,' Efforts']));
ylabel('Mean Pupil diameter (micrometers)');
legend('0-20','20-40','40-60','60-80','80-100','Location','northwest');
legend('boxoff')
set(findobj(gca,'type','line'),'linew',2);
set(gca,'fontsize',20,'FontWeight','bold');
set(findall(gcf,'type','text'),'FontSize',20,'fontWeight','bold');
% reward
subplot(1,3,2);
plot(allsubs_mean_pupil_R_lock_RT_0_20,'k-');
hold on
plot(allsubs_mean_pupil_R_lock_RT_20_40,'g-');
plot(allsubs_mean_pupil_R_lock_RT_40_60,'r-');
plot(allsubs_mean_pupil_R_lock_RT_60_80,'b-');
plot(allsubs_mean_pupil_R_lock_RT_80_100,'m-');
xlim([min_RT max_RT]);
ylim([min_pupil max_pupil]);
xlabel(sprintf([xlab_rt, ' \n ',task,' Reward (text)']));
ylabel('Mean Pupil diameter (micrometers)');
legend('0-20','20-40','40-60','60-80','80-100','Location','northwest');
legend('boxoff')
set(findobj(gca,'type','line'),'linew',2);
set(gca,'fontsize',20,'FontWeight','bold');
set(findall(gcf,'type','text'),'FontSize',20,'fontWeight','bold');
% reward images
subplot(1,3,3);
plot(allsubs_mean_pupil_Rim_lock_RT_0_20,'k-');
hold on
plot(allsubs_mean_pupil_Rim_lock_RT_20_40,'g-');
plot(allsubs_mean_pupil_Rim_lock_RT_40_60,'r-');
plot(allsubs_mean_pupil_Rim_lock_RT_60_80,'b-');
plot(allsubs_mean_pupil_Rim_lock_RT_80_100,'m-');
xlim([min_RT max_RT]);
ylim([min_pupil max_pupil]);
xlabel(sprintf([xlab_rt, ' \n ',task,' Reward (images)']));
ylabel('Mean Pupil diameter (micrometers)');
legend('0-20','20-40','40-60','60-80','80-100','Location','northwest');
legend('boxoff')
set(findobj(gca,'type','line'),'linew',2);
set(gca,'fontsize',20,'FontWeight','bold');
set(findall(gcf,'type','text'),'FontSize',20,'fontWeight','bold');
% save image
cd([root,'\behavior_summary']);
img_name = ['1a_eye_pupil_diameter_',task,'_bins_mean_only_',num2str(NS),'subs_RT',dec,'_lock.png'];
if exist(img_name,'file') ~= 2
    set(gcf,'PaperPosition',[0 0 1 1]);
    set(gcf,'PaperPositionMode','auto');
    saveas(gcf,img_name);
end

%% RT validation for rating task
if whichtask == 0
    min_RTval = timelength - 2500;
    max_RTval = timelength;
    % display mean +- SEM (very noisy)
    figure
    drawnow;
    set(get(handle(gcf),'JavaFrame'),'Maximized',1); % maximize window size
    % effort
    subplot(1,3,1);
    errorbar(allsubs_mean_pupil_E_lock_RTval_0_20, allsubs_sem_pupil_E_lock_RTval_0_20,'k-');
    hold on
    errorbar(allsubs_mean_pupil_E_lock_RTval_20_40, allsubs_sem_pupil_E_lock_RTval_20_40,'g-');
    errorbar(allsubs_mean_pupil_E_lock_RTval_40_60, allsubs_sem_pupil_E_lock_RTval_40_60,'r-');
    errorbar(allsubs_mean_pupil_E_lock_RTval_60_80, allsubs_sem_pupil_E_lock_RTval_60_80,'b-');
    errorbar(allsubs_mean_pupil_E_lock_RTval_80_100, allsubs_sem_pupil_E_lock_RTval_80_100,'m-');
    xlim([min_RTval max_RTval]);
    ylim([min_pupil max_pupil]);
    xlabel(sprintf(['Time before validation (ms) \n ',task,' Efforts']));
    ylabel('Pupil diameter (micrometers) +/- SEM');
    legend('0-20','20-40','40-60','60-80','80-100','Location','northwest');
    legend('boxoff')
    set(findobj(gca,'type','line'),'linew',2);
    set(gca,'fontsize',20,'FontWeight','bold');
    set(findall(gcf,'type','text'),'FontSize',20,'fontWeight','bold');
    % reward
    subplot(1,3,2);
    errorbar(allsubs_mean_pupil_R_lock_RTval_0_20, allsubs_sem_pupil_R_lock_RTval_0_20,'k-');
    hold on
    errorbar(allsubs_mean_pupil_R_lock_RTval_20_40, allsubs_sem_pupil_R_lock_RTval_20_40,'g-');
    errorbar(allsubs_mean_pupil_R_lock_RTval_40_60, allsubs_sem_pupil_R_lock_RTval_40_60,'r-');
    errorbar(allsubs_mean_pupil_R_lock_RTval_60_80, allsubs_sem_pupil_R_lock_RTval_60_80,'b-');
    errorbar(allsubs_mean_pupil_R_lock_RTval_80_100, allsubs_sem_pupil_R_lock_RTval_80_100,'m-');
    xlim([min_RTval max_RTval]);
    ylim([min_pupil max_pupil]);
    xlabel(sprintf(['Time before validation (ms) \n ',task,' Reward (text)']));
    ylabel('Pupil diameter (micrometers) +/- SEM');
    legend('0-20','20-40','40-60','60-80','80-100','Location','northwest');
    legend('boxoff')
    set(findobj(gca,'type','line'),'linew',2);
    set(gca,'fontsize',20,'FontWeight','bold');
    set(findall(gcf,'type','text'),'FontSize',20,'fontWeight','bold');
    % reward images
    subplot(1,3,3);
    errorbar(allsubs_mean_pupil_Rim_lock_RTval_0_20, allsubs_sem_pupil_Rim_lock_RTval_0_20,'k-');
    hold on
    errorbar(allsubs_mean_pupil_Rim_lock_RTval_20_40, allsubs_sem_pupil_Rim_lock_RTval_20_40,'g-');
    errorbar(allsubs_mean_pupil_Rim_lock_RTval_40_60, allsubs_sem_pupil_Rim_lock_RTval_40_60,'r-');
    errorbar(allsubs_mean_pupil_Rim_lock_RTval_60_80, allsubs_sem_pupil_Rim_lock_RTval_60_80,'b-');
    errorbar(allsubs_mean_pupil_Rim_lock_RTval_80_100, allsubs_sem_pupil_Rim_lock_RTval_80_100,'m-');
    xlim([min_RTval max_RTval]);
    ylim([min_pupil max_pupil]);
    xlabel(sprintf(['Time before validation (ms) \n ',task,' Reward (images)']));
    ylabel('Pupil diameter (micrometers) +/- SEM');
    legend('0-20','20-40','40-60','60-80','80-100','Location','northwest');
    legend('boxoff')
    set(findobj(gca,'type','line'),'linew',2);
    set(gca,'fontsize',20,'FontWeight','bold');
    set(findall(gcf,'type','text'),'FontSize',20,'fontWeight','bold');
    % save image
    cd([root,'\behavior_summary']);
    img_name = ['1a_eye_pupil_diameter_',task,'_bins_',num2str(NS),'subs_RTval_lock.png'];
    if exist(img_name,'file') ~= 2
        set(gcf,'PaperPosition',[0 0 1 1]);
        set(gcf,'PaperPositionMode','auto');
        saveas(gcf,img_name);
    end
    
    % display only the means
    figure
    drawnow;
    set(get(handle(gcf),'JavaFrame'),'Maximized',1); % maximize window size
    % effort
    subplot(1,3,1);
    plot(allsubs_mean_pupil_E_lock_RTval_0_20,'k-');
    hold on
    plot(allsubs_mean_pupil_E_lock_RTval_20_40,'g-');
    plot(allsubs_mean_pupil_E_lock_RTval_40_60,'r-');
    plot(allsubs_mean_pupil_E_lock_RTval_60_80,'b-');
    plot(allsubs_mean_pupil_E_lock_RTval_80_100,'m-');
    xlim([min_RTval max_RTval]);
    ylim([min_pupil max_pupil]);
    xlabel(sprintf(['Time before validation (ms) \n ',task,' Efforts']));
    ylabel('Mean Pupil diameter (micrometers)');
    legend('0-20','20-40','40-60','60-80','80-100','Location','northwest');
    legend('boxoff')
    set(findobj(gca,'type','line'),'linew',2);
    set(gca,'fontsize',20,'FontWeight','bold');
    set(findall(gcf,'type','text'),'FontSize',20,'fontWeight','bold');
    % reward
    subplot(1,3,2);
    plot(allsubs_mean_pupil_R_lock_RTval_0_20,'k-');
    hold on
    plot(allsubs_mean_pupil_R_lock_RTval_20_40,'g-');
    plot(allsubs_mean_pupil_R_lock_RTval_40_60,'r-');
    plot(allsubs_mean_pupil_R_lock_RTval_60_80,'b-');
    plot(allsubs_mean_pupil_R_lock_RTval_80_100,'m-');
    xlim([min_RTval max_RTval]);
    ylim([min_pupil max_pupil]);
    xlabel(sprintf(['Time before validation (ms) \n ',task,' Reward (text)']));
    ylabel('Mean Pupil diameter (micrometers)');
    legend('0-20','20-40','40-60','60-80','80-100','Location','northwest');
    legend('boxoff')
    set(findobj(gca,'type','line'),'linew',2);
    set(gca,'fontsize',20,'FontWeight','bold');
    set(findall(gcf,'type','text'),'FontSize',20,'fontWeight','bold');
    % reward images
    subplot(1,3,3);
    plot(allsubs_mean_pupil_Rim_lock_RTval_0_20,'k-');
    hold on
    plot(allsubs_mean_pupil_Rim_lock_RTval_20_40,'g-');
    plot(allsubs_mean_pupil_Rim_lock_RTval_40_60,'r-');
    plot(allsubs_mean_pupil_Rim_lock_RTval_60_80,'b-');
    plot(allsubs_mean_pupil_Rim_lock_RTval_80_100,'m-');
    xlim([min_RTval max_RTval]);
    ylim([min_pupil max_pupil]);
    xlabel(sprintf(['Time before validation (ms) \n ',task,' Reward (images)']));
    ylabel('Mean Pupil diameter (micrometers)');
    legend('0-20','20-40','40-60','60-80','80-100','Location','northwest');
    legend('boxoff')
    set(findobj(gca,'type','line'),'linew',2);
    set(gca,'fontsize',20,'FontWeight','bold');
    set(findall(gcf,'type','text'),'FontSize',20,'fontWeight','bold');
    % save image
    cd([root,'\behavior_summary']);
    img_name = ['1a_eye_pupil_diameter_',task,'_bins_mean_only_',num2str(NS),'subs_RTval_lock.png'];
    if exist(img_name,'file') ~= 2
        set(gcf,'PaperPosition',[0 0 1 1]);
        set(gcf,'PaperPositionMode','auto');
        saveas(gcf,img_name);
    end
end
% end