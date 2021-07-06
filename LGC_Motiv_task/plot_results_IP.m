% Script to extract various values from pilots and plot them.

%% clean workspace before startingclearvars;
clearvars;
close all;
clc;

%% working directories
% launch within the folder where scripts are stored or will not work
cd ..
main_folder                 = [pwd filesep]; % you have to be sure that you are in the correct path when you launch the script
main_task_folder            = [main_folder, 'LGC_Motiv_task' filesep];
results_folder              = [main_folder, 'LGC_Motiv_results' filesep 'pilots_v1_IP_Nback1' filesep];

main_task_folder = 'D:\LGC_motiv\LGC_Motiv_task\';

% add personal functions (needed for PTB opening at least)
addpath(genpath(results_folder));

cd(results_folder)

%% find the number of participants and initialize parameters
fstruct = dir('IP_pilot_data*_sub_*.mat')
name_tmp = fstruct(1).name;


nb_pilots = str2double(name_tmp(21));
% nb_pilots = 1;

pilot_ID = [1 2 3 4];

% % go back to folder with scripts
% cd(main_task_folder);

%% extract relevant features
% IP matrice has as columns : 
%  Effort lvl 1/2, Physical/Mental, Reward/Punishment, iteration 1/2
%  1      2      3      4      5     6      7       8
% E1PR1,E1PR2, E1MR1, E1MR2, E1PP1, E1PP2, E1MP1, E1MP2
%  9      10     11     12     13     14     15     16
% E2PR1,E2PR2, E2MR1, E2MR2, E2PP1, E2PP2, E2MP1, E2MP2

for i_pilot = 1:nb_pilots
    loading_file_nm_tmp = ['*',num2str(pilot_ID(i_pilot)),'.mat'];
    load(dir(loading_file_nm_tmp).name)
    delta_IP(i_pilot,1) = -1.5 + all.physical.EffortLvl_1.session_nb1.repeat_nb1.perfSummary.IP;
    delta_IP(i_pilot,2) = -1.5 + all.physical.EffortLvl_1.session_nb1.repeat_nb2.perfSummary.IP;
    delta_IP(i_pilot,3) = -1.5 + all.mental.EffortLvl_1.session_nb1.repeat_nb1.perfSummary.IP;
    delta_IP(i_pilot,4) = -1.5 + all.mental.EffortLvl_1.session_nb1.repeat_nb2.perfSummary.IP;
    delta_IP(i_pilot,5) = 1.5 - all.physical.EffortLvl_1.session_nb2.repeat_nb1.perfSummary.IP;
    delta_IP(i_pilot,6) = 1.5 - all.physical.EffortLvl_1.session_nb2.repeat_nb2.perfSummary.IP;
    delta_IP(i_pilot,7) = 1.5 - all.mental.EffortLvl_1.session_nb2.repeat_nb1.perfSummary.IP;
    delta_IP(i_pilot,8) = 1.5 - all.mental.EffortLvl_1.session_nb2.repeat_nb2.perfSummary.IP;
    delta_IP(i_pilot,9) = -1.5 + all.physical.EffortLvl_2.session_nb1.repeat_nb1.perfSummary.IP;
    delta_IP(i_pilot,10) = -1.5 + all.physical.EffortLvl_2.session_nb1.repeat_nb2.perfSummary.IP;
    delta_IP(i_pilot,11) = -1.5 + all.mental.EffortLvl_2.session_nb1.repeat_nb1.perfSummary.IP;
    delta_IP(i_pilot,12) = -1.5 + all.mental.EffortLvl_2.session_nb1.repeat_nb2.perfSummary.IP;
    delta_IP(i_pilot,13) = 1.5 - all.physical.EffortLvl_2.session_nb2.repeat_nb1.perfSummary.IP;
    delta_IP(i_pilot,14) = 1.5 - all.physical.EffortLvl_2.session_nb2.repeat_nb2.perfSummary.IP;
    delta_IP(i_pilot,15) = 1.5 - all.mental.EffortLvl_2.session_nb2.repeat_nb1.perfSummary.IP;
    delta_IP(i_pilot,16) = 1.5 - all.mental.EffortLvl_2.session_nb2.repeat_nb2.perfSummary.IP;
    if exist('initial_MVC')
        delta_MVC(i_pilot) = (initial_MVC.MVC - last_MVC.MVC)/ initial_MVC.MVC * 100;
    end

    if exist('t_min_calib')
        delta_MVM(i_pilot) = t_min_lastCalib - t_min_calib;
    end
       
    E1_P_R(i_pilot) = mean(delta_IP(i_pilot,1:2),2);
    E1_M_R(i_pilot) = mean(delta_IP(i_pilot,3:4),2);
    E1_P_P(i_pilot) = mean(delta_IP(i_pilot,5:6),2);
    E1_M_P(i_pilot) = mean(delta_IP(i_pilot,7:8),2);
    E2_P_R(i_pilot) = mean(delta_IP(i_pilot,9:10),2);
    E2_M_R(i_pilot) = mean(delta_IP(i_pilot,11:12),2);
    E2_P_P(i_pilot) = mean(delta_IP(i_pilot,13:14),2);
    E2_M_P(i_pilot) = mean(delta_IP(i_pilot,15:16),2);
end


%% plots
figure()
bar(delta_MVC)

figure()
bar(delta_MVM)

figure()

bar([E1_P_R; E1_M_R; E1_P_P; E1_M_P;E2_P_R; E2_M_R; E2_P_P; E2_M_P])
set(gca,'XTick',[1 2 3 4 5 6 7 8])
set(gca,'XTickLabel',{'E1PR','E1MR','E1PP','E1MP','E2PR','E2MR','E2PP','E2MP'})



% % go back to folder with scripts
cd(main_task_folder);

