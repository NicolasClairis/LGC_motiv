% Script to extract various values from pilots and plot them.

%% clean workspace before startingclearvars;
clearvars;
close all;
clc;

%% working directories
% launch within the folder where scripts are stored or will not work
cd ..
main_folder                 = [pwd filesep]; % you have to be sure that you are in the correct path when you launch the script
main_task_folder_analysis            = [main_folder, 'LGC_Motiv_task' filesep];
folder_to_open              = ['pilots_v6_IP_Nback2_NOtaskSwitching_NOriskRepeatAfterFail'];
results_folder              = [main_folder, 'LGC_Motiv_results' filesep folder_to_open filesep];

main_task_folder_analysis = 'D:\LGC_motiv\LGC_Motiv_task\';

% add personal functions (needed for PTB opening at least)
addpath(genpath(results_folder));

cd(results_folder)

%% find the number of participants and initialize parameters
fstruct = dir('IP_pilot_data*_sub_*.mat');

nb_files_per_pilot = 19;
% find the number of participants in the folder
nb_pilots = round(size(fstruct,1)/nb_files_per_pilot);
% pilot_ID = linspace(1,nb_pilots,nb_pilots);

% find the ID of pilots
for i = 1:length(fstruct)
    name_tmp = fstruct(i).name;
    if name_tmp(22) ~= '1' || name_tmp(22) ~= '2' || name_tmp(22) ~= '3' || name_tmp(22) ~= '4' || name_tmp(22) ~= '5' || name_tmp(22) ~= '6' || name_tmp(22) ~= '7' || name_tmp(22) ~= '8' || name_tmp(22) ~= '9'
    ID_tmp(i) = str2double(name_tmp(21:22));
    else
        ID_tmp(i) = str2double(name_tmp(21));
    end
end
pilot_ID = unique(ID_tmp);
% pilot_ID = [1 2 3 4];

%% extract relevant features
% IP matrice has as columns :
%  Effort lvl 1/2, Physical/Mental, Reward/Punishment, iteration 1/2
%  1      2      3      4      5     6      7       8
% E1PR1,E1PR2, E1MR1, E1MR2, E1PP1, E1PP2, E1MP1, E1MP2
%  9      10     11     12     13     14     15     16
% E2PR1,E2PR2, E2MR1, E2MR2, E2PP1, E2PP2, E2MP1, E2MP2

for i_pilot = 1:nb_pilots
    loading_file_nm_tmp = ['*_',num2str(pilot_ID(i_pilot)),'.mat'];
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
    if exist('initial_MVC') && exist('last_MVC') 
        delta_MVC(i_pilot) = (last_MVC.MVC - initial_MVC.MVC)/ initial_MVC.MVC * 100;
        init_MVC(i_pilot) = initial_MVC.MVC;
    end
    
    if exist('t_min_calib') && exist('t_min_lastCalib')
        delta_MVM(i_pilot) = t_min_lastCalib - t_min_calib;
        MVM(i_pilot) = t_min_calib;
    end
    
    for i_mean = 1:length(delta_IP)/2
        mean_IP(i_pilot, i_mean) = mean(delta_IP(i_pilot,(i_mean-1)*2 + 1:i_mean*2),2);
        std_IP(i_pilot, i_mean) = std(delta_IP(i_pilot,(i_mean-1)*2 + 1:i_mean*2));
       
    end
end


%% plots
figure()
bar(delta_MVC)

figure()
bar(delta_MVM)

figure()

b = bar(mean_IP', 'grouped');
hold on
set(gca,'XTick',[1 2 3 4 5 6 7 8])
set(gca,'XTickLabel',{'E1PR','E1MR','E1PP','E1MP','E2PR','E2MR','E2PP','E2MP'})

% Calculate the number of groups and number of bars in each group
[ngroups,nbars] = size(mean_IP);
% Get the x coordinate of the bars
x = nan(ngroups,nbars);
for i = 1:ngroups
    x(i,:) = b(i).XEndPoints;
end
errorbar(x, mean_IP,std_IP,'k','linestyle','none')
er.Color = [0 0 0];
er.LineStyle = 'none';
hold off
hold on
set(gca,'XTick',[1 2 3 4 5 6 7 8])
set(gca,'XTickLabel',{'E1-Ph-Re','E1-Me-Re','E1-Ph-Pu','E1-Me-Pu','E2-Ph-Re','E2-Me-Re','E2-Ph-Pu','E2-Me-Pu'})
ylabel('Indifference Point (IP) CHF')
title('E1/E2 low/high effort level, Ph/Me physical mental, Re/Pu reward punishment')



% % go back to folder with scripts
cd(main_task_folder_analysis);

