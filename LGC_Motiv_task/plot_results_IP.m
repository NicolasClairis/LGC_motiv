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
    
    if i_pilot > 3 && i_pilot ~= 10
        for q = 1:90
            totalTime(i_pilot,q) = all.mental.learning.extendedLearning.(strcat('trial_',num2str(q))).totalTime_success;
            errorsMade(i_pilot,q) = all.mental.learning.extendedLearning.(strcat('trial_',num2str(q))).n_errorsMade;
        end
            for m = 1:18
        meanTotalTime(m) = mean(totalTime(i_pilot,1 + (m-1)*5:m*5));
    end
    end
    

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
        end_MVC(i_pilot) = last_MVC.MVC;
    end
    
    if exist('t_min_calib') && exist('t_min_lastCalib')
        delta_MVM(i_pilot) = t_min_lastCalib - t_min_calib;
        init_MVM(i_pilot) = t_min_calib;
        end_MVM(i_pilot) = t_min_lastCalib;
    end
    
    for i_mean = 1:length(delta_IP)/2
        mean_IP(i_pilot, i_mean) = mean(delta_IP(i_pilot,(i_mean-1)*2 + 1:i_mean*2),2);
        std_IP(i_pilot, i_mean) = std(delta_IP(i_pilot,(i_mean-1)*2 + 1:i_mean*2));
        %  1       2       3       4       5      6       7        8
        % E1PhRe, E1MeRe, E1PhPu, E1MePu, E2PhRe, E2MeRe, E2PhPu, E2MePu
    end
        % prepare data for the difference in reward and punishment for low effort
    t_test_re_pu_E1(i_pilot,1) = mean(mean_IP(i_pilot,1:2),2);
    t_test_re_pu_E1(i_pilot,2) = mean(mean_IP(i_pilot,3:4),2);
    
    % prepare data for the difference in physical and mental low efforts
    t_test_ph_me_E1(i_pilot,1) = mean(mean_IP(i_pilot,[1,3]),2);
    t_test_ph_me_E1(i_pilot,2) = mean(mean_IP(i_pilot,[2,4]),2);
    
    % prepare data for the difference in reward and punishment for high effort
    t_test_re_pu_E2(i_pilot,1) = mean(mean_IP(i_pilot,5:6),2);
    t_test_re_pu_E2(i_pilot,2) = mean(mean_IP(i_pilot,7:8),2);
    
    % prepare data for the difference in physical and mental high effort
    t_test_ph_me_E2(i_pilot,1) = mean(mean_IP(i_pilot,[5,7]),2);
    t_test_ph_me_E2(i_pilot,2) = mean(mean_IP(i_pilot,[6,8]),2);
    
    % prepare data for the difference between low and high efforts, indep of other conditions
    t_test_E2(i_pilot,1) = mean(mean_IP(i_pilot,1:4),2);
    t_test_E2(i_pilot,2) = mean(mean_IP(i_pilot,5:8),2);
    
    % prepare data for the difference in reward and punishment, indep of other conditions
    t_test_re_pu(i_pilot,1) = mean(mean_IP(i_pilot,[1,2,5,6]),2);
    t_test_re_pu(i_pilot,2) = mean(mean_IP(i_pilot,[3,4,7,8]),2);
    
    % prepare data for the difference in physical and mental effort, indep of other conditions
    t_test_ph_me(i_pilot,1) = mean(mean_IP(i_pilot,[1,3,5,7]),2);
    t_test_ph_me(i_pilot,2) = mean(mean_IP(i_pilot,[2,4,6,8]),2);
%     % prepare data for the difference in reward and punishment for low effort
%     t_test_re_pu_E1(i_pilot,1) = mean(delta_IP(i_pilot,1:2),2);
%     t_test_re_pu_E1(i_pilot,2) = mean(delta_IP(i_pilot,3:4),2);
%     
%     % prepare data for the difference in physical and mental low efforts
%     t_test_ph_me_E1(i_pilot,1) = mean(delta_IP(i_pilot,[1,3]),2);
%     t_test_ph_me_E1(i_pilot,2) = mean(delta_IP(i_pilot,[2,4]),2);
%     
%     % prepare data for the difference in reward and punishment for high effort
%     t_test_re_pu_E2(i_pilot,1) = mean(delta_IP(i_pilot,5:6),2);
%     t_test_re_pu_E2(i_pilot,2) = mean(delta_IP(i_pilot,7:8),2);
%     
%     % prepare data for the difference in physical and mental high effort
%     t_test_ph_me_E2(i_pilot,1) = mean(delta_IP(i_pilot,[5,7]),2);
%     t_test_ph_me_E2(i_pilot,2) = mean(delta_IP(i_pilot,[6,8]),2);
%     
%     % prepare data for the difference between low and high efforts, indep of other conditions
%     t_test_E2(i_pilot,1) = mean(delta_IP(i_pilot,1:4),2);
%     t_test_E2(i_pilot,2) = mean(delta_IP(i_pilot,5:8),2);
%     
%     % prepare data for the difference in reward and punishment, indep of other conditions
%     t_test_re_pu(i_pilot,1) = mean(delta_IP(i_pilot,[1,2,5,6]),2);
%     t_test_re_pu(i_pilot,2) = mean(delta_IP(i_pilot,[3,4,7,8]),2);
%     
%     % prepare data for the difference in physical and mental effort, indep of other conditions
%     t_test_ph_me(i_pilot,1) = mean(delta_IP(i_pilot,[1,3,5,7]),2);
%     t_test_ph_me(i_pilot,2) = mean(delta_IP(i_pilot,[2,4,6,8]),2);

    end
    

    %% remove outlier
    i_outlier = 9;
    remove_outlier = false;
    if remove_outlier == true
    delta_IP(i_outlier,:) = [];
    delta_MVC(i_outlier) = [];
    delta_MVM(i_outlier) = [];
    mean_IP(i_outlier,:) = [];
    std_IP(i_outlier,:) = [];
    init_MVC(i_outlier) = [];
    end_MVC(i_outlier) = [];
    init_MVM(i_outlier) = [];
    end_MVM(i_outlier) = [];
    

    t_test_re_pu_E1(i_outlier,:) = [];
    t_test_ph_me_E1(i_outlier,:) = [];
    t_test_re_pu_E2(i_outlier,:) = [];
    t_test_ph_me_E2(i_outlier,:) = [];
    t_test_E2(i_outlier,:) = [];
    t_test_re_pu(i_outlier,:) = [];
    t_test_ph_me(i_outlier,:) = [];
    end
%% t-test on our pilots. between conditions and calibrations
figure()
% test with a t test the difference in reward and punishment for low effort
[h_re_pu_E1, p_re_pu_E1] = ttest(t_test_re_pu_E1(:,1),t_test_re_pu_E1(:,2))

% test with a t test the difference in physical and mental low efforts
[h_ph_me_E1, p_ph_me_E1] = ttest(t_test_re_pu_E1(:,1),t_test_re_pu_E1(:,2))

% test with a t test the difference in reward and punishment for high effort
[h_re_pu_E2, p_re_pu_E2] = ttest(t_test_re_pu_E2(:,1),t_test_re_pu_E2(:,2))
% test with a t test the difference in physical and mental high effort
[h_ph_me_E2, p_ph_me_E2] = ttest(t_test_re_pu_E2(:,1),t_test_re_pu_E2(:,2))

% test with a t test the difference in reward and punishment, indep of other conditions
[h_re_pu,p_re_pu] = ttest(t_test_re_pu(:,1),t_test_re_pu(:,2))
subplot(1,3,1)
boxplot([t_test_re_pu(:,1) t_test_re_pu(:,2)],'Labels',{'Reward','Punishment'})
   ax = gca; 
   ax.FontSize = 16;
   ylabel('Indifference point (IP)')
% test with a t test the difference in physical and mental effort, indep of other conditions
[h_ph_me,p_ph_me] = ttest(t_test_ph_me(:,1),t_test_ph_me(:,2))
subplot(1,3,2)
boxplot([t_test_ph_me(:,1) t_test_ph_me(:,2)],'Labels',{'Physical','Mental'})
   ax = gca; 
   ax.FontSize = 16;
      ylabel('Indifference point (IP)')
% test with a t test the difference between low and high efforts, indep of other conditions
[h_E2, p_E2] = ttest(t_test_E2(:,1),t_test_E2(:,2))
subplot(1,3,3)
boxplot([t_test_E2(:,1) t_test_E2(:,2)],'Labels',{'Low effort','High effort'})
   ax = gca; 
   ax.FontSize = 16;
      ylabel('Indifference point (IP)')
      
figure()
% test with t test if MVC is reduced in participants
[h_MVC, p_MVC] = ttest(init_MVC,end_MVC)
subplot(1,2,1)
boxplot([init_MVC; end_MVC]','Labels',{'Initial MVC','Final MVC'})
   ax = gca; 
   ax.FontSize = 16;
      ylabel('Output (Volt)')

% test with t test if MVM is better in participants
[h_MVM, p_MVM] = ttest(init_MVM,end_MVM)
subplot(1,2,2)
boxplot([init_MVM; end_MVM]','Labels',{'Initial MVM','Final MVM'})
   ax = gca; 
   ax.FontSize = 16;
      ylabel('Time to finish (sec)')
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


%% Import fmax_theoritical
opts = spreadsheetImportOptions("NumVariables", 4);

% Specify sheet and range
opts.Sheet = "Sheet1";
opts.DataRange = "C2:F17";

% Specify column names and types
opts.VariableNames = ["Pliantrieur", "Plipostriur", "Longueurdelavantbras", "Circonfrence"];
opts.VariableTypes = ["double", "double", "double", "double"];

% Import the data
Fmaxtheorique1 = readtable("D:\LGC_motiv\LGC_Motiv_results\pilots_v6_IP_Nback2_NOtaskSwitching_NOriskRepeatAfterFail\Fmax_theorique.xlsx", opts, "UseExcel", false);
Fmaxtheorique1 = table2array(Fmaxtheorique1);

% Clear temporary variables
clear opts

% compute fmax_theoritical
for j = 1:nb_pilots
    predictedForce(j) = Emax_morpho(Fmaxtheorique1(j,1),Fmaxtheorique1(j,2),Fmaxtheorique1(j,4),Fmaxtheorique1(j,3));
end

% transform MVC into Newtons (you have Volts. Divide by nominal output (782 microV/kgf *gain of the
% machine 200) then transform kgf into netwton by dividing by 0.1019716
init_MVC = ((init_MVC/(0.1564))/0.1019716 )* 3.128

%plot and compute correlation between initial/final MVC and theoritical strength
init_MVC = ((init_MVC/(0.1564))/0.1019716 )* 3.128

% init_MVC(17:19) = [381.7337 , 294.1015 , 766.8410];
% predictedForce(17:19) = [1.9463,1.4995,3.9098 ];
figure()
scatter(init_MVC,predictedForce)
corrcoef(init_MVC,predictedForce)

figure()
scatter(end_MVC,predictedForce)
corrcoef(end_MVC,predictedForce)
% % go back to folder with scripts
cd(main_task_folder_analysis);

