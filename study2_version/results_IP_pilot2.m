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

nb_sessions=4;

% Encode physical or manual

mentalFirst = [0 0 1 0 1 0 1 0 1];
cd ..
% manually put the number of pilots
for i_pilot = 1:length(mentalFirst)
    
subjectID = ['CID',num2str(i_pilot+200)];
% go to subject path

rootPath = [pwd, filesep];
subPath = [rootPath,filesep,subjectID,filesep,'behavior',filesep];

    mentalCalib_nm = [subjectID,'_mentalCalib.mat'];
    physicalCalib_nm = [subjectID,'_physicalCalib.mat'];
    
    training_nm = ['training_data_',subjectID,'.mat'];
%     start_mental = getfield(getfield(load([subPath,training_nm],'all'),'all'),'mentalfirst');
    
    for i_session=1:nb_sessions
        if mod(mentalFirst(i_pilot)+i_session -1,2) == 0
       test_nm.(['session_',num2str(i_session)]) = [subjectID,'_session',num2str(i_session),'_physical_task_messyAllStuff.mat'];
        else
       test_nm.(['session_',num2str(i_session)]) = [subjectID,'_session',num2str(i_session),'_mental_task_messyAllStuff.mat'];
        end
    end
    
    %% loading part
    mentalCalib(i_pilot) = load([subPath,mentalCalib_nm],'NMP');
    physicalCalib(i_pilot) = load([subPath,physicalCalib_nm],'MVC');
    
    pilot_nm = ['pilot',num2str(i_pilot)];
    all.(pilot_nm).ses1 = load([subPath,test_nm.session_1]);
    all.(pilot_nm).ses2 = load([subPath,test_nm.session_2]);
    all.(pilot_nm).ses3 = load([subPath,test_nm.session_3]);
    all.(pilot_nm).ses4 = load([subPath,test_nm.session_4]);
    
    if mentalFirst == 1
    percentagePerf_phys1(i_pilot,:) = all.(pilot_nm).ses2.perfSummary.percentagePerf;
    percentagePerf_phys2(i_pilot,:) = all.(pilot_nm).ses4.perfSummary.percentagePerf;
    percentagePerf_ment1(i_pilot,:) = all.(pilot_nm).ses1.perfSummary.percentagePerf;
    percentagePerf_ment2(i_pilot,:) = all.(pilot_nm).ses3.perfSummary.percentagePerf;
    else
    percentagePerf_phys1(i_pilot,:) = all.(pilot_nm).ses1.perfSummary.percentagePerf;
    percentagePerf_phys2(i_pilot,:) = all.(pilot_nm).ses3.perfSummary.percentagePerf;
    percentagePerf_ment1(i_pilot,:) = all.(pilot_nm).ses2.perfSummary.percentagePerf;
    percentagePerf_ment2(i_pilot,:) = all.(pilot_nm).ses4.perfSummary.percentagePerf;
    end     
end
%% AVG HERE
for i = 1:6
avgPhysical1(i) = mean(percentagePerf_phys1((i-1)*9+1:i*9));
avgPhysical2(i) = mean(percentagePerf_phys2((i-1)*9+1:i*9));
avgMental1(i) = mean(percentagePerf_ment1((i-1)*9+1:i*9));
avgMental2(i) = mean(percentagePerf_ment2((i-1)*9+1:i*9));
end

%% CORRELATION HERE

%% PLOT HERE

figure()

plot(avgPhysical1);
hold on
plot(avgPhysical2);
plot(avgMental1);
plot(avgMental2);
legend('phys1','phys2','ment1','ment2')
xlabel('bins of avg of 9 trials')
ylabel('performance')

%     if i_pilot > 3 && i_pilot ~= 10
%         for q = 1:90
%             totalTime(i_pilot,q) = all.mental.learning.extendedLearning.(strcat('trial_',num2str(q))).totalTime_success;
%             errorsMade(i_pilot,q) = all.mental.learning.extendedLearning.(strcat('trial_',num2str(q))).n_errorsMade;
%         end
%             for m = 1:18
%         meanTotalTime(m) = mean(totalTime(i_pilot,1 + (m-1)*5:m*5));
%     end
%     end
%     
% 
%     delta_IP(i_pilot,1) = -1.5 + all.physical.EffortLvl_1.session_nb1.repeat_nb1.perfSummary.IP;
%     delta_IP(i_pilot,2) = -1.5 + all.physical.EffortLvl_1.session_nb1.repeat_nb2.perfSummary.IP;
%     delta_IP(i_pilot,3) = -1.5 + all.mental.EffortLvl_1.session_nb1.repeat_nb1.perfSummary.IP;
%     delta_IP(i_pilot,4) = -1.5 + all.mental.EffortLvl_1.session_nb1.repeat_nb2.perfSummary.IP;
%     delta_IP(i_pilot,5) = 1.5 - all.physical.EffortLvl_1.session_nb2.repeat_nb1.perfSummary.IP;
%     delta_IP(i_pilot,6) = 1.5 - all.physical.EffortLvl_1.session_nb2.repeat_nb2.perfSummary.IP;
%     delta_IP(i_pilot,7) = 1.5 - all.mental.EffortLvl_1.session_nb2.repeat_nb1.perfSummary.IP;
%     delta_IP(i_pilot,8) = 1.5 - all.mental.EffortLvl_1.session_nb2.repeat_nb2.perfSummary.IP;
%     delta_IP(i_pilot,9) = -1.5 + all.physical.EffortLvl_2.session_nb1.repeat_nb1.perfSummary.IP;
%     delta_IP(i_pilot,10) = -1.5 + all.physical.EffortLvl_2.session_nb1.repeat_nb2.perfSummary.IP;
%     delta_IP(i_pilot,11) = -1.5 + all.mental.EffortLvl_2.session_nb1.repeat_nb1.perfSummary.IP;
%     delta_IP(i_pilot,12) = -1.5 + all.mental.EffortLvl_2.session_nb1.repeat_nb2.perfSummary.IP;
%     delta_IP(i_pilot,13) = 1.5 - all.physical.EffortLvl_2.session_nb2.repeat_nb1.perfSummary.IP;
%     delta_IP(i_pilot,14) = 1.5 - all.physical.EffortLvl_2.session_nb2.repeat_nb2.perfSummary.IP;
%     delta_IP(i_pilot,15) = 1.5 - all.mental.EffortLvl_2.session_nb2.repeat_nb1.perfSummary.IP;
%     delta_IP(i_pilot,16) = 1.5 - all.mental.EffortLvl_2.session_nb2.repeat_nb2.perfSummary.IP;
%     if exist('initial_MVC') && exist('last_MVC') 
%         delta_MVC(i_pilot) = (last_MVC.MVC - initial_MVC.MVC)/ initial_MVC.MVC * 100;
%         init_MVC(i_pilot) = initial_MVC.MVC;
%         end_MVC(i_pilot) = last_MVC.MVC;
%     end
%     
%     if exist('t_min_calib') && exist('t_min_lastCalib')
%         delta_MVM(i_pilot) = t_min_lastCalib - t_min_calib;
%         init_MVM(i_pilot) = t_min_calib;
%         end_MVM(i_pilot) = t_min_lastCalib;
%     end
%     
%     for i_mean = 1:length(delta_IP)/2
%         mean_IP(i_pilot, i_mean) = mean(delta_IP(i_pilot,(i_mean-1)*2 + 1:i_mean*2),2);
%         std_IP(i_pilot, i_mean) = std(delta_IP(i_pilot,(i_mean-1)*2 + 1:i_mean*2));
%         %  1       2       3       4       5      6       7        8
%         % E1PhRe, E1MeRe, E1PhPu, E1MePu, E2PhRe, E2MeRe, E2PhPu, E2MePu
%     end
%         % prepare data for the difference in reward and punishment for low effort
%     t_test_re_pu_E1(i_pilot,1) = mean(mean_IP(i_pilot,1:2),2);
%     t_test_re_pu_E1(i_pilot,2) = mean(mean_IP(i_pilot,3:4),2);
%     
%     % prepare data for the difference in physical and mental low efforts
%     t_test_ph_me_E1(i_pilot,1) = mean(mean_IP(i_pilot,[1,3]),2);
%     t_test_ph_me_E1(i_pilot,2) = mean(mean_IP(i_pilot,[2,4]),2);
%     
%     % prepare data for the difference in reward and punishment for high effort
%     t_test_re_pu_E2(i_pilot,1) = mean(mean_IP(i_pilot,5:6),2);
%     t_test_re_pu_E2(i_pilot,2) = mean(mean_IP(i_pilot,7:8),2);
%     
%     % prepare data for the difference in physical and mental high effort
%     t_test_ph_me_E2(i_pilot,1) = mean(mean_IP(i_pilot,[5,7]),2);
%     t_test_ph_me_E2(i_pilot,2) = mean(mean_IP(i_pilot,[6,8]),2);
%     
%     % prepare data for the difference between low and high efforts, indep of other conditions
%     t_test_E2(i_pilot,1) = mean(mean_IP(i_pilot,1:4),2);
%     t_test_E2(i_pilot,2) = mean(mean_IP(i_pilot,5:8),2);
%     
%     % prepare data for the difference in reward and punishment, indep of other conditions
%     t_test_re_pu(i_pilot,1) = mean(mean_IP(i_pilot,[1,2,5,6]),2);
%     t_test_re_pu(i_pilot,2) = mean(mean_IP(i_pilot,[3,4,7,8]),2);
%     
%     % prepare data for the difference in physical and mental effort, indep of other conditions
%     t_test_ph_me(i_pilot,1) = mean(mean_IP(i_pilot,[1,3,5,7]),2);
%     t_test_ph_me(i_pilot,2) = mean(mean_IP(i_pilot,[2,4,6,8]),2);
% %     % prepare data for the difference in reward and punishment for low effort
% %     t_test_re_pu_E1(i_pilot,1) = mean(delta_IP(i_pilot,1:2),2);
% %     t_test_re_pu_E1(i_pilot,2) = mean(delta_IP(i_pilot,3:4),2);
% %     
% %     % prepare data for the difference in physical and mental low efforts
% %     t_test_ph_me_E1(i_pilot,1) = mean(delta_IP(i_pilot,[1,3]),2);
% %     t_test_ph_me_E1(i_pilot,2) = mean(delta_IP(i_pilot,[2,4]),2);
% %     
% %     % prepare data for the difference in reward and punishment for high effort
% %     t_test_re_pu_E2(i_pilot,1) = mean(delta_IP(i_pilot,5:6),2);
% %     t_test_re_pu_E2(i_pilot,2) = mean(delta_IP(i_pilot,7:8),2);
% %     
% %     % prepare data for the difference in physical and mental high effort
% %     t_test_ph_me_E2(i_pilot,1) = mean(delta_IP(i_pilot,[5,7]),2);
% %     t_test_ph_me_E2(i_pilot,2) = mean(delta_IP(i_pilot,[6,8]),2);
% %     
% %     % prepare data for the difference between low and high efforts, indep of other conditions
% %     t_test_E2(i_pilot,1) = mean(delta_IP(i_pilot,1:4),2);
% %     t_test_E2(i_pilot,2) = mean(delta_IP(i_pilot,5:8),2);
% %     
% %     % prepare data for the difference in reward and punishment, indep of other conditions
% %     t_test_re_pu(i_pilot,1) = mean(delta_IP(i_pilot,[1,2,5,6]),2);
% %     t_test_re_pu(i_pilot,2) = mean(delta_IP(i_pilot,[3,4,7,8]),2);
% %     
% %     % prepare data for the difference in physical and mental effort, indep of other conditions
% %     t_test_ph_me(i_pilot,1) = mean(delta_IP(i_pilot,[1,3,5,7]),2);
% %     t_test_ph_me(i_pilot,2) = mean(delta_IP(i_pilot,[2,4,6,8]),2);
% 
% %     end
%     
% 
%     %% remove outlier
%     i_outlier = 9;
%     remove_outlier = false;
%     if remove_outlier == true
%     delta_IP(i_outlier,:) = [];
%     delta_MVC(i_outlier) = [];
%     delta_MVM(i_outlier) = [];
%     mean_IP(i_outlier,:) = [];
%     std_IP(i_outlier,:) = [];
%     init_MVC(i_outlier) = [];
%     end_MVC(i_outlier) = [];
%     init_MVM(i_outlier) = [];
%     end_MVM(i_outlier) = [];
%     
% 
%     t_test_re_pu_E1(i_outlier,:) = [];
%     t_test_ph_me_E1(i_outlier,:) = [];
%     t_test_re_pu_E2(i_outlier,:) = [];
%     t_test_ph_me_E2(i_outlier,:) = [];
%     t_test_E2(i_outlier,:) = [];
%     t_test_re_pu(i_outlier,:) = [];
%     t_test_ph_me(i_outlier,:) = [];
%     end
% %% t-test on our pilots. between conditions and calibrations
% % test with a t test the difference in reward and punishment for low effort
% [h_re_pu_E1, p_re_pu_E1] = ttest(t_test_re_pu_E1(:,1),t_test_re_pu_E1(:,2))
% % test with a t test the difference in physical and mental low efforts
% [h_ph_me_E1, p_ph_me_E1] = ttest(t_test_re_pu_E1(:,1),t_test_re_pu_E1(:,2))
% 
% % test with a t test the difference in reward and punishment for high effort
% [h_re_pu_E2, p_re_pu_E2] = ttest(t_test_re_pu_E2(:,1),t_test_re_pu_E2(:,2))
% % test with a t test the difference in physical and mental high effort
% [h_ph_me_E2, p_ph_me_E2] = ttest(t_test_re_pu_E2(:,1),t_test_re_pu_E2(:,2))
% 
% % test with a t test the difference in reward and punishment, indep of other conditions
% [h_re_pu,p_re_pu] = ttest(t_test_re_pu(:,1),t_test_re_pu(:,2))
% % test with a t test the difference in physical and mental effort, indep of other conditions
% [h_ph_me,p_ph_me] = ttest(t_test_ph_me(:,1),t_test_ph_me(:,2))
% % test with a t test the difference between low and high efforts, indep of other conditions
% [h_E2, p_E2] = ttest(t_test_E2(:,1),t_test_E2(:,2))
% 
% % test with t test if MVC is reduced in participants
% [h_MVC, p_MVC] = ttest(init_MVC,end_MVC)
% % test with t test if MVM is better in participants
% [h_MVM, p_MVM] = ttest(init_MVM,end_MVM)
% 
% 
% %% plots
% figure()
% bar(delta_MVC)
% 
% figure()
% bar(delta_MVM)
% 
% figure()
% 
% b = bar(mean_IP', 'grouped');
% hold on
% set(gca,'XTick',[1 2 3 4 5 6 7 8])
% set(gca,'XTickLabel',{'E1PR','E1MR','E1PP','E1MP','E2PR','E2MR','E2PP','E2MP'})
% 
% % Calculate the number of groups and number of bars in each group
% [ngroups,nbars] = size(mean_IP);
% % Get the x coordinate of the bars
% x = nan(ngroups,nbars);
% for i = 1:ngroups
%     x(i,:) = b(i).XEndPoints;
% end
% errorbar(x, mean_IP,std_IP,'k','linestyle','none')
% er.Color = [0 0 0];
% er.LineStyle = 'none';
% hold off
% hold on
% set(gca,'XTick',[1 2 3 4 5 6 7 8])
% set(gca,'XTickLabel',{'E1-Ph-Re','E1-Me-Re','E1-Ph-Pu','E1-Me-Pu','E2-Ph-Re','E2-Me-Re','E2-Ph-Pu','E2-Me-Pu'})
% ylabel('Indifference Point (IP) CHF')
% title('E1/E2 low/high effort level, Ph/Me physical mental, Re/Pu reward punishment')
% 
% 
% %% Import fmax_theoritical
% opts = spreadsheetImportOptions("NumVariables", 4);
% 
% % Specify sheet and range
% opts.Sheet = "Sheet1";
% opts.DataRange = "C2:F17";
% 
% % Specify column names and types
% opts.VariableNames = ["Pliantrieur", "Plipostriur", "Longueurdelavantbras", "Circonfrence"];
% opts.VariableTypes = ["double", "double", "double", "double"];
% 
% % Import the data
% Fmaxtheorique1 = readtable("D:\LGC_motiv\LGC_Motiv_results\pilots_v6_IP_Nback2_NOtaskSwitching_NOriskRepeatAfterFail\Fmax_theorique.xlsx", opts, "UseExcel", false);
% Fmaxtheorique1 = table2array(Fmaxtheorique1);
% 
% % Clear temporary variables
% clear opts
% 
% % compute fmax_theoritical
% for j = 1:nb_pilots
%     predictedForce(j) = Emax_morpho(Fmaxtheorique1(j,1),Fmaxtheorique1(j,2),Fmaxtheorique1(j,4),Fmaxtheorique1(j,3));
% end
% %plot and compute correlation between initial/final MVC and theoritical strength
% init_MVC = ((init_MVC/(0.1564))/0.1019716 )* 3.128
% 
% % init_MVC(17:19) = [381.7337 , 294.1015 , 766.8410];
% % predictedForce(17:19) = [1.9463,1.4995,3.9098 ];
% figure()
% scatter(init_MVC,predictedForce)
% corrcoef(init_MVC,predictedForce)
% 
% figure()
% scatter(end_MVC,predictedForce)
% corrcoef(end_MVC,predictedForce)
% % % go back to folder with scripts
cd(main_task_folder_analysis);

