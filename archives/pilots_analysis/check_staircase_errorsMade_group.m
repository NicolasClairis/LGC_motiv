
% root = fullfile('C:','Users','Loco','Documents','GitHub','LGC_motiv','LGC_Motiv_results');
root = fullfile('C:','Users','clairis','Desktop','GitHub','LGC_motiv','LGC_Motiv_results');

groups = {'pilots_v0',...
    'pilots_v1_IP_Nback1',...
    'pilots_v2_IP_Nback2',...
    'pilots_v3_IP_Nback2',...
    'pilots_v4_IP_Nback2_NOtaskSwitching',...
    'pilots_v5_IP_Nback2_NOtaskSwitching'};
groupSelection = listdlg('ListString',groups);
groupFolderName = groups{groupSelection};

switch groupFolderName
    case 'pilots_v0'
        subInit = {'MC','AC','HC','AR'};
    case 'pilots_v1_IP_Nback1'
        subInit = {'RC','DK','MC','AM'};
    case 'pilots_v2_IP_Nback2'
        subInit = {'SR'};
    case 'pilots_v3_IP_Nback2'
        subInit = {'DU'};
    case 'pilots_v4_IP_Nback2_NOtaskSwitching'
        subInit = {'SN','AD','JA'};
    case 'pilots_v5_IP_Nback2_NOtaskSwitching'
        subInit = {''};
end
NS = length(subInit);
groupFolder = [root, filesep, groupFolderName, filesep];

% initialize variable
[n_AVGerrorsEffort.Em.E1, n_AVGerrorsEffort.Em.E2, n_AVGerrorsEffort.Em.E3,...
    n_success.Em.E1, n_success.Em.E2, n_success.Em.E3,...
    nChoice.Em.E1, nChoice.Em.E2, nChoice.Em.E3] = deal(NaN(1,NS));

for iSub = 1:NS
    subInitials = subInit{iSub};
    % load data
    loadStruct = load([groupFolder,'IP_pilot_data',subInitials,'_sub_',num2str(iSub),'.mat']);
    
    [n_AVGerrorsEffort_tmp, n_success_tmp, nChoice_tmp] = check_staircase_errorsMade_perSub(groupFolder, subInitials, iSub);
    n_AVGerrorsEffort.Em.E1(iSub) = n_AVGerrorsEffort_tmp.E1;
    n_AVGerrorsEffort.Em.E2(iSub) = n_AVGerrorsEffort_tmp.E2;
    n_AVGerrorsEffort.Em.E3(iSub) = n_AVGerrorsEffort_tmp.E3;
    n_success.Em.E1(iSub) = n_success_tmp.E1;
    n_success.Em.E2(iSub) = n_success_tmp.E2;
    n_success.Em.E3(iSub) = n_success_tmp.E3;
    nChoice.Em.E1(iSub) = nChoice_tmp.Em.E1;
    nChoice.Em.E2(iSub) = nChoice_tmp.Em.E2;
    nChoice.Em.E3(iSub) = nChoice_tmp.Em.E3;
end % subject loop

%% average across subjects
[m_n_AVGerrorsEffort.Em.E1, sem_n_AVGerrorsEffort.Em.E1] = mean_sem_sd(n_AVGerrorsEffort.Em.E1, 2);
[m_n_AVGerrorsEffort.Em.E2, sem_n_AVGerrorsEffort.Em.E2] = mean_sem_sd(n_AVGerrorsEffort.Em.E2, 2);
[m_n_AVGerrorsEffort.Em.E3, sem_n_AVGerrorsEffort.Em.E3] = mean_sem_sd(n_AVGerrorsEffort.Em.E3, 2);
[m_n_success.Em.E1, sem_n_success.Em.E1] = mean_sem_sd(n_success.Em.E1, 2);
[m_n_success.Em.E2, sem_n_success.Em.E2] = mean_sem_sd(n_success.Em.E2, 2);
[m_n_success.Em.E3, sem_n_success.Em.E3] = mean_sem_sd(n_success.Em.E3, 2);
[m_nChoice.Em.E1, sem_nChoice.Em.E1] = mean_sem_sd(nChoice.Em.E1, 2);
[m_nChoice.Em.E2, sem_nChoice.Em.E2] = mean_sem_sd(nChoice.Em.E2, 2);
[m_nChoice.Em.E3, sem_nChoice.Em.E3] = mean_sem_sd(nChoice.Em.E3, 2);

%% figure

% errors
fig;
bar(1:3,...
    [m_n_AVGerrorsEffort.Em.E1, m_n_AVGerrorsEffort.Em.E2, m_n_AVGerrorsEffort.Em.E3]);
hold on;
errorbar(1:3,...
    [m_n_AVGerrorsEffort.Em.E1, m_n_AVGerrorsEffort.Em.E2, m_n_AVGerrorsEffort.Em.E3],...
    [sem_n_AVGerrorsEffort.Em.E1, sem_n_AVGerrorsEffort.Em.E2, sem_n_AVGerrorsEffort.Em.E3],...
    ' ','LineWidth',3);
ylabel('errors');
xticks(1:3);
xticklabels({'E1','E2','E3'});
legend_size(50);

% success
fig;
bar(1:3,...
    [m_n_success.Em.E1, m_n_success.Em.E2, m_n_success.Em.E3]);
hold on;
errorbar(1:3,...
    [m_n_success.Em.E1, m_n_success.Em.E2, m_n_success.Em.E3],...
    [sem_n_success.Em.E1, sem_n_success.Em.E2, sem_n_success.Em.E3],...
    ' ','LineWidth',3);
ylabel('success');
xticks(1:3);
xticklabels({'E1','E2','E3'});
legend_size(50);

% choices
fig;
bar(1:3,...
    [m_nChoice.Em.E1, m_nChoice.Em.E2, m_nChoice.Em.E3]);
hold on;
errorbar(1:3,...
    [m_nChoice.Em.E1, m_nChoice.Em.E2, m_nChoice.Em.E3],...
    [sem_nChoice.Em.E1, sem_nChoice.Em.E2, sem_nChoice.Em.E3],...
    ' ','LineWidth',3);
ylabel('choice');
xticks(1:3);
xticklabels({'E1','E2','E3'});
legend_size(50);