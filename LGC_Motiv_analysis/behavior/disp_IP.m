function[mean_deltaIP_Ep, mean_deltaIP_Em, sd_deltaIP_Ep, sd_deltaIP_Em, pval] = disp_IP(computerRoot, study_nm, showIndiv)
% [mean_IP_Ep, mean_IP_Em, sd_IP_Ep, sd_IP_Em, pval] = disp_IP(study_nm, showIndiv)
% disp_IP will display indifference point for physical and mental effort
% task.
%
% INPUTS
% computerRoot: pathway where data is
%
% study_nm: study name
%
% showIndiv:
% (0) only show mean data on the graph
% (1) link the mental and physical IP for each individual (by default)
%
% OUTPUT
% mean_IP_Ep: mean indifference point physical effort
% 
% mean_IP_Em: mean indifference point mental effort
%
% sd_IP_Ep: standard deviation indifference point physical effort
% 
% sd_IP_Em: standard deviation indifference point mental effort
%
% pval: p.value for t.test comparing physical and mental indifference point

%% if root not defined => ask for it
if ~exist('computerRoot','var') || isempty(computerRoot)
    computerRoot = LGCM_root_paths;
end

%% study names
if ~exist('study_nm','var') || isempty(study_nm)
    study_nm = 'study1';
end

%% link between indifference points
if ~exist('showIndiv','var') || isempty(showIndiv) || ~ismember(showIndiv,[0,1])
    showIndiv = 1;
end

%% working directories
dataRoot = [computerRoot, filesep, study_nm, filesep];

%% extract subjects
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

[mean_deltaIP_perSub_Ep, mean_deltaIP_perSub_Em] = deal(NaN(1,NS));

for iS = 1:NS
    sub_nm = ['CID',subject_id{iS}];
    subFolder = [dataRoot, sub_nm, filesep, 'behavior', filesep];
    IPdata = getfield(load([subFolder,'delta_IP_',sub_nm,'.mat'],'IP_variables'),'IP_variables');
    mean_deltaIP_perSub_Ep(iS) = IPdata.physicalDeltaIP;
    mean_deltaIP_perSub_Em(iS) = IPdata.mentalDeltaIP;
end % subject loop

%% average and SD
mean_deltaIP_Ep = mean(mean_deltaIP_perSub_Ep,2);
mean_deltaIP_Em = mean(mean_deltaIP_perSub_Em,2);
sd_deltaIP_Ep = mean(mean_deltaIP_perSub_Ep,2);
sd_deltaIP_Em = mean(mean_deltaIP_perSub_Em,2);

%% compare both IP
[~,pval] = ttest(mean_deltaIP_perSub_Ep, mean_deltaIP_perSub_Em);

%% figure
pSize = 30;
lSize = 2;

fig;

% show average data across subjects
hdl = bar(1:2, [mean_deltaIP_Ep, mean_deltaIP_Em],'LineWidth',lSize);
hdl.FaceColor = [143 143 143]./255;
hold on;
errorbar(1, mean_deltaIP_Ep, sd_deltaIP_Ep, 'k','LineWidth',lSize);
errorbar(2, mean_deltaIP_Em, sd_deltaIP_Em, 'k','LineWidth',lSize);

% show individual data
if showIndiv == 1
    scatter(1, mean_deltaIP_perSub_Ep,...
        'LineWidth',lSize,'MarkerFaceColor','none',...
        'MarkerEdgeColor','k','Marker','o');
    scatter(2, mean_deltaIP_perSub_Em,...
        'LineWidth',lSize,'MarkerFaceColor','none',...
        'MarkerEdgeColor','k','Marker','o');
    for iS = 1:NS
        plot(1:2, [mean_deltaIP_perSub_Ep(iS),...
            mean_deltaIP_perSub_Em(iS)],...
            'k','LineWidth',lSize);
    end
end

% legend
xticks(1:2);
xticklabels({'Ep','Em'});
ylabel('delta indifference point');
legend_size(pSize);

end % function