
subs = {'201','202','203','204','205'};
NS = length(subs);

[mean_IP_perSub_Ep, mean_IP_perSub_Em] = deal(NaN(1,NS));

for iS = 1:NS
    sub_nm = subs{iS};
    IPdata = getfield(load(fullfile(['CID',sub_nm],'behavior',['delta_IP_CID',sub_nm,'.mat']),'IP_variables'),'IP_variables');
    mean_IP_perSub_Ep(iS) = IPdata.physicalDeltaIP;
    mean_IP_perSub_Em(iS) = IPdata.mentalDeltaIP;
end % subject loop

%% average and SD
mean_IP_Ep = mean(mean_IP_perSub_Ep,2);
mean_IP_Em = mean(mean_IP_perSub_Em,2);
sd_IP_Ep = mean(mean_IP_perSub_Ep,2);
sd_IP_Em = mean(mean_IP_perSub_Em,2);

%% figure
pSize = 30;
lSize = 2;

fig;
bar(1:2, [mean_IP_Ep, mean_IP_Em]);
hold on;
errorbar(1, mean_IP_Ep, sd_IP_Ep, 'k','LineWidth',lSize);
errorbar(2, mean_IP_Em, sd_IP_Em, 'k','LineWidth',lSize);
scatter(1, mean_IP_perSub_Ep, 'LineWidth',lSize,'MarkerFaceColor','none','MarkerEdgeColor','k','Marker','o');
scatter(2, mean_IP_perSub_Em, 'LineWidth',lSize,'MarkerFaceColor','none','MarkerEdgeColor','k','Marker','o');
for iS = 1:NS
    plot(1:2, [mean_IP_perSub_Ep(iS), mean_IP_perSub_Em(iS)],'k','LineWidth',lSize);
end
xticks(1:2);
xticklabels({'Ep','Em'});
ylabel('delta indifference point');
legend_size(pSize);