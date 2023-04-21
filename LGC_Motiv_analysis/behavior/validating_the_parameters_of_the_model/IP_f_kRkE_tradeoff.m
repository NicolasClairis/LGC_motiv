% check at correlation between indifference point in each task and
% parameters extracted with the model. In particular, it should correlate
% with the reward/effort tradeoff.

%% subject selection
study_nm = 'study1';
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% working directories
computerRoot = LGCM_root_paths;
dataRoot = [computerRoot, filesep, study_nm, filesep];

%% load parameters
[mdlType, mdlN] = behavioral_model_selection;
[prm, mdlType, mdlN] = prm_extraction(study_nm, subject_id, mdlType, mdlN);
% extract parameters of interest
kR = prm.kR;
kEp = prm.kEp;
kEm = prm.kEm;

%% load indifference point
[IP.Ep, IP.Em] = deal(NaN(1,NS));
for iS = 1:NS
    sub_nm = ['CID',subject_id{iS}];
    subFolder = [dataRoot, sub_nm, filesep, 'behavior', filesep];
    IPdata = getfield(load([subFolder,'delta_IP_',sub_nm,'.mat'],...
        'IP_variables'),'IP_variables');
    IP.Ep(iS) = IPdata.baselineR + IPdata.physicalDeltaIP;
    IP.Em(iS) = IPdata.baselineR + IPdata.mentalDeltaIP;
end % subject loop

%% perform the correlation
% physical effort
kRkEp_tradeoff = nanzscore(kR)./nanzscore(kEp);
goodSubsEp = ~isnan(kRkEp_tradeoff);
% goodSubsEp = ~isnan(kRkEp_tradeoff) &...
%     (kRkEp_tradeoff <= mean(kRkEp_tradeoff,2,'omitnan') + 3*std(kRkEp_tradeoff,0,2,'omitnan')) &...
%     (kRkEp_tradeoff >= mean(kRkEp_tradeoff,2,'omitnan') - 3*std(kRkEp_tradeoff,0,2,'omitnan'));
[beta.Ep.kRkEtradeoff, ~,stats.Ep.kRkEtradeoff] = glmfit(kRkEp_tradeoff(goodSubsEp), IP.Ep(goodSubsEp),'normal');
kRkEp_sorted = sort(kRkEp_tradeoff(goodSubsEp));
IP_Ep_kRkE_fit = glmval(beta.Ep.kRkEtradeoff, kRkEp_sorted, 'identity');

kR_min_kEp = nanzscore(kR)-nanzscore(kEp);
goodSubsEp_bis = ~isnan(kR_min_kEp);
% goodSubsEp_bis = ~isnan(kR_min_kEp) &...
%     (kR_min_kEp <= mean(kR_min_kEp,2,'omitnan') + 3*std(kR_min_kEp,0,2,'omitnan')) &...
%     (kR_min_kEp >= mean(kR_min_kEp,2,'omitnan') - 3*std(kR_min_kEp,0,2,'omitnan'));
[beta.Ep.kR_min_kE, ~,stats.Ep.kR_min_kE] = glmfit(kR_min_kEp(goodSubsEp_bis), IP.Ep(goodSubsEp_bis),'normal');
kR_min_kEp_sorted = sort(kR_min_kEp(goodSubsEp_bis));
IP_Ep_kRminkE_fit = glmval(beta.Ep.kR_min_kE, kR_min_kEp_sorted, 'identity');

% mental effort
kRkEm_tradeoff = nanzscore(kR)./nanzscore(kEm);
goodSubsEm = ~isnan(kRkEm_tradeoff);
% goodSubsEm = ~isnan(kRkEm_tradeoff) &...
%     (kRkEm_tradeoff <= mean(kRkEm_tradeoff,2,'omitnan') + 3*std(kRkEm_tradeoff,0,2,'omitnan')) &...
%     (kRkEm_tradeoff >= mean(kRkEm_tradeoff,2,'omitnan') - 3*std(kRkEm_tradeoff,0,2,'omitnan'));
[beta.Em.kRkEtradeoff, ~,stats.Em.kRkEtradeoff] = glmfit(kRkEm_tradeoff(goodSubsEm), IP.Em(goodSubsEm),'normal');
kRkEm_sorted = sort(kRkEm_tradeoff(goodSubsEm));
IP_Em_fit = glmval(beta.Em.kRkEtradeoff, kRkEm_sorted, 'identity');

kR_min_kEm = nanzscore(kR)-nanzscore(kEm);
goodSubsEm_bis = ~isnan(kR_min_kEm);
% goodSubsEm_bis = ~isnan(kR_min_kEm) &...
%     (kR_min_kEm <= mean(kR_min_kEm,2,'omitnan') + 3*std(kR_min_kEm,0,2,'omitnan')) &...
%     (kR_min_kEm >= mean(kR_min_kEm,2,'omitnan') - 3*std(kR_min_kEm,0,2,'omitnan'));
[beta.Em.kR_min_kE, ~,stats.Em.kR_min_kE] = glmfit(kR_min_kEm(goodSubsEm_bis), IP.Em(goodSubsEm_bis),'normal');
kR_min_kEm_sorted = sort(kR_min_kEm(goodSubsEm_bis));
IP_Em_kRminkE_fit = glmval(beta.Em.kR_min_kE, kR_min_kEm_sorted, 'identity');

%% figure
% general figure infos
pSize = 30;
lSize = 2;
lWidth = 3;
grey = [143 143 143]./255;

% kR/kE
fig;
% physical effort
subplot(1,2,1);
hold on;
data_hdl = scatter(kRkEp_tradeoff(goodSubsEp),...
    IP.Ep(goodSubsEp));
data_hdl.LineWidth = lWidth;
data_hdl.MarkerEdgeColor = 'k';
fit_hdl = plot(kRkEp_sorted, IP_Ep_kRkE_fit);
fit_hdl.LineWidth = lWidth;
fit_hdl.Color = grey;
fit_hdl.LineStyle = '--';
xlabel('kR/kEp');
ylabel('IP');
legend_size(pSize);

% mental effort
subplot(1,2,2);
hold on;
data_hdl = scatter(kRkEm_tradeoff(goodSubsEm),...
    IP.Em(goodSubsEm));
data_hdl.LineWidth = lWidth;
data_hdl.MarkerEdgeColor = 'k';
fit_hdl = plot(kRkEm_sorted, IP_Em_fit);
fit_hdl.LineWidth = lWidth;
fit_hdl.Color = grey;
fit_hdl.LineStyle = '--';
xlabel('kR/kEm');
ylabel('IP');
legend_size(pSize);

% kR - kE
fig;
% physical effort
subplot(1,2,1);
hold on;
data_hdl = scatter(kR_min_kEp(goodSubsEp_bis),...
    IP.Ep(goodSubsEp_bis));
data_hdl.LineWidth = lWidth;
data_hdl.MarkerEdgeColor = 'k';
fit_hdl = plot(kR_min_kEp_sorted, IP_Ep_kRminkE_fit);
fit_hdl.LineWidth = lWidth;
fit_hdl.Color = grey;
fit_hdl.LineStyle = '--';
xlabel('kR-kEp');
ylabel('IP');
legend_size(pSize);

% mental effort
subplot(1,2,2);
hold on;
data_hdl = scatter(kR_min_kEm(goodSubsEm_bis),...
    IP.Em(goodSubsEm_bis));
data_hdl.LineWidth = lWidth;
data_hdl.MarkerEdgeColor = 'k';
fit_hdl = plot(kR_min_kEm_sorted, IP_Em_kRminkE_fit);
fit_hdl.LineWidth = lWidth;
fit_hdl.Color = grey;
fit_hdl.LineStyle = '--';
xlabel('kR-kEm');
ylabel('IP');
legend_size(pSize);