function[] = plot_model_comparison2(out_BMC,AIC,MAE,BIC,R2, RMSE, model_names)
% [] = plot_model_comparison(out_BMC,AIC,MAE,BIC,R2, RMSE, model_names)
% plot_model_comparison will plot one main figure summarizing the different
% model quality metrics across the different models and metrics entered.
% Similar to plot_model_comparison.m but also includes RMSE
%
% INPUTS
% out_BMC: output from VBA_groupBMC.m containing exceedance probability and
% estimated model frequency
%
% AIC: Akaike's Information Criterion for each model
%
% MAE: mean absolute error (mean absolute error between model predictions and
% actual data mean(abs(observed-predicted)))
%
% BIC: Bayesian Information Criterion for each model
%
% R²: correlation coefficient for each model
%
% RMSE: average root-mean squared error: square root (sqrt) of the mean squared error between actual data and model
% predictions (sqrt(mean(observed-predicted)^2)))
%
% model_names: variable with model names
%
% Adapted by N.Clairis from A.Barakat

%% define names manually for figures
% model_names = {'kI','kI+kE','kI+kE+kB','kI+kE+kB+kTime'};

%% extract SEM or SD for each metric
% sem_R2 = nanstd(R2,0,2)/sqrt(size(R2,2));
% sem_AIC = nanstd(AIC,0,2)/sqrt(size(AIC,2));
% sem_BIC = nanstd(BIC,0,2)/sqrt(size(BIC,2));
% sem_MAE = nanstd(MAE,0,2)/sqrt(size(MAE,2));
sd_R2 = std(R2,0,2,'omitnan');
sd_AIC = std(AIC,0,2,'omitnan');
sd_BIC = std(BIC,0,2,'omitnan');
sd_MAE = std(MAE,0,2,'omitnan');
sd_RMSE = std(RMSE,0,2,'omitnan');

%% open figure
fig;

%% exceedance probability subplot
subplot(3,3,1)
bar(out_BMC.ep,0.7,'FaceColor',[0.75 .75 .75],'EdgeColor',[0 0 0],'LineWidth',1.5)
ylim([0 1.1])
ylabel(compose('Exceedance\n probability'))
ax = gca;
ax.XTick      = linspace(1,length(model_names),length(model_names));
ax.XTickLabel = model_names;
ax.XTickLabelRotation = 30;
set(gca,'FontSize',15);

%% estimated model frequency subplot
subplot(3,3,2)
bar(out_BMC.Ef,0.7,'FaceColor',[0.75 .75 .75],'EdgeColor',[0 0 0],'LineWidth',1.5)
ylabel(compose('Estimated model\n frequency'))
ax = gca;
ax.XTick      = linspace(1,length(model_names),length(model_names));
ax.XTickLabel = model_names;
ax.XTickLabelRotation = 30;
set(gca,'FontSize',15);

%% R² subplot
subplot(3,3,3);
b = bar(mean(R2,2,'omitnan'),0.7,'FaceColor',[0.75 .75 .75],'EdgeColor',[0 0 0],'LineWidth',1.5);
ylim([0 0.75])
[ngroups,nbars] = size(mean(R2,2));
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
hold on
errorbar(x',mean(R2,2,'omitnan'),sd_R2,'k','linestyle','none','HandleVisibility','off','LineWidth',1.5);
hold off
ylabel('R2')
ax = gca;
ax.XTick      = linspace(1,length(model_names),length(model_names));
ax.XTickLabel = model_names;
ax.XTickLabelRotation = 30;
set(gca,'FontSize',15);

%% AIC subplot
subplot(3,3,4)
ax = gca;
set(gca,'FontSize',18)
b = bar(150-mean(AIC,2,'omitnan'),0.7,'FaceColor',[0.75 .75 .75],'EdgeColor',[0 0 0],'LineWidth',1.5);
[ngroups,nbars] = size(mean(AIC,2));
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
hold on
errorbar(x',150-mean(AIC,2,'omitnan'),sd_AIC,'k','linestyle','none','HandleVisibility','off','LineWidth',1.5);
hold off
ylabel('AIC')
ax = gca;
ax.XTick      = linspace(1,length(model_names),length(model_names));
ax.XTickLabel = model_names;
ax.XTickLabelRotation = 30;
set(gca,'FontSize',15);

%% BIC subplot
subplot(3,3,5);
b = bar(150-mean(BIC,2,'omitnan'),0.7,'FaceColor',[0.75 .75 .75],'EdgeColor',[0 0 0],'LineWidth',1.5);
[ngroups,nbars] = size(mean(BIC,2));
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
hold on
errorbar(x',150-mean(BIC,2,'omitnan'),sd_BIC,'k','linestyle','none','HandleVisibility','off','LineWidth',1.5);
hold off

ylabel('BIC')
ax = gca;
ax.XTick      = linspace(1,length(model_names),length(model_names));
ax.XTickLabel = model_names;
ax.XTickLabelRotation = 30;
set(gca,'FontSize',15);

%% MAE subplot
subplot(3,3,6);
b = bar(mean(MAE,2,'omitnan'),0.7,'FaceColor',[0.75 .75 .75],'EdgeColor',[0 0 0],'LineWidth',1.5);
ylim([0 0.5])
[ngroups,nbars] = size(mean(MAE,2));
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
hold on
errorbar(x',mean(MAE,2,'omitnan'),sd_MAE,'k','linestyle','none','HandleVisibility','off','LineWidth',1.5);
hold off

ylabel('MAE')
ax = gca;
ax.XTick      = linspace(1,length(model_names),length(model_names));
ax.XTickLabel = model_names;
ax.XTickLabelRotation = 30;
set(gca,'FontSize',15);

%% RMSE subplot
subplot(3,3,7);
b = bar(mean(RMSE,2,'omitnan'),0.7,'FaceColor',[0.75 .75 .75],'EdgeColor',[0 0 0],'LineWidth',1.5);
ylim([0 0.5]);
[ngroups,nbars] = size(mean(RMSE,2));
% Get the x coordinate of the bars
x = nan(nbars, ngroups);
for i = 1:nbars
    x(i,:) = b(i).XEndPoints;
end
% Plot the errorbars
hold on
errorbar(x',mean(RMSE,2,'omitnan'),sd_RMSE,'k','linestyle','none','HandleVisibility','off','LineWidth',1.5);
hold off

ylabel('RMSE')
ax = gca;
ax.XTick      = linspace(1,length(model_names),length(model_names));
ax.XTickLabel = model_names;
ax.XTickLabelRotation = 30;
set(gca,'FontSize',15)


end % function