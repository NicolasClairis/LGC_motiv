%% check_eNAMPT_correlates is a script to look at correlations between 
% several factors of interest (BMI, age, date of acquisition, gender, etc.)
% and the eNAMPT concentrations in the plasma in study 1.

%% load data
load('M:\Nicolas_Clairis\study1\subs_physio_variables.mat',...
    'BMI','age','gender','eNAMPT','datesExp','weight','height');

%% look at correlation between BMI and eNAMPT concentrations

figure;
curveHdl = scatter(BMI,eNAMPT);
curveHdl.MarkerEdgeColor = [0 0 0];
curveHdl.LineWidth = 2;
curveHdl.SizeData = 100;
hold on;

% add thresholds for anorexia/obesity
xThresholds = [18.5, 25, 30];
for iBMI = 1:length(xThresholds)
    jBMI = xThresholds(iBMI);
    hdl = line([jBMI jBMI],[0 30]);
    hdl.LineWidth = 3; hdl.Color = [0 0 0];
end

% add legends
xlabel('BMI');
ylabel('eNAMPT (ng/mL)');
legend_size(30);

%% look at correlation between weight and eNAMPT concentrations

figure;
curveHdl = scatter(weight,eNAMPT);
curveHdl.MarkerEdgeColor = [0 0 0];
curveHdl.LineWidth = 2;
curveHdl.SizeData = 100;

% add legends
xlabel('weight (kg)');
ylabel('eNAMPT (ng/mL)');
ylim([0 30])
legend_size(30);

% stats
goodSubs_weight = (~isnan(weight) & ~isnan(eNAMPT));
% do the stats
weight_corrcoeff = corr2(weight(goodSubs_weight), eNAMPT(goodSubs_weight));
[weight_betas, ~, weight_stats] = glmfit(weight(goodSubs_weight), eNAMPT(goodSubs_weight),'normal');

%% look at correlation between height and eNAMPT concentrations

figure;
curveHdl = scatter(height,eNAMPT);
curveHdl.MarkerEdgeColor = [0 0 0];
curveHdl.LineWidth = 2;
curveHdl.SizeData = 100;

% add legends
xlabel('height (m)');
ylabel('eNAMPT (ng/mL)');
ylim([0 30])
legend_size(30);

% stats
goodSubs_weight = (~isnan(weight) & ~isnan(eNAMPT));
% do the stats
weight_corrcoeff = corr2(weight(goodSubs_weight), eNAMPT(goodSubs_weight));
[weight_betas, ~, weight_stats] = glmfit(weight(goodSubs_weight), eNAMPT(goodSubs_weight),'normal');

%% look at correlation between gender and eNAMPT
figure;
curveHdl = scatter(gender,eNAMPT);
curveHdl.MarkerEdgeColor = [0 0 0];
curveHdl.LineWidth = 2;
curveHdl.SizeData = 100;

% add legends
xlabel('gender');
ylabel('eNAMPT (ng/mL)');
legend_size(30);

%% look at correlation between age and eNAMPT
figure;
curveHdl = scatter(age,eNAMPT);
curveHdl.MarkerEdgeColor = [0 0 0];
curveHdl.LineWidth = 2;
curveHdl.SizeData = 100;

% add legends
xlabel('age (y.o.)');
ylabel('eNAMPT (ng/mL)');
legend_size(30);
ylim([0 30]);

%% look correlation between date of experiment and eNAMPT

% convert dates of experiment into seconds
% datesExp = NaN(1,length(eNAMPT));
% NS = length(eNAMPT);
% for iS = 1:NS
%     datesExp(iS) = datenum(dates{iS},'DD.mm.yyyy');
% end

% display results
figure;
curveHdl = scatter(datesExp, eNAMPT);
curveHdl.MarkerEdgeColor = [0 0 0];
curveHdl.LineWidth = 2;
curveHdl.SizeData = 100;

% add legends
xlabel('date experiment (s)');
ylabel('eNAMPT (ng/mL)');
legend_size(30);
ylim([0 30]);

%% look correlation between CTQ and eNAMPT
ctq_scores = {'Emotional abuse','Physical abuse','Emotional neglect',...
    'Physical neglect','Sexual abuse','Minimization denial'};
ctq_fields = {'ea','pa','en','pn','sa','md'};
n_ctq_scores = length(ctq_scores);
figure;
for iCTQ = 1:n_ctq_scores
    ctq_field = ctq_fields{iCTQ};
    subplot(2,3,iCTQ);
    CTQ_var = CTQ.(ctq_field);
    curveHdl = scatter(CTQ_var, eNAMPT);
    curveHdl.MarkerEdgeColor = [0 0 0];
    curveHdl.LineWidth = 2;
    curveHdl.SizeData = 100;

    % add legends
    xlabel(ctq_scores{iCTQ});
    ylabel('eNAMPT (ng/mL)');
    legend_size(30);
    ylim([0 30]);
    if ~strcmp(ctq_field, 'md')
        xlim([5 25]);
    end
    
    % do some stats
    % 1) remove NaN
    goodSubs = (~isnan(CTQ_var) & ~isnan(eNAMPT));
    % do the stats
    CTQ_vs_eNAMPT_stats.(ctq_field).corrcoeff = corr2(CTQ_var(goodSubs), eNAMPT(goodSubs));
    [CTQ_vs_eNAMPT_stats.(ctq_field).betas,...
        ~,...
        CTQ_vs_eNAMPT_stats.(ctq_field).stats] = glmfit(CTQ_var(goodSubs),...
        eNAMPT(goodSubs),'normal');
%     [CTQ_vs_eNAMPT_stats.(ctq_field).betas,...
%         ~,...
%         CTQ_vs_eNAMPT_stats.(ctq_field).stats] = glmfit(CTQ_var(goodSubs),...
%         eNAMPT(goodSubs),'normal','constant','off');
end