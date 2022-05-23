%% look at correlation between BMI and eNAMPT concentrations

figure;
curveHdl = scatter(BMI,eNAMPT);
curveHdl.MarkerEdgeColor = [0 0 0];
curveHdl.LineWidth = 2;
curveHdl.SizeData = 100;
hold on;

% add thresholds for anorexia/obesity
xThresholds = [18.5, 25, 30];
for iWeight = 1:length(xThresholds)
    jWeight = xThresholds(iWeight);
    hdl = line([jWeight jWeight],[0 30]);
    hdl.LineWidth = 3; hdl.Color = [0 0 0];
end

% add legends
xlabel('BMI');
ylabel('eNAMPT (ng/mL)');
legend_size(30);


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
datesExp = NaN(1,length(eNAMPT));
NS = length(eNAMPT);
for iS = 1:NS
    datesExp(iS) = datenum(dates{iS},'DD.mm.yyyy');
end
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
    CTQ_eNAMPT_corr.(ctq_field) = corrcoef(CTQ_var,eNAMPT);
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
end