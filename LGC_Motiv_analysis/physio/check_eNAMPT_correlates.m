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