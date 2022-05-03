function [betas, conf_bins, pval] = confidence_inferred_vs_rated(study_nm, sub_nm,...
    iModel, n_conf_bins, figDisp)
% [betas, conf_bins, pval] = confidence_inferred_vs_rated(study_nm, sub_nm,...
%     iModel, n_conf_bins, figDisp)
% confidence_inferred_vs_rated will compare confidence inferred by the
% model 'iModel' for the subject sub_nm in the study study_nm and will
% eventually display the corresponding result (if figDisp = 1) 
%
% INPUTS
% study_nm: study name
% 'study1': dmPFC/aINS MRS study
%
% sub_nm: subject name
%
% iModel: number of the model to check
%
% n_conf_bins: number of bins for confidence
%
% figDisp:
% (0) don't display anything
% (1) display corresponding figure
%
% OUTPUTS
% betas: betas for the test confidence rated = f(confidence rated)
%
% conf_bins: structure with bins of confidence for inferred and rated
% confidence
%
% pval: structure with p.value for correlation between inferred and rated
% confidence

%% if root not defined => ask for it
if ~exist('computerRoot','var') || isempty(computerRoot)
    computerRoot = 'E:\';
%     computerRoot = LGCM_root_paths;
end

%% working directories
subBehaviorFolder = [computerRoot, filesep, study_nm, filesep,...
    'CID',sub_nm, filesep, 'behavior', filesep];

%% by default, display individual figure
if ~exist('figDisp','var') || isempty(figDisp)
    figDisp = 1;
    disp(['figDisp was not defined in the inputs so that by default ',...
        'figures are displayed for each individual.']);
end

%% extract runs
[runsStruct] = runs_definition(study_nm, sub_nm, 'behavior');
nRuns = length(runsStruct.tasks);
nTrialsPerRun = 54;
nTotalTrials = nRuns*nTrialsPerRun;

%% extract confidence inferred
[~, dataInferred] = logitfit_choices(computerRoot, study_nm, sub_nm,...
    0, 'levels', 6, 6);

%% extract rated and inferred confidence
[confRatedAllTrials, confInferredAllTrials] = deal(NaN(1,nTotalTrials));
for iRun = 1:nRuns
    run_nm = ['run',num2str(iRun)];
    run_task_id = runsStruct.tasks{iRun};
    switch run_task_id
        case 'Ep'
            task_fullName = 'physical';
        case 'Em'
            task_fullName = 'mental';
    end
    run_trials_idx = (1:nTrialsPerRun) + nTrialsPerRun*(iRun - 1);
    
    confInferredAllTrials(run_trials_idx) = dataInferred.confidenceFitted.(['mdl_',num2str(iModel)]).(run_nm);
    
    % load rated confidence
    behaviorStruct_tmp = load([subBehaviorFolder,...
        'CID',sub_nm,'_session',num2str(iRun),'_',task_fullName,...
        '_task.mat']);
    switch run_task_id
        case 'Ep'
            choice_LR_tmp = behaviorStruct_tmp.physicalPerf.choice;
        case 'Em'
            choice_LR_tmp = behaviorStruct_tmp.mentalE_perf.choice;
    end
    conf_tmp = choice_LR_tmp;
    conf_tmp(abs(choice_LR_tmp) == 2) = 1;
    conf_tmp(abs(choice_LR_tmp) == 1) = 0;
    conf_tmp(abs(choice_LR_tmp) == 0) = NaN;
    confRatedAllTrials(run_trials_idx) = conf_tmp;
end % run loop

%% exclude NaN trials (ie those where no choice was made)
goodTrials = ~isnan(confRatedAllTrials);
confRatedAllTrials = confRatedAllTrials(goodTrials);
confInferredAllTrials = confInferredAllTrials(goodTrials);

%% extract bins
[confRated_bin,...
    confInferred_bin] = do_bin2(confRatedAllTrials, confInferredAllTrials,...
    n_conf_bins, 0);

%% correlate inferred and rated confidence
% [betas, ~, stats] = glmfit(confInferredAllTrials',...
%     confRatedAllTrials', 'normal');
% pval = stats.p;
% fittedConf = glmval(betas, confInferred_bin', 'identity');
% xData = 0:0.001:1;
[betas, ~, stats] = glmfit(confInferredAllTrials',...
    confRatedAllTrials', 'binomial','Link','logit');
pval = stats.p;
fittedConf = glmval(betas, confInferred_bin', 'logit');

%% display result
if figDisp == 1
    pSize = 50;
    lWidth = 3;
    
    % graphical representation
    fig;
    
    % represent rated confidence = f(inferred confidence)
    scatHdl = scatter(confInferred_bin, confRated_bin);
    scatHdl.MarkerFaceColor = [143 143 143]./255;
    scatHdl.MarkerEdgeColor = scatHdl.MarkerFaceColor;
    scatHdl.SizeData = 50;
    
    % add fitted data
    hold on;
    plotHdl = plot(confInferred_bin', fittedConf);
    plotHdl.Color = 'k';
    plotHdl.LineStyle = '--';
    plotHdl.LineWidth = lWidth;
    
    % define thresholds
    xlim([0 1]);
    ylim([0 1]);
    
    ylabel('rated confidence');
    xlabel('inferred confidence');
    legend_size(pSize);
end % display figure

%% extract bins
conf_bins.confInferred = confInferred_bin;
conf_bins.confRated_bin = confRated_bin;
conf_bins.fittedConf = fittedConf;

end % function