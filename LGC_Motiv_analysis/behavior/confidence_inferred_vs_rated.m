function [betas, conf_bins, pval] = confidence_inferred_vs_rated(study_nm, sub_nm, condition,...
    n_conf_bins, figDisp, mdlType, mdl_nm)
% [betas, conf_bins, pval] = confidence_inferred_vs_rated(study_nm, sub_nm,condition,...
%     n_conf_bins, figDisp, mdlType, mdl_nm)
% confidence_inferred_vs_rated will compare confidence inferred by the
% model for the subject sub_nm in the study study_nm and will
% eventually display the corresponding result (if figDisp = 1) 
%
% INPUTS
% study_nm: study name
% 'study1': dmPFC/aINS MRS study
%
% sub_nm: subject name
%
% condition: string indicating which runs to include for each subject
%
% n_conf_bins: number of bins for confidence
%
% figDisp:
% (0) don't display anything
% (1) display corresponding figure
%
% mdlType: indication about which model type to use 'bayesian'/'simple'
%
% mdl_nm: string in the form of 'mdl_X' containing the indication of the
% model number used
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
    computerRoot = '\\svfas5\sandi-lab\human_data_private\raw_data_subject\';
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
[runsStruct] = runs_definition(study_nm, sub_nm, condition);
nRuns = length(runsStruct.tasks);
nTrialsPerRun = 54;
nTotalTrials = nRuns*nTrialsPerRun;

%% extract confidence inferred
switch mdlType
    case 'simple'
        [~, dataInferred] = logitfit_choices(computerRoot, study_nm, sub_nm, condition,...
            0, 'levels', 6, 6);
    case 'bayesian'
        gitResultsFolder = [fullfile('C:','Users','clairis','Desktop',...
                'GitHub','LGC_motiv','LGC_Motiv_results',study_nm,'bayesian_modeling'),filesep];
        gitResultsFolder = [fullfile('C:','Users','Nicolas Clairis','Documents',...
                'GitHub','LGC_motiv','LGC_Motiv_results',study_nm,'bayesian_modeling'),filesep];
end

%% extract rated and inferred confidence
[confRatedAllTrials, confInferredAllTrials] = deal(NaN(1,nTotalTrials));
for iRun = 1:nRuns
    kRun = runsStruct.runsToKeep(iRun);
    run_nm = ['run',num2str(kRun)];
    run_task_id = runsStruct.tasks{iRun};
    switch run_task_id
        case 'Ep'
            task_fullName = 'physical';
        case 'Em'
            task_fullName = 'mental';
    end
    run_trials_idx = (1:nTrialsPerRun) + nTrialsPerRun*(iRun - 1);
    
    %% load confidence inferred from the model
    switch mdlType
        case 'simple'
            confInferredAllTrials(run_trials_idx) = dataInferred.confidenceFitted.(mdl_nm).(run_nm);
        case 'bayesian'
            [~, ~, confMdl_tmp] = extract_bayesian_mdl(gitResultsFolder, subBehaviorFolder,...
                sub_nm, num2str(kRun), task_fullName, mdl_nm);
            confInferredAllTrials(run_trials_idx) = confMdl_tmp;
    end

    %% load rated confidence
    [conf_tmp] = extract_confidence_rating(subBehaviorFolder, sub_nm, run_nm, task_fullName);
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