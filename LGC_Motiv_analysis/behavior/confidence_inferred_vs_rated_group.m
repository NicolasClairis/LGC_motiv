function[] = confidence_inferred_vs_rated_group(study_nm, iModel, n_conf_bins, figGroupDisp)
% [] = confidence_inferred_vs_rated_group(study_nm, iModel, n_conf_bins, figGroupDisp)
% confidence_inferred_vs_rated_group will compare confidence based on
% ratings given by the subjects during the choice to the confidence
% inferred by the model "iModel" according to the formula (pChoice-0.5)^2.
%
% INPUTS
% study_nm: name of the study to look at
% 'study1': first study (dmPFC + AI)
% 'study2': second study (clinical trial)
% 
% iModel: number of the model to check
%
% n_conf_bins: number of confidence bins
%
% figGroupDisp: display group average (1) or not (0)?
%
% OUTPUTS
% betas: betas of the confidence rating = f(confidence inferred) models
%
% pval: p.value for the t.test checking if the betas are significantly
% different from zero


%% define subjects
[subject_id, NS] = LGCM_subject_selection(study_nm);

%% general parameters
% avoid displaying individual figures
figIndivDisp = 0;
% betas to extract (where test will be applied)
[beta_zero, beta_confInferred] = deal(NaN(1,NS));
% bins to extract (for the figure)
[confRated_bin,...
    confInferred_bin,...
    confRatedVsInferredFitted_bin] = deal(NaN(n_conf_bins, NS));

%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    [betas_tmp, conf_bins] = confidence_inferred_vs_rated(study_nm, sub_nm,...
        iModel, n_conf_bins, figIndivDisp);
    beta_zero(iS) = betas_tmp(1);
    beta_confInferred(iS) = betas_tmp(2);
    % extract bins
    confRated_bin(:,iS)    = conf_bins.confRated_bin;
    confInferred_bin(:,iS) = conf_bins.confInferred;
    confRatedVsInferredFitted_bin(:,iS) = conf_bins.fittedConf;
end % subject loop

%% average across subjects
[avg_confInferred, sem_confInferred] = mean_sem_sd(confInferred_bin, 2);
[avg_confRated, sem_confRated] = mean_sem_sd(confRated_bin, 2);
[avg_fittedConf] = mean_sem_sd(confRatedVsInferredFitted_bin, 2);

%% display result
if figGroupDisp == 1
    pSize = 50;
    lWidth = 3;
    
    % graphical representation
    fig;
    
    % represent rated confidence = f(inferred confidence)
    avgHdl = errorbar(avg_confInferred, avg_confRated,...
        avg_confRated-sem_confRated, avg_confRated+sem_confRated,...
        avg_confInferred-sem_confInferred, avg_confInferred+sem_confInferred);
    avgHdl.Color = 'k';
    avgHdl.LineStyle = 'none';
    avgHdl.LineWidth = lWidth;
    
    % add fitted data
    hold on;
    plotHdl = plot(avg_confInferred', avg_fittedConf);
    plotHdl.Color = 'r';
    plotHdl.LineStyle = '--';
    plotHdl.LineWidth = lWidth;
    
    % define thresholds
    xlim([-0.5 0.5]);
    ylim([-0.5 0.5]);
    
    xlabel('inferred confidence');
    ylabel('rated confidence');
    legend_size(pSize);
end % display figure