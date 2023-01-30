% check correlation between decrease of performance in max perf
% before/after each run and fatigue parameter

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
% extract parameter of interest
prmToTest = {'kFp','kEp'};
nPrmToTest = length(prmToTest);

%% load fatigue before/after each run
figDisp = 0;
maxPerf = NaN(4, NS);
for iS = 1:NS
    sub_nm = subject_id{iS};
     maxPerf_tmp = maxPerfEvolutionAcrossRuns(computerRoot, study_nm, sub_nm, figDisp);
     maxPerf(:,iS) = maxPerf_tmp.Ep;
end % subject loop
maxPerfdelta_run1 = maxPerf(2,:) - maxPerf(1,:);
maxPerfdelta_run2 =  maxPerf(4,:) - maxPerf(3,:);
maxPerfdelta_avg =  mean([maxPerfdelta_run1; maxPerfdelta_run2],1,'omitnan');
maxPerfdelta_run2_min_run1 =  maxPerf(4,:) - maxPerf(1,:);

for iPrm = 1:nPrmToTest
    prm_nm = prmToTest{iPrm};
    prm_of_interest = prm.(prm_nm);
    
    %% correlate for each run and averaging across runs
    goodSubs_r1 = (~isnan(prm_of_interest)).*(~isnan(maxPerfdelta_run1)) == 1;
    [beta.r1,~,stats.r1] = glmfit(prm_of_interest(goodSubs_r1), maxPerfdelta_run1(goodSubs_r1), 'normal');
    bhvPrm_sorted_r1 = sort(prm_of_interest(goodSubs_r1));
    maxPerf_delta_r1_fit = glmval(beta.r1, bhvPrm_sorted_r1, 'identity');
    
    goodSubs_r2 = (~isnan(prm_of_interest)).*(~isnan(maxPerfdelta_run2)) == 1;
    [beta.r2,~,stats.r2] = glmfit(prm_of_interest(goodSubs_r2), maxPerfdelta_run2(goodSubs_r2), 'normal');
    bhvPrm_sorted_r2 = sort(prm_of_interest(goodSubs_r2));
    maxPerf_delta_r2_fit = glmval(beta.r2, bhvPrm_sorted_r2, 'identity');
    
    goodSubs_avg = (~isnan(prm_of_interest)).*(~isnan(maxPerfdelta_avg)) == 1;
    [beta.avg,~,stats.avg] = glmfit(prm_of_interest(goodSubs_avg), maxPerfdelta_avg(goodSubs_avg), 'normal');
    bhvPrm_sorted_avg = sort(prm_of_interest(goodSubs_avg));
    maxPerf_delta_avg_fit = glmval(beta.avg, bhvPrm_sorted_avg, 'identity');
    
    goodSubs_r2_min_r1 = (~isnan(prm_of_interest)).*(~isnan(maxPerfdelta_run2_min_run1)) == 1;
    [beta.r2_min_r1,~,stats.r2_min_r1] = glmfit(prm_of_interest(goodSubs_r2_min_r1), maxPerfdelta_run2_min_run1(goodSubs_r2_min_r1), 'normal');
    bhvPrm_sorted_r2_min_r1 = sort(prm_of_interest(goodSubs_r2_min_r1));
    maxPerf_delta_r2_min_r1_fit = glmval(beta.r2_min_r1, bhvPrm_sorted_r2_min_r1, 'identity');
    
    %% display figure
    % general figure infos
    pSize = 30;
    lSize = 2;
    lWidth = 3;
    grey = [143 143 143]./255;
    
    fig;
    % run 1
    subplot(2,2,1);
    hold on;
    scat_hdl = scatter(prm_of_interest(goodSubs_r1), maxPerfdelta_run1(goodSubs_r1));
    scat_hdl.LineWidth = lWidth;
    scat_hdl.MarkerEdgeColor = 'k';
    fit_hdl = plot(bhvPrm_sorted_r1, maxPerf_delta_r1_fit);
    fit_hdl.LineWidth = lWidth;
    fit_hdl.Color = grey;
    fit_hdl.LineStyle = '--';
    xlabel(prm_nm);
    ylabel({'max perf change (%)';'run 1'})
    legend_size(pSize);
    
    % run 2
    subplot(2,2,2);
    hold on;
    scat_hdl = scatter(prm_of_interest(goodSubs_r2), maxPerfdelta_run2(goodSubs_r2));
    scat_hdl.LineWidth = lWidth;
    scat_hdl.MarkerEdgeColor = 'k';
    fit_hdl = plot(bhvPrm_sorted_r2, maxPerf_delta_r2_fit);
    fit_hdl.LineWidth = lWidth;
    fit_hdl.Color = grey;
    fit_hdl.LineStyle = '--';
    xlabel(prm_nm);
    ylabel({'max perf change (%)';'run 2'})
    legend_size(pSize);
    
    % (run 2 end)-(run1 start)
    subplot(2,2,3);
    hold on;
    scat_hdl = scatter(prm_of_interest(goodSubs_r2_min_r1), maxPerfdelta_run2_min_run1(goodSubs_r2_min_r1));
    scat_hdl.LineWidth = lWidth;
    scat_hdl.MarkerEdgeColor = 'k';
    fit_hdl = plot(bhvPrm_sorted_r2_min_r1, maxPerf_delta_r2_min_r1_fit);
    fit_hdl.LineWidth = lWidth;
    fit_hdl.Color = grey;
    fit_hdl.LineStyle = '--';
    xlabel(prm_nm);
    ylabel({'max perf change (%)';'task end - task start'})
    legend_size(pSize);
    
    % average
    subplot(2,2,4);
    hold on;
    scat_hdl = scatter(prm_of_interest(goodSubs_avg), maxPerfdelta_run1(goodSubs_avg));
    scat_hdl.LineWidth = lWidth;
    scat_hdl.MarkerEdgeColor = 'k';
    fit_hdl = plot(bhvPrm_sorted_avg, maxPerf_delta_avg_fit);
    fit_hdl.LineWidth = lWidth;
    fit_hdl.Color = grey;
    fit_hdl.LineStyle = '--';
    xlabel(prm_nm);
    ylabel({'max perf change (%)';'average 2 runs'})
    legend_size(pSize);
    
end % parameter loop