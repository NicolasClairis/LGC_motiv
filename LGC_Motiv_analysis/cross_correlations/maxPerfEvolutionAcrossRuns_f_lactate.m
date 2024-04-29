function[betas, pval, r_corr] = maxPerfEvolutionAcrossRuns_f_lactate(fig_disp, rmv_outliers_yn)
% [betas, pval, rho] = maxPerfEvolutionAcrossRuns_f_lactate(fig_disp, rmv_outliers_yn)
% maxPerfEvolutionAcrossRuns_f_lactate will compare plasma and brain levels of lactate to
% maximal performance before/after each run.
%
% INPUTS
% fig_disp: display figure (1) or not (0)? By default will display it
%
% rmv_outliers_yn: remove median +/- 3*SD outliers yes (1) or no (0)? Yes by
% default
%
% OUTPUTS
% betas: regression estimates for correlations
%
% pval: structure with p.value for correlations
%
% rho: correlation coefficient for correlations

%% default inputs
if ~exist('fig_disp','var') || isempty(fig_disp) || ~ismember(fig_disp,[0,1])
    fig_disp = 1;
end
% remove outliers by default
if ~exist('rmv_outliers_yn','var') || isempty(rmv_outliers_yn) || ~ismember(rmv_outliers_yn,[0,1])
    rmv_outliers_yn = 1;
end

%% subject selection
[study_nm, ~, ~, subject_id, NS] = sub_id;

%% working directory
if ~exist('computerRoot','var') || isempty(computerRoot)
    computerRoot = LGCM_root_paths;
end

%% initialize variable of interest
[plasma_Lac, dmPFC_Lac, aIns_Lac] = deal(NaN(1,NS));
n_maxPerf = 4;
[maxPerf.Ep, maxPerf.Em] = deal(NaN(n_maxPerf, NS));
n_sessionsPerTask = 2;
[deltaMaxPerf.Ep, deltaMaxPerf.Em] = deal(NaN(n_sessionsPerTask, NS));

%% load max perf for each runs and each task
%loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    
    % extract individual data
    [maxPerf_tmp] = maxPerfEvolutionAcrossRuns(computerRoot, study_nm, sub_nm, 0);
    
    for iPM = 1:2
        switch iPM
            case 1
                task_id = 'Ep';
            case 2
                task_id = 'Em';
        end
        % extract data
        maxPerf.(task_id)(:, iS) = maxPerf_tmp.(task_id);
        deltaMaxPerf.(task_id)(1, iS) = maxPerf.(task_id)(2, iS) - maxPerf.(task_id)(1, iS);
        deltaMaxPerf.(task_id)(2, iS) = maxPerf.(task_id)(4, iS) - maxPerf.(task_id)(3, iS);
    end % loop over physical/mental task
end % subject loop

%% load brain metabolites
[metabolites] = metabolite_load(subject_id);
dmPFC_Lac(:) = metabolites.dmPFC.Lac;
aIns_Lac(:) = metabolites.aIns.Lac;

%% load plasma lactate
[plasma_Lac_struct] = load_plasma_Lac(subject_id);
plasma_Lac(:) = plasma_Lac_struct.Lac;

%% remove outliers
if rmv_outliers_yn == 1
    for iPM = 1:2
        switch iPM
            case 1
                task_id = 'Ep';
            case 2
                task_id = 'Em';
        end
        [~,~,maxPerf.(task_id)(1,:)] = rmv_outliers_3sd(maxPerf.(task_id)(1,:));
        [~,~,maxPerf.(task_id)(2,:)] = rmv_outliers_3sd(maxPerf.(task_id)(2,:));
        [~,~,maxPerf.(task_id)(3,:)] = rmv_outliers_3sd(maxPerf.(task_id)(3,:));
        [~,~,maxPerf.(task_id)(4,:)] = rmv_outliers_3sd(maxPerf.(task_id)(4,:));
        [~,~,deltaMaxPerf.(task_id)(1,:)] = rmv_outliers_3sd(deltaMaxPerf.(task_id)(1,:));
        [~,~,deltaMaxPerf.(task_id)(2,:)] = rmv_outliers_3sd(deltaMaxPerf.(task_id)(2,:));
    end % loop physical/mental task
    [~,~,dmPFC_Lac] = rmv_outliers_3sd(dmPFC_Lac);
    [~,~,aIns_Lac] = rmv_outliers_3sd(aIns_Lac);
    [~,~,plasma_Lac] = rmv_outliers_3sd(plasma_Lac);
end % remove outliers

%% test correlations
for iPM = 1:2
    switch iPM
        case 1
            task_id = 'Ep';
        case 2
            task_id = 'Em';
    end
    %% test cor max performance for each run
    for iTimePoint = 1:n_maxPerf
        timePoint_nm = ['timePoint_',num2str(iTimePoint)];
        [r_corr.maxPerf.(task_id).plasma_Lac.(timePoint_nm),...
            betas.maxPerf.(task_id).plasma_Lac.(timePoint_nm),...
            pval.maxPerf.(task_id).plasma_Lac.(timePoint_nm),...
            ~, plasma_Lac_sorted,...
            maxPerf_fit.maxPerf.(task_id).plasma_Lac.(timePoint_nm)] = glm_package(plasma_Lac, maxPerf.(task_id)(iTimePoint,:), 'normal');
        
        [r_corr.maxPerf.(task_id).dmPFC_Lac.(timePoint_nm),...
            betas.maxPerf.(task_id).dmPFC_Lac.(timePoint_nm),...
            pval.maxPerf.(task_id).dmPFC_Lac.(timePoint_nm),...
            ~, dmPFC_Lac_sorted,...
            maxPerf_fit.maxPerf.dmPFC_Lac.(task_id).(timePoint_nm)] = glm_package(dmPFC_Lac, maxPerf.(task_id)(iTimePoint,:), 'normal');
        
        [r_corr.maxPerf.(task_id).aIns_Lac.(timePoint_nm),...
            betas.maxPerf.(task_id).aIns_Lac.(timePoint_nm),...
            pval.maxPerf.(task_id).aIns_Lac.(timePoint_nm),...
            ~, aIns_Lac_sorted,...
            maxPerf_fit.maxPerf.(task_id).aIns_Lac.(timePoint_nm)] = glm_package(aIns_Lac, maxPerf.(task_id)(iTimePoint,:), 'normal');
    end % loop over time points (before/after each run)
    
    %% test how the delta is affected
    for iR = 1:n_sessionsPerTask
        run_nm = ['run_',num2str(iR)];
        [r_corr.deltaMaxPerf.(task_id).plasma_Lac.(run_nm),...
            betas.deltaMaxPerf.(task_id).plasma_Lac.(run_nm),...
            pval.deltaMaxPerf.(task_id).plasma_Lac.(run_nm),...
            ~, plasma_Lac_sorted,...
            deltaMaxPerf_fit.deltaMaxPerf.(task_id).plasma_Lac.(run_nm)] = glm_package(plasma_Lac, deltaMaxPerf.(task_id)(iR,:), 'normal');
        
        [r_corr.deltaMaxPerf.(task_id).dmPFC_Lac.(run_nm),...
            betas.deltaMaxPerf.(task_id).dmPFC_Lac.(run_nm),...
            pval.deltaMaxPerf.(task_id).dmPFC_Lac.(run_nm),...
            ~, dmPFC_Lac_sorted,...
            deltaMaxPerf_fit.deltaMaxPerf.(task_id).dmPFC_Lac.(run_nm)] = glm_package(dmPFC_Lac, deltaMaxPerf.(task_id)(iR,:), 'normal');
        
        [r_corr.deltaMaxPerf.(task_id).aIns_Lac.(run_nm),...
            betas.deltaMaxPerf.(task_id).aIns_Lac.(run_nm),...
            pval.deltaMaxPerf.(task_id).aIns_Lac.(run_nm),...
            ~, aIns_Lac_sorted,...
            deltaMaxPerf_fit.deltaMaxPerf.(task_id).aIns_Lac.(run_nm)] = glm_package(aIns_Lac, deltaMaxPerf.(task_id)(iR,:), 'normal');
    end % loop over sessions
end % task loop

%% figure
if fig_disp == 1
    %% aggregate data in one big matrix for the heatmap
    for iPM = 1:2
        switch iPM
            case 1
                task_id = 'Ep';
            case 2
                task_id = 'Em';
        end
        
        %% correlation matrix and pvalue with all correlations
        dmPFC_r_vector.(task_id) = [r_corr.maxPerf.(task_id).dmPFC_Lac.timePoint_1;...
            r_corr.maxPerf.(task_id).dmPFC_Lac.timePoint_2;...
            r_corr.maxPerf.(task_id).dmPFC_Lac.timePoint_3;...
            r_corr.maxPerf.(task_id).dmPFC_Lac.timePoint_4;...
            r_corr.deltaMaxPerf.(task_id).dmPFC_Lac.run_1;...
            r_corr.deltaMaxPerf.(task_id).dmPFC_Lac.run_2];
        aIns_r_vector.(task_id) = [r_corr.maxPerf.(task_id).aIns_Lac.timePoint_1;...
            r_corr.maxPerf.(task_id).aIns_Lac.timePoint_2;...
            r_corr.maxPerf.(task_id).aIns_Lac.timePoint_3;...
            r_corr.maxPerf.(task_id).aIns_Lac.timePoint_4;...
            r_corr.deltaMaxPerf.(task_id).aIns_Lac.run_1;...
            r_corr.deltaMaxPerf.(task_id).aIns_Lac.run_2];
        plasma_r_vector.(task_id) = [r_corr.maxPerf.(task_id).plasma_Lac.timePoint_1;...
            r_corr.maxPerf.(task_id).plasma_Lac.timePoint_2;...
            r_corr.maxPerf.(task_id).plasma_Lac.timePoint_3;...
            r_corr.maxPerf.(task_id).plasma_Lac.timePoint_4;...
            r_corr.deltaMaxPerf.(task_id).plasma_Lac.run_1;...
            r_corr.deltaMaxPerf.(task_id).plasma_Lac.run_2];
        corr_mtrx.(task_id) = [dmPFC_r_vector.(task_id),...
            aIns_r_vector.(task_id),...
            plasma_r_vector.(task_id)];
        
        % same with p.value
        dmPFC_pval_vector.(task_id) = [pval.maxPerf.(task_id).dmPFC_Lac.timePoint_1(2);...
            pval.maxPerf.(task_id).dmPFC_Lac.timePoint_2(2);...
            pval.maxPerf.(task_id).dmPFC_Lac.timePoint_3(2);...
            pval.maxPerf.(task_id).dmPFC_Lac.timePoint_4(2);...
            pval.deltaMaxPerf.(task_id).dmPFC_Lac.run_1(2);...
            pval.deltaMaxPerf.(task_id).dmPFC_Lac.run_2(2)];
        aIns_pval_vector.(task_id) = [pval.maxPerf.(task_id).aIns_Lac.timePoint_1(2);...
            pval.maxPerf.(task_id).aIns_Lac.timePoint_2(2);...
            pval.maxPerf.(task_id).aIns_Lac.timePoint_3(2);...
            pval.maxPerf.(task_id).aIns_Lac.timePoint_4(2);...
            pval.deltaMaxPerf.(task_id).aIns_Lac.run_1(2);...
            pval.deltaMaxPerf.(task_id).aIns_Lac.run_2(2)];
        plasma_pval_vector.(task_id) = [pval.maxPerf.(task_id).plasma_Lac.timePoint_1(2);...
            pval.maxPerf.(task_id).plasma_Lac.timePoint_2(2);...
            pval.maxPerf.(task_id).plasma_Lac.timePoint_3(2);...
            pval.maxPerf.(task_id).plasma_Lac.timePoint_4(2);...
            pval.deltaMaxPerf.(task_id).plasma_Lac.run_1(2);...
            pval.deltaMaxPerf.(task_id).plasma_Lac.run_2(2)];
        pval_mtrx.(task_id) = [dmPFC_pval_vector.(task_id),...
            aIns_pval_vector.(task_id),...
            plasma_pval_vector.(task_id)];
        
        %% simplified correlation matrix and pvalue focused on delta
        dmPFC_deltaMP_r_vector.(task_id) = [r_corr.deltaMaxPerf.(task_id).dmPFC_Lac.run_1;...
            r_corr.deltaMaxPerf.(task_id).dmPFC_Lac.run_2];
        aIns_deltaMP_r_vector.(task_id) = [r_corr.deltaMaxPerf.(task_id).aIns_Lac.run_1;...
            r_corr.deltaMaxPerf.(task_id).aIns_Lac.run_2];
        plasma_deltaMP_r_vector.(task_id) = [r_corr.deltaMaxPerf.(task_id).plasma_Lac.run_1;...
            r_corr.deltaMaxPerf.(task_id).plasma_Lac.run_2];
        corr_mtrx_deltaMP.(task_id) = [dmPFC_deltaMP_r_vector.(task_id),...
            aIns_deltaMP_r_vector.(task_id),...
            plasma_deltaMP_r_vector.(task_id)];
        
        % same with p.value
        dmPFC_deltaMP_pval_vector.(task_id) = [pval.deltaMaxPerf.(task_id).dmPFC_Lac.run_1(2);...
            pval.deltaMaxPerf.(task_id).dmPFC_Lac.run_2(2)];
        aIns_deltaMP_pval_vector.(task_id) = [pval.deltaMaxPerf.(task_id).aIns_Lac.run_1(2);...
            pval.deltaMaxPerf.(task_id).aIns_Lac.run_2(2)];
        plasma_deltaMP_pval_vector.(task_id) = [pval.deltaMaxPerf.(task_id).plasma_Lac.run_1(2);...
            pval.deltaMaxPerf.(task_id).plasma_Lac.run_2(2)];
        pval_mtrx_deltaMP.(task_id) = [dmPFC_deltaMP_pval_vector.(task_id),...
            aIns_deltaMP_pval_vector.(task_id),...
            plasma_deltaMP_pval_vector.(task_id)];
    end % loop over physical/mental tasks
    %% general figure parameters
    [pSize, ~, col] = general_fig_prm;
    % define which colormap you want to use (see full list here if you are not
    % happy with the selection:
    % https://ch.mathworks.com/help/matlab/ref/colormap.html)
    % color_range_choices = 'hot';
    % color_range_choices = 'turbo';
    % color_range_choices = 'jet';
    color_range_choices = redblue(45);
    
    % correlation range
    corr_range = [-0.65 0.65];
    
    % x/y-axis ratio
    ax_ratio = [1.5 1 1];
    
    %% heatmap with correlation coefficients for each maximal performance measure + for delta between end-maxPerf and start-maxPerf
    fig;
    for iPM = 1:2
        switch iPM
            case 1
                task_id = 'Ep';
            case 2
                task_id = 'Em';
        end
        subP = subplot(1,2,iPM);
        imagesc(corr_mtrx.(task_id), corr_range);
        colormap(subP, color_range_choices);
        cbar = colorbar;
        cbar.Label.String = 'r';
        xticks(1:3);
        xticklabels({'dmPFC/dACC','aIns','plasma'});
        yticks(1:size(corr_mtrx.(task_id),1));
        yticklabels({'MP1','MP2','MP3','MP4','dMP1','dMP2'});
        % add stars in the graph if some correlations are significant
        for iLac_measure = 1:size(corr_mtrx.(task_id),2)
            for iFatigueRtg = 1:size(corr_mtrx.(task_id),1)
                if pval_mtrx.(task_id)(iFatigueRtg, iLac_measure) <= 0.05
                    if pval_mtrx.(task_id)(iFatigueRtg, iLac_measure) > 0.01 && pval_mtrx.(task_id)(iFatigueRtg, iLac_measure) <= 0.05
                        pval_hdl = text(iLac_measure, iFatigueRtg, '*');
                    elseif pval_mtrx.(task_id)(iFatigueRtg, iLac_measure) > 0.001 && pval_mtrx.(task_id)(iFatigueRtg, iLac_measure) <= 0.01
                        pval_hdl = text(iLac_measure, iFatigueRtg, '**');
                    elseif pval_mtrx.(task_id)(iFatigueRtg, iLac_measure) <= 0.001
                        pval_hdl = text(iLac_measure, iFatigueRtg, '***');
                    end % p.value
                    % adjust p.value parameters
                    pval_hdl.Color = col.white;
                    pval_hdl.FontSize = 70;
                    pval_hdl.FontWeight = 'bold';
                    pval_hdl.HorizontalAlignment = 'center'; % center text on x-axis
                    pval_hdl.VerticalAlignment = 'middle'; % center text on y-axis
                end % when p.value is significant
            end % loop over Y variables
        end % loop over X variables
    end % loop over physical/mental tasks
        
        
    %% simplified heatmap with correlation coefficients for delta between end-maxPerf and start-maxPerf
    fig;
    for iPM = 1:2
        switch iPM
            case 1
                task_id = 'Ep';
            case 2
                task_id = 'Em';
        end
        subP = subplot(1,2,iPM);
        imagesc(corr_mtrx_deltaMP.(task_id), corr_range);
        colormap(subP, color_range_choices);
        cbar = colorbar;
        cbar.Label.String = 'r';
        xticks(1:3);
        xticklabels({'dmPFC/dACC','aIns','plasma'});
        yticks(1:size(corr_mtrx_deltaMP.(task_id),1));
        yticklabels({'R1','R2'});
        % add stars in the graph if some correlations are significant
        for iLac_measure = 1:size(corr_mtrx_deltaMP.(task_id),2)
            for iFatigueRtg = 1:size(corr_mtrx_deltaMP.(task_id),1)
                if pval_mtrx_deltaMP.(task_id)(iFatigueRtg, iLac_measure) <= 0.05
                    if pval_mtrx_deltaMP.(task_id)(iFatigueRtg, iLac_measure) > 0.01 && pval_mtrx_deltaMP.(task_id)(iFatigueRtg, iLac_measure) <= 0.05
                        pval_hdl = text(iLac_measure, iFatigueRtg, '*');
                    elseif pval_mtrx_deltaMP.(task_id)(iFatigueRtg, iLac_measure) > 0.001 && pval_mtrx_deltaMP.(task_id)(iFatigueRtg, iLac_measure) <= 0.01
                        pval_hdl = text(iLac_measure, iFatigueRtg, '**');
                    elseif pval_mtrx_deltaMP.(task_id)(iFatigueRtg, iLac_measure) <= 0.001
                        pval_hdl = text(iLac_measure, iFatigueRtg, '***');
                    end % p.value
                    % adjust p.value parameters
                    pval_hdl.Color = col.white;
                    pval_hdl.FontSize = 70;
                    pval_hdl.FontWeight = 'bold';
                    pval_hdl.HorizontalAlignment = 'center'; % center text on x-axis
                    pval_hdl.VerticalAlignment = 'middle'; % center text on y-axis
                end % when p.value is significant
            end % loop over Y variables
        end % loop over X variables
    end % loop over physical/mental tasks
end % figure display
end % function