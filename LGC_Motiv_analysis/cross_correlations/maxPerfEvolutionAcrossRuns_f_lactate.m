function[betas, pval, r_corr] = maxPerfEvolutionAcrossRuns_f_lactate()
% [betas, pval, rho] = maxPerfEvolutionAcrossRuns_f_lactate()
% maxPerfEvolutionAcrossRuns_f_lactate will compare plasma and brain levels of lactate to
% maximal performance before/after each run.
%
% OUTPUTS
% betas: regression estimates for correlations
%
% pval: structure with p.value for correlations
%
% rho: correlation coefficient for correlations

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
        run_nm = ['run_',num2str(iTimePoint)];
        [r_corr.maxPerf.(task_id).plasma_Lac.(run_nm),...
            betas.maxPerf.(task_id).plasma_Lac.(run_nm),...
            pval.maxPerf.(task_id).plasma_Lac.(run_nm),...
            ~, plasma_Lac_sorted,...
            maxPerf_fit.maxPerf.(task_id).plasma_Lac.(run_nm)] = glm_package(plasma_Lac, maxPerf.(task_id)(iTimePoint,:), 'normal');
        
        [r_corr.maxPerf.(task_id).dmPFC_Lac.(run_nm),...
            betas.maxPerf.(task_id).dmPFC_Lac.(run_nm),...
            pval.maxPerf.(task_id).dmPFC_Lac.(run_nm),...
            ~, dmPFC_Lac_sorted,...
            maxPerf_fit.maxPerf.dmPFC_Lac.(task_id).(run_nm)] = glm_package(dmPFC_Lac, maxPerf.(task_id)(iTimePoint,:), 'normal');
        
        [r_corr.maxPerf.(task_id).aIns_Lac.(run_nm),...
            betas.maxPerf.(task_id).aIns_Lac.(run_nm),...
            pval.maxPerf.(task_id).aIns_Lac.(run_nm),...
            ~, aIns_Lac_sorted,...
            maxPerf_fit.maxPerf.(task_id).aIns_Lac.(run_nm)] = glm_package(aIns_Lac, maxPerf.(task_id)(iTimePoint,:), 'normal');
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

end % function