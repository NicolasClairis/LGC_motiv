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
    for iTimePoint = 1:n_maxPerf
        maxPerf_f_plasma_nm = ['maxPerf_',num2str(iTimePoint),'_f_plasma_Lac'];
        [r_corr.(task_id).(maxPerf_f_plasma_nm), betas.(task_id).(maxPerf_f_plasma_nm),...
            pval.(task_id).(maxPerf_f_plasma_nm),...
            ~, plasma_Lac_sorted,...
            maxPerf_fit.(task_id).(maxPerf_f_plasma_nm)] = glm_package(plasma_Lac, maxPerf.(task_id)(iTimePoint,:), 'normal');
        
        maxPerf_f_dmPFC_nm = ['maxPerf_',num2str(iTimePoint),'_f_dmPFC_Lac'];
        [r_corr.(task_id).(maxPerf_f_dmPFC_nm), betas.(task_id).(maxPerf_f_dmPFC_nm),...
            pval.(task_id).(maxPerf_f_dmPFC_nm),...
            ~, dmPFC_Lac_sorted,...
            maxPerf_fit.(task_id).(maxPerf_f_dmPFC_nm)] = glm_package(dmPFC_Lac, maxPerf.(task_id)(iTimePoint,:), 'normal');
        
        maxPerf_f_aIns_nm = ['maxPerf_',num2str(iTimePoint),'_f_aIns_Lac'];
        [r_corr.(task_id).(maxPerf_f_aIns_nm), betas.(task_id).(maxPerf_f_aIns_nm),...
            pval.(task_id).(maxPerf_f_aIns_nm),...
            ~, aIns_Lac_sorted,...
            maxPerf_fit.(task_id).(maxPerf_f_aIns_nm)] = glm_package(aIns_Lac, maxPerf.(task_id)(iTimePoint,:), 'normal');
    end % loop over time points (before/after each run)

end % task loop

end % function