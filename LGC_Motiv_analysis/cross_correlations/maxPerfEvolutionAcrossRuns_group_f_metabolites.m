function[maxPerf, pval] = maxPerfEvolutionAcrossRuns_group_f_metabolites(computerRoot, study_nm, figGroupDisp, figIndivDisp)
% [maxPerf, pval] = maxPerfEvolutionAcrossRuns_group_f_metabolites(computerRoot, study_nm, figGroupDisp, figIndivDisp)
% maxPerfEvolutionAcrossRuns_group_f_metabolites will look at maximal
% performance of participants before and after each run for each task
% depending on their metabolite levels (based on the metabolite selected by
% the user).
%
% INPUTS
% computerRoot: pathway where data is
% 
% study_nm: study name
%
% sub_nm: subject number id 'XXX'
%
% figGroupDisp: display group figure (1) or not (0)
%
% figIndivDisp: display individual figure (1) or not (0)
%
% OUTPUTS
% maxPerf: structure with maximal performance individually and averaged
%
% pval: structure with p.value for different tests     

%% if root not defined => ask for it
if ~exist('computerRoot','var') || isempty(computerRoot)
    computerRoot = LGCM_root_paths;
end

%% study names
if ~exist('study_nm','var') || isempty(study_nm)
    study_nm = 'study1';
end

%% working directories
studyBehaviorFolder = [computerRoot, filesep, study_nm, filesep];
resultFolder_a = [studyBehaviorFolder,'results',filesep];
if ~exist(resultFolder_a,'dir')
    mkdir(resultFolder_a);
end
resultFolder = [resultFolder_a,'figures',filesep];
if ~exist(resultFolder,'dir')
    mkdir(resultFolder);
end

%% subject selection
[subject_id, NS] = LGCM_subject_selection(study_nm, 'behavior');

%% by default, display group figure
if ~exist('figGroupDisp','var') || isempty(figGroupDisp)
    figGroupDisp = 1;
    disp(['figGroupDisp was not defined in the inputs so that by default ',...
        'figures are displayed for the group.']);
end

%% by default, do not display individual figure
if ~exist('figIndivDisp','var') || isempty(figIndivDisp)
    figIndivDisp = 0;
    disp(['figGroupDisp was not defined in the inputs so that by default ',...
        'figures are displayed for each individual.']);
end

%% define metabolite and ROI you want to focus on
[low_met_subs, high_met_subs,...
    metabolite_nm, MRS_ROI_nm, metabolite_allSubs] = medSplit_metabolites(study_nm, subject_id);

%% initialize variables of interest
n_maxPerf = 4;
[maxPerf.Ep.allData, maxPerf.Em.allData] = deal(NaN(n_maxPerf, NS));
%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    
    % extract individual data
    [maxPerf_tmp] = maxPerfEvolutionAcrossRuns(computerRoot, study_nm, sub_nm, figIndivDisp);
    
    for iPM = 1:2
        switch iPM
            case 1
                task_id = 'Ep';
                task_fullName = 'physical';
            case 2
                task_id = 'Em';
                task_fullName = 'mental';
        end
        % extract data
        maxPerf.(task_id).allData(:, iS) = maxPerf_tmp.(task_id);
        
    end % physical/mental
end % subject loop

%% average data according to metabolite levels
for iPM = 1:2
        switch iPM
            case 1
                task_id = 'Ep';
                task_fullName = 'physical';
            case 2
                task_id = 'Em';
                task_fullName = 'mental';
        end
        % average data
        [maxPerf.(task_id).mean.low,...
            maxPerf.(task_id).sem.low,...
            maxPerf.(task_id).sd.low] = mean_sem_sd(maxPerf.(task_id).allData(:,low_met_subs),2);
        [maxPerf.(task_id).mean.high,...
            maxPerf.(task_id).sem.high,...
            maxPerf.(task_id).sd.high] = mean_sem_sd(maxPerf.(task_id).allData(:,high_met_subs),2);
        % perform ttest for each group of before vs after run
        [~,pval.end_vs_start.(task_id).low] = ttest(maxPerf.(task_id).allData(1,low_met_subs),...
            maxPerf.(task_id).allData(end,low_met_subs));
        [~,pval.end_vs_start.(task_id).high] = ttest(maxPerf.(task_id).allData(1,high_met_subs),...
            maxPerf.(task_id).allData(end,high_met_subs));
%         % perform GLM to include within and between session effects
%         run_effect = [0; 0; 1; 1];
%         withinRun_effect = [0; 1; 0; 1];
%         [~,~,stats] = glmfit([run_effect, withinRun_effect], maxPerf.(task_id).mean,'normal');
%         pval.glm.(task_id).run_effect = stats.p(2);
%         pval.glm.(task_id).withinRun_effect = stats.p(3);

% perform GLM to test metabolite level main effect:
% Performance = b0 + bM*metabolite + bS
% xGLM = [metabolite_allSubs, ];
        %% display figure
        if figGroupDisp == 1
            pSize = 30;
            bWidth = 0.4;
            bDist = 0.2;
            lowColor = ([67 162 202])./255;
            highColor = [161 217 155]./255;
            % display figure
            fig;
            hold on;
            lowHdl = bar((1-bDist):(4-bDist), maxPerf.(task_id).mean.low,...
                'BarWidth',bWidth,'FaceColor',lowColor);
            highHdl = bar((1+bDist):(4+bDist), maxPerf.(task_id).mean.high,...
                'BarWidth',bWidth,'FaceColor',highColor);
            errorbarHdl = errorbar((1-bDist):(4-bDist),...
                maxPerf.(task_id).mean.low, maxPerf.(task_id).sem.low);
            errorbarHdl.LineStyle = 'none';
            errorbarHdl.LineWidth = 3;
            errorbarHdl.Color = [0 0 0];
            errorbarHdl = errorbar((1+bDist):(4+bDist),...
                maxPerf.(task_id).mean.high, maxPerf.(task_id).sem.high);
            errorbarHdl.LineStyle = 'none';
            errorbarHdl.LineWidth = 3;
            errorbarHdl.Color = [0 0 0];
            legend([lowHdl, highHdl],...
                {['low ',metabolite_nm],...
                ['high ',metabolite_nm]},...
                'Location','NorthWest');
            legend('boxoff');
            ylim([60 120]);
            ylabel(['Max Performance ',task_fullName, ' (%)']);
            xticks(1:4);
            xticklabels({'pre_R_1','post_R_1','pre_R_2','post_R_2'})
            legend_size(pSize);
        end % figure display

end % physical/mental

end % function