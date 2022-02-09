function[maxPerf] = maxPerfEvolutionAcrossRuns(computerRoot, study_nm, sub_nm,...
    figDisp)
% [maxPerf] = maxPerfEvolutionAcrossRuns(computerRoot, study_nm, sub_nm,...
%     figDisp)
%
% INPUTS
% computerRoot: pathway where data is
% 
% study_nm: study name
%
% sub_nm: subject number id 'XXX'
%
% figDisp: display individual figure (1) or not (0)
%
% OUTPUTS
% maxPerf: structure with maximal performance

%% if root not defined => ask for it
if ~exist('computerRoot','var') || isempty(computerRoot)
    computerRoot = LGCM_root_paths;
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
runs.Ep = strcmp(runsStruct.tasks,'Ep');
runs.Em = strcmp(runsStruct.tasks,'Em');
n_maxPerf = 4;

%% loop through physical/mental tasks
for iPM = 1:2
    switch iPM
        case 1
            task_id = 'Ep';
            task_fullName = 'physical';
        case 2
            task_id = 'Em';
            task_fullName = 'mental';
    end
    maxPerf.(task_id) = NaN(1,n_maxPerf);
    
    %% load initial maximum
    switch task_id
        case 'Ep'
            maxCalib = getfield(load([subBehaviorFolder,'CID',sub_nm,'_',task_fullName,'Calib.mat']),'MVC');
        case 'Em'
            maxCalib = getfield(load([subBehaviorFolder,'CID',sub_nm,'_',task_fullName,'Calib.mat']),'NMP');
    end
    
    %% loop through runs
    jRun = 0;
    for iRun = find(runs.(task_id))
        jRun = jRun + 1;
        % load max perf beginning and end of each run
        behaviorStruct_tmp = load([subBehaviorFolder,...
                'CID',sub_nm,'_session',num2str(iRun),'_',task_fullName,...
                '_task.mat']);
            maxPerf_startRun_idx    = 1 + 2*(jRun - 1);
            maxPerf_endRun_idx      = 2 + 2*(jRun - 1);
            switch task_id
                case 'Ep'
                    maxPerf.(task_id)(maxPerf_startRun_idx) = behaviorStruct_tmp.start_maxPerf.MVC.MVC*(100/maxCalib);
                    maxPerf.(task_id)(maxPerf_endRun_idx)   = behaviorStruct_tmp.end_maxPerf.MVC.MVC*(100/maxCalib);
                case 'Em'
                    maxPerf.(task_id)(maxPerf_startRun_idx) = behaviorStruct_tmp.start_maxPerf.n_maxPerf*(100/maxCalib);
                    maxPerf.(task_id)(maxPerf_endRun_idx) = behaviorStruct_tmp.end_maxPerf.n_maxPerf*(100/maxCalib);
            end
    end % run loop
    
    %% display figure
    if figDisp == 1
        pSize = 30;
        % display figure
        fig;
        bar(1:4, maxPerf.(task_id));
        ylim([0 130]);
        ylabel(['Max Performance ',task_fullName, ' (%)']);
        xticks(1:4);
        xticklabels({'pre_R_1','post_R_1','pre_R_2','post_R_2'})
        legend_size(pSize);
    end
end % physical/mental

end % function