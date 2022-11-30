function [choiceND_perRun, saturationSubs] = choiceNDproportion_perRun_group(figGrpDisp)
%[choiceND_perRun, saturationSubs] = choiceNDproportion_perRun_group()
% choiceNDproportion_perRun_group will extract the average proportion of
% non-default choices per run across all subjects.
%
% INPUTS
% figGrpDisp : display group figure (1) or not (0)? (yes by default)
%
% OUTPUTS
% choiceND_perRun: structure with average, SD and SEM of percentage of 
% choosing the non-default option across subjects
%
% saturationSubs: structure indicating which subjects saturated and for
% which run
%   .allSat: list with all the subject who had a saturation in either side
%   (non
%   .fullDef: list with all the subjects who saturated by always getting 
%
% See also choiceNDproportion_perRun

%% subject identification
[study_nm, ~,~,subject_id, NS] = sub_id;

%% display group figure
if ~exist('figGrpDisp','var') || isempty(figGrpDisp)
    figGrpDisp = 1; % by default
end

%% working directories
computerRoot = LGCM_root_paths();
% computerRoot = ['E:',filesep];
switch study_nm
    case 'fMRI_pilots'
        studyRoot = fullfile(computerRoot,'fMRI_pilots');
    case 'study1'
        studyRoot = fullfile(computerRoot,'study1');
    case 'study2'
        studyRoot = fullfile(computerRoot,'study2');
end

%% general parameters
% display individual figures?
figIndivDisp = 0;
% prepare data to extract
[choiceND_perRun.run1, choiceND_perRun.run2,...
    choiceND_perRun.run3, choiceND_perRun.run4,...
    choiceND_perRun.Ep.run1, choiceND_perRun.Ep.run2,...
    choiceND_perRun.Em.run1, choiceND_perRun.Em.run2] = deal(NaN(1,NS));

nRuns = 4;
nRunPerTask = 2;
nTasks = 2; % physical/mental
task_names = {'Ep','Em'};
%% extract data per subject
[saturationSubs.allSat.subList,...
    saturationSubs.fullNonDef.subList,...
    saturationSubs.fullDef.subList] = deal(cell(1,NS));
for iS = 1:NS
    sub_nm = subject_id{iS};
    sub_nm_bis = ['CID',sub_nm];
    cd([studyRoot, filesep, sub_nm_bis, filesep,'behavior']);
    [choiceND_percentage_perRun_tmp] = choiceNDproportion_perRun(sub_nm, figIndivDisp);
    for iRun = 1:nRuns
        run_nm = ['run',num2str(iRun)];
        choiceND_perRun.(run_nm)(iS) = choiceND_percentage_perRun_tmp.(run_nm);
        
        % extract name of subject if subject saturated
        if choiceND_percentage_perRun_tmp.(run_nm) >= 95 % always took non-default option
            if sum(strcmp(saturationSubs.allSat.subList, sub_nm)) == 0 % add subject to the list
                saturationSubs.allSat.subList{iS} = sub_nm;
                saturationSubs.fullNonDef.subList{iS} = sub_nm;
            end
            saturationSubs.satRunsPerSub.(sub_nm_bis).(run_nm) = 'all_ND';
        elseif choiceND_percentage_perRun_tmp.(run_nm) <= 5 % always took the default option
            if sum(strcmp(saturationSubs.allSat.subList, sub_nm)) == 0 % add subject to the list
                saturationSubs.allSat.subList{iS} = sub_nm;
                saturationSubs.fullDef.subList{iS} = sub_nm;
            end
            saturationSubs.satRunsPerSub.(sub_nm_bis).(run_nm) = 'all_Def';
        end
    end % run loop
    
    for iPM = 1:nTasks
        task_nm = task_names{iPM};
        for iRun = 1:nRunPerTask
            runPerTask_nm = ['run',num2str(iRun)];
            choiceND_perRun.(task_nm).(runPerTask_nm)(iS) = choiceND_percentage_perRun_tmp.(task_nm).(runPerTask_nm);
        end % run/task loop
    end % physical/mental loop
end % subject loop
cd(studyRoot);

% clean subjects without any saturation
[subjectToRemoveAllSat,...
   subjectToRemoveNDSat,...
   subjectToRemoveDefSat] = deal(false(1,NS));
for iS = 1:NS
    subjectToRemoveAllSat(iS) = isempty(saturationSubs.allSat.subList{iS});
    subjectToRemoveDefSat(iS) = isempty(saturationSubs.fullDef.subList{iS});
    subjectToRemoveNDSat(iS) = isempty(saturationSubs.fullNonDef.subList{iS});
end
saturationSubs.allSat.subList(subjectToRemoveAllSat) = [];
saturationSubs.fullDef.subList(subjectToRemoveDefSat) = [];
saturationSubs.fullNonDef.subList(subjectToRemoveNDSat) = [];

%% average across subjects
for iRun = 1:nRuns
    run_nm = ['run',num2str(iRun)];
    [choiceND_perRun.avg.(run_nm),...
        choiceND_perRun.sem.(run_nm),...
        choiceND_perRun.std.(run_nm)] = mean_sem_sd(choiceND_perRun.(run_nm), 2);
end % run loop

for iPM = 1:nTasks
    task_nm = task_names{iPM};
    for iRun = 1:nRunPerTask
        runPerTask_nm = ['run',num2str(iRun)];
        [choiceND_perRun.avg.(task_nm).(runPerTask_nm),...
           choiceND_perRun.sem.(task_nm).(runPerTask_nm),...
            choiceND_perRun.std.(task_nm).(runPerTask_nm)] = mean_sem_sd(choiceND_perRun.(task_nm).(runPerTask_nm), 2);
    end % run/task loop
end % physical/mental loop

%% group figure
if figGrpDisp == 1
    pSize = 40;
    lSize = 3;
    
    % choice = non-default = f(run)
    fig;
    barHdl = bar(1:4,...
        [choiceND_perRun.avg.run1, choiceND_perRun.avg.run2,...
        choiceND_perRun.avg.run3, choiceND_perRun.avg.run4]);
    hold on;
    %     erHdl = errorbar(1:4,...
    %         [choiceND_perRun.avg.run1, choiceND_perRun.avg.run2,...
    %         choiceND_perRun.avg.run3, choiceND_perRun.avg.run4],...
    %         [choiceND_perRun.sem.run1, choiceND_perRun.sem.run2,...
    %         choiceND_perRun.sem.run3, choiceND_perRun.sem.run4]);
    erHdl = errorbar(1:4,...
        [choiceND_perRun.avg.run1, choiceND_perRun.avg.run2,...
        choiceND_perRun.avg.run3, choiceND_perRun.avg.run4],...
        [choiceND_perRun.std.run1, choiceND_perRun.std.run2,...
        choiceND_perRun.std.run3, choiceND_perRun.std.run4]);
    erHdl.LineStyle = 'none';
    erHdl.Color = 'k';
    erHdl.LineWidth = lSize;
    ylabel('Choice non-default option (%)');
    legend_size(pSize);
    xticks(1:4);
    xlabel('Run');
    ylim([0 100]);
    
    % choice = non-default = f(run/task)
    fig;
    barHdl = bar(1:4,...
        [choiceND_perRun.avg.Ep.run1, choiceND_perRun.avg.Ep.run2,...
        choiceND_perRun.avg.Em.run1, choiceND_perRun.avg.Em.run2]);
    hold on;
%     erHdl = errorbar(1:4,...
%         [choiceND_perRun.avg.Ep.run1, choiceND_perRun.avg.Ep.run2,...
%         choiceND_perRun.avg.Em.run1, choiceND_perRun.avg.Em.run2],...
%         [choiceND_perRun.sem.Ep.run1, choiceND_perRun.sem.Ep.run2,...
%         choiceND_perRun.sem.Em.run1, choiceND_perRun.sem.Em.run2]);
    erHdl = errorbar(1:4,...
        [choiceND_perRun.avg.Ep.run1, choiceND_perRun.avg.Ep.run2,...
        choiceND_perRun.avg.Em.run1, choiceND_perRun.avg.Em.run2],...
        [choiceND_perRun.std.Ep.run1, choiceND_perRun.std.Ep.run2,...
        choiceND_perRun.std.Em.run1, choiceND_perRun.std.Em.run2]);
    erHdl.LineStyle = 'none';
    erHdl.Color = 'k';
    erHdl.LineWidth = lSize;
    ylabel('Choice non-default option (%)');
    legend_size(pSize);
    xticks(1:4);
    xticklabels({'Ep 1','Ep 2','Em 1','Em 2'});
    xlabel('Run');
    ylim([0 100]);
end % display group figure

end % function