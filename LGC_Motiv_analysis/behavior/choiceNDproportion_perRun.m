function[choiceND_percentage_perRun, choiceND_perRun] = choiceNDproportion_perRun(sub_nm, figDisp, sub_folder)
%[choiceND_percentage_perRun, choiceND_perRun] = choiceNDproportion_perRun(sub_nm, figDisp, sub_folder)
% choiceNDproportion_perRun = function to check the proportion of
% non-default choices in each run
%
% INPUTS
% sub_nm: subject name
%
% figDisp: figure display (1) or not (0)
%
% sub_folder: subject behavior folder. If left empty, script will consider
% that it's the current (pwd) folder by default.
%
% OUTPUTS
% choiceND_percentage_perRun: structure with indication about the proportion of
% non-default choices compared to the total of choices done (ignoring
% trials where no choice was performed)
%
% choiceND_perRun: structure with indication about the sum of
% non-default and of default choices in the run

%% define study name
if ~exist('study_nm','var')
    %     study_names = {'study1','study2','fMRI_pilots'};
    %     study_nm_idx = listdlg('ListString',study_names);
    %     study_nm = study_names{study_nm_idx};
    study_nm = 'study1'; % by default
end

%% path check
if ~exist('sub_folder', 'var') || isempty(sub_folder)
    sub_folder = [pwd, filesep];
elseif ~strcmp(sub_folder(end),filesep)
    sub_folder = [sub_folder,filesep];
end

%% extract information about runs
[runs_task, nRuns] = runs_definition(study_nm, sub_nm, 'behavior');

[choiceND_percentage_perRun.run1, choiceND_percentage_perRun.run2,...
    choiceND_percentage_perRun.run3, choiceND_percentage_perRun.run4,...
    choiceND_perRun.run1.hEchosen, choiceND_perRun.run2.hEchosen,...
    choiceND_perRun.run3.hEchosen, choiceND_perRun.run4.hEchosen,...
    choiceND_perRun.run1.lEchosen, choiceND_perRun.run2.lEchosen,...
    choiceND_perRun.run3.lEchosen, choiceND_perRun.run4.lEchosen] = deal(NaN);
[choiceND_percentage_perRun.Ep.run1,...
    choiceND_percentage_perRun.Ep.run2,...
    choiceND_percentage_perRun.Em.run1,...
    choiceND_percentage_perRun.Em.run2,...
    choiceND_perRun.Ep.run1.hEchosen,...
    choiceND_perRun.Ep.run2.hEchosen,...
    choiceND_perRun.Em.run1.hEchosen,...
    choiceND_perRun.Em.run2.hEchosen,...
    choiceND_perRun.Ep.run1.lEchosen,...
    choiceND_perRun.Ep.run2.lEchosen,...
    choiceND_perRun.Em.run1.lEchosen,...
    choiceND_perRun.Em.run2.lEchosen] = deal(NaN);
nTrials = 54;
for iRun = 1:nRuns
    run_nm = ['run',num2str(iRun)];
    [jRun] = run_conversion(iRun);
    run_nm_bis = ['run',num2str(jRun)];
    filenm = ls([sub_folder,'CID',sub_nm,'_session',num2str(iRun),'_*_task_messyAllStuff.mat']);
    runData = load([sub_folder,filenm]);
    % side of default option
    defaultSide_tmp = runData.choiceOptions.default_LR;
    % side of chosen option
    choice_tmp = runData.perfSummary.choice;
    choiceSide_tmp = choice_tmp;
    choiceSide_tmp(choice_tmp >= 1) = 1;
    choiceSide_tmp(choice_tmp <= -1) = -1;
    % non-default (high-effort) choices
    choiceND_tmp = NaN(1,nTrials);
    jTrials = 0;
    for iTrial = 1:nTrials
        % did the participant choose the default or the non-default option
        % for the current trial?
        if choiceSide_tmp(iTrial) == -defaultSide_tmp(iTrial) % choice non-default option
            choiceND_tmp(iTrial) = 1;
            jTrials = jTrials + 1;
        elseif choiceSide_tmp(iTrial) == defaultSide_tmp(iTrial) % choice default option
            choiceND_tmp(iTrial) = 0;
            jTrials = jTrials + 1;
        else % no choice was made = ignore the trial
        end % choice = default or not?


        % % if you want to reproduce Arthur's filtering method, you also have
        % % to base the filter on the proportion of high effort*high
        % % confidence trials:
        % if choiceSide_tmp(iTrial) == -defaultSide_tmp(iTrial) && runData.perfSummary.confidence.lowOrHigh(iTrial) == 2 % choice non-default option
        %     choiceND_tmp(iTrial) = 1;
        %     jTrials = jTrials + 1;
        % elseif choiceSide_tmp(iTrial) == defaultSide_tmp(iTrial) || runData.perfSummary.confidence.lowOrHigh(iTrial) == 1 % choice default option
        %     choiceND_tmp(iTrial) = 0;
        %     jTrials = jTrials + 1;
        % else % no choice was made = ignore the trial
        % end % choice = default or not?
    end % trial loop
    
    % extract proportion of non-default choices per run (ignoring trials
    % where no choice was made)
    choiceND_percentage_perRun.(run_nm) = (sum(choiceND_tmp,'omitnan')./jTrials).*100;
    % extract proportion of non-default choices per run per task type
    choiceND_percentage_perRun.(runs_task.tasks{iRun}).(run_nm_bis) = (sum(choiceND_tmp,'omitnan')./jTrials).*100;
    
    %  extract sum for each
    choiceND_perRun.(run_nm).hEchosen = sum(choiceND_tmp == 1,'omitnan');
    choiceND_perRun.(run_nm).lEchosen = sum(choiceND_tmp == 0,'omitnan');
    choiceND_perRun.(runs_task.tasks{iRun}).(run_nm_bis).hEchosen = sum(choiceND_tmp == 1,'omitnan');
    choiceND_perRun.(runs_task.tasks{iRun}).(run_nm_bis).lEchosen = sum(choiceND_tmp == 0,'omitnan');
end % run loop

%% display figure
if figDisp == 1
    pSize = 40;
    figure;
    bar(1:4, [choiceND_percentage_perRun.run1, choiceND_percentage_perRun.run2,...
        choiceND_percentage_perRun.run3, choiceND_percentage_perRun.run4]);
    ylabel('Choice non-default option (%)');
    legend_size(pSize);
    xticks(1:4);
    switch nRuns
        case 4
            xticklabels({[runs_task.tasks{1},' 1'],...
                [runs_task.tasks{2},' 2'],...
                [runs_task.tasks{3},' 3'],...
                [runs_task.tasks{4},' 4']});
        case 3
            xticklabels({[runs_task.tasks{1},' 1'],...
                [runs_task.tasks{2},' 2'],...
                [runs_task.tasks{3},' 3']});
        case 2
            xticklabels({[runs_task.tasks{1},' 1'],...
                [runs_task.tasks{2},' 2']});
    end
    xlabel('Run');
    ylim([0 100]);
end % figure display
end % function