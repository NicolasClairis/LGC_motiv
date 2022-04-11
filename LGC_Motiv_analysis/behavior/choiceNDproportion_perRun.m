function[choiceND_perRun] = choiceNDproportion_perRun(sub_nm, figDisp)
%[choiceND_perRun] = choiceNDproportion_perRun(sub_nm, figDisp)
% choiceNDproportion_perRun = function to check the proportion of 
% non-default choices in each run
%
% INPUTS
% sub_nm: subject name
%
% figDisp: figure display (1) or not (0)
%
% OUTPUTS
% choiceND_perRun: structure with indication about the proportion of
% non-default choices

%% define study name
if ~exist('study_nm','var') || isempty(study_nm)
    %     study_names = {'study1','study2','fMRI_pilots'};
    %     study_nm_idx = listdlg('ListString',study_names);
    %     study_nm = study_names{study_nm_idx};
    study_nm = 'study1'; % by default
end

%% extract information about runs
[runs_task, nRuns_bis] = runs_definition(study_nm, sub_nm, 'behavior');

nRuns = 4;
nTrials = 54;
for iRun = 1:nRuns
    run_nm = ['run',num2str(iRun)];
    switch iRun
        case {1, 2}
            jRun = 1;
        case {3, 4}
            jRun = 2;
    end
    run_nm_bis = ['run',num2str(jRun)];
    filenm = ls(['CID',sub_nm,'_session',num2str(iRun),'_*_task_messyAllStuff.mat']);
    runData = load(filenm);
    % side of default option
    defaultSide_tmp = runData.choiceOptions.default_LR;
    % side of chosen option
    choice_tmp = runData.perfSummary.choice;
    choiceSide_tmp = choice_tmp;
    choiceSide_tmp(choice_tmp >= 1) = 1;
    choiceSide_tmp(choice_tmp <= -1) = -1;
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
    end % trial loop
    
    % extract proportion of non-default choices per run
    choiceND_perRun.(run_nm) = (sum(choiceND_tmp,'omitnan')./jTrials).*100;
    % extract proportion of non-default choices per run per task type
    choiceND_perRun.(runs_task.tasks{iRun}).(run_nm_bis) = (sum(choiceND_tmp,'omitnan')./jTrials).*100;
end % run loop

%% display figure
if figDisp == 1
    pSize = 40;
    figure;
    bar(1:4, [choiceND_perRun.run1, choiceND_perRun.run2,...
        choiceND_perRun.run3, choiceND_perRun.run4]);
    ylabel('Choice non-default option (%)');
    legend_size(pSize);
    xticks(1:4);
    xticklabels({[runs_task.tasks{1},' 1'],...
        [runs_task.tasks{2},' 2'],...
        [runs_task.tasks{3},' 3'],...
        [runs_task.tasks{4},' 4']});
    xlabel('Run');
    ylim([0 100]);
end % figure display
end % function