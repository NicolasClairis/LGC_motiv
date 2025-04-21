function[choice_hE_percentage_perRun, choice_hE_perRun] = choiceNDproportion_perRun_f_GL(sub_nm, figDisp, sub_folder)
%[choiceND_percentage_perRun, choiceND_perRun] = choiceNDproportion_perRun_f_GL(sub_nm, figDisp, sub_folder)
% choiceNDproportion_perRun_f_GL = function to check the proportion of
% non-default choices in each run separating between gain and loss trials
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
% trials where no choice was performed) splitting gains and loss trials.
%
% choiceND_perRun: structure with indication about the sum of
% non-default and of default choices in the run splitting gains and loss trials.

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

%% initialize variables of interest
% define conditions on which the loop will operate
GL = {'G','L'};
n_GL = length(GL);
run_names = {'run1','run2','run3','run4'};
nPotentialRuns = length(run_names);
tasks = {'Ep','Em'};
nTasks = length(tasks);
for iGL = 1:n_GL
    GL_nm = GL{iGL};
    for iR = 1:nPotentialRuns
        run_nm = run_names{iR};
        [choice_hE_percentage_perRun.(run_nm).(GL_nm),...
            choice_hE_perRun.(run_nm).hEchosen.(GL_nm),...
            choice_hE_perRun.(run_nm).lEchosen.(GL_nm)] = deal(NaN);
    end % run loop
    for iT = 1:nTasks
        task_nm = tasks{iT};
        [choice_hE_percentage_perRun.(task_nm).run1.(GL_nm),...
            choice_hE_percentage_perRun.(task_nm).run2.(GL_nm),...
            choice_hE_perRun.(task_nm).run1.hEchosen.(GL_nm),...
            choice_hE_perRun.(task_nm).run2.hEchosen.(GL_nm),...
            choice_hE_perRun.(task_nm).run1.lEchosen.(GL_nm),...
            choice_hE_perRun.(task_nm).run2.lEchosen.(GL_nm)] = deal(NaN);
    end % task loop
end % gain/loss loop

%% extract the proportion of high-effort choices per gain/loss condition and per session
nTotalTrials = 54; % total number of trials/run
nTrials_GL = nTotalTrials/2; % total number of gain or loss trials/run
for iRun = 1:nRuns
    run_nm = ['run',num2str(iRun)];
    [jRun] = run_conversion(iRun);
    run_nm_bis = ['run',num2str(jRun)];
    filenm = ls([sub_folder,'CID',sub_nm,'_session',num2str(iRun),'_*_task_messyAllStuff.mat']);
    runData = load([sub_folder,filenm]);
    
    % identify gain and loss trials
    G_trials = strcmp(runData.choiceOptions.R_or_P,'R');
    L_trials = strcmp(runData.choiceOptions.R_or_P,'P');
    
    for iGL = 1:n_GL
        GL_nm = GL{iGL};
        switch GL_nm
            case 'G'
                trial_idx = G_trials;
            case 'L'
                trial_idx = L_trials;
            otherwise
                error('WTF?');
        end
        
        % side of default option
        defaultSide_tmp = runData.choiceOptions.default_LR(trial_idx);
        % side of chosen option
        choice_tmp = runData.perfSummary.choice(trial_idx);
        choiceSide_tmp = choice_tmp;
        choiceSide_tmp(choice_tmp >= 1) = 1;
        choiceSide_tmp(choice_tmp <= -1) = -1;
        % high-effort choices
        choice_hE_tmp = NaN(1,nTrials_GL);
        jTrials = 0;
        for iTrial = 1:nTrials_GL
            % did the participant choose the default or the non-default option
            % for the current trial?
            if choiceSide_tmp(iTrial) == -defaultSide_tmp(iTrial) % choice non-default option
                choice_hE_tmp(iTrial) = 1;
                jTrials = jTrials + 1;
            elseif choiceSide_tmp(iTrial) == defaultSide_tmp(iTrial) % choice default option
                choice_hE_tmp(iTrial) = 0;
                jTrials = jTrials + 1;
            else % no choice was made = ignore the trial
            end % choice = default or not?
        end % trial loop
        
        % extract proportion of high-effort choices per run (ignoring trials
        % where no choice was made)
        choice_hE_percentage_perRun.(run_nm).(GL_nm) = (sum(choice_hE_tmp,'omitnan')./jTrials).*100;
        % extract proportion of non-default choices per run per task type
        choice_hE_percentage_perRun.(runs_task.tasks{iRun}).(run_nm_bis).(GL_nm) = (sum(choice_hE_tmp,'omitnan')./jTrials).*100;
        
        %  extract sum for each
        choice_hE_perRun.(run_nm).hEchosen.(GL_nm) = sum(choice_hE_tmp == 1,'omitnan');
        choice_hE_perRun.(run_nm).lEchosen.(GL_nm) = sum(choice_hE_tmp == 0,'omitnan');
        choice_hE_perRun.(runs_task.tasks{iRun}).(run_nm_bis).hEchosen.(GL_nm) = sum(choice_hE_tmp == 1,'omitnan');
        choice_hE_perRun.(runs_task.tasks{iRun}).(run_nm_bis).lEchosen.(GL_nm) = sum(choice_hE_tmp == 0,'omitnan');
    end % gain/loss loop
end % run loop

%% display figure
if figDisp == 1
    pSize = 40;
    
    figure; hold on;
    % gains
    G_hdl = bar(0.8:3.8, [choice_hE_percentage_perRun.run1.G, choice_hE_percentage_perRun.run2.G,...
        choice_hE_percentage_perRun.run3.G, choice_hE_percentage_perRun.run4.G]);
    % losses
    L_hdl = bar(1.2:4.2, [choice_hE_percentage_perRun.run1.L, choice_hE_percentage_perRun.run2.L,...
        choice_hE_percentage_perRun.run3.L, choice_hE_percentage_perRun.run4.L]);
    % color for gains vs losses
    G_hdl.FaceColor = [241 163 64]./255;
    L_hdl.FaceColor = [153 142 195]./255;
    % bar width
    G_hdl.BarWidth = 0.4;
    L_hdl.BarWidth = 0.4;
    ylabel('High Effort Choices (%)');
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
    legend([G_hdl,L_hdl],{'Gain','Loss'});
    legend('boxoff');
end % figure display
end % function