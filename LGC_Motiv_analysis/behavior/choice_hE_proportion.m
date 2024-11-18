function[choice_hE] = choice_hE_proportion(study_nm, condition, subject_id, fig_disp)
% [choice_hE] = choice_hE_proportion(study_nm, condition, subject_id, fig_disp)
% choice_hE_proportion will extract the proportion of effortful choices per
%individual and per task (physical/mental) and also across the two tasks.
%
% INPUTS
% study_nm: study name 'study1'/'study2'
%
% condition: condition to use for filtering subjects and runs
%
% subject_id: subject identification list
%
% fig_disp: figure display (1) or not (0)?
%
% OUTPUTS
% choice_hE: structure with information regarding average across
% participants and individual proportion of high effort choices per task
% and across tasks

%% define subjects to include
if ~exist('study_nm','var') || ~ismember(study_nm,{'study1','study2'})
    study_nm = 'study1';
end
if ~exist('condition','var') || isempty(condition)
    condition = subject_condition();
end
if ~exist('subject_id','var') || isempty(subject_id)
    [subject_id, NS] = LGCM_subject_selection(study_nm, condition, 'all');
else
    NS = length(subject_id);
end

%% working directories
which_pc = 'home';
switch which_pc
    case 'lab'
        rootPath = [fullfile('E:',study_nm),filesep];
    case 'home'
        serverRoot = fullfile(filesep,filesep,'sv-nas1.rcp.epfl.ch',filesep,'sandi-lab');
        rootPath = [fullfile(serverRoot,'human_data_private','raw_data_subject',study_nm),filesep];
end

%% display figure = yes by default
if ~exist('fig_disp','var') || ~ismember(fig_disp,[0,1])
    fig_disp = 1;
end
    
%% extract the data
figDispIndiv = 0;
% initialize variables of interest
[choice_hE.Ep,...
    choice_hE.Em,...
    choice_hE.EpEm] = deal(NaN(1,NS));
for iS = 1:NS
    sub_nm = subject_id{iS};
    sub_folder = [rootPath,'CID',sub_nm,filesep,'behavior'];
    % select ok runs
    runs = runs_definition(study_nm, sub_nm, condition);
    runs_ok_Ep_tmp = run_conversion(runs.Ep.runsToKeep);
    runs_ok_Em_tmp = run_conversion(runs.Em.runsToKeep);
    
    [choiceND_percentage_perRun_tmp,...
        ~] = choiceNDproportion_perRun(sub_nm, figDispIndiv, sub_folder);
    
    %% physical task
    if ~isempty(runs_ok_Ep_tmp)
        if ismember(1,runs_ok_Ep_tmp)
            run1_Ep = choiceND_percentage_perRun_tmp.Ep.run1;
        else
            run1_Ep = [];
        end
        if ismember(2,runs_ok_Ep_tmp)
            run2_Ep = choiceND_percentage_perRun_tmp.Ep.run2;
        else
            run2_Ep = [];
        end
        choice_hE.Ep(iS) = mean([run1_Ep, run2_Ep],2,'omitnan');
    end % filter case where no physical effort run
    
    %% mental task
    if ~isempty(runs_ok_Ep_tmp)
        if ismember(1,runs_ok_Em_tmp)
            run1_Em = choiceND_percentage_perRun_tmp.Em.run1;
        else
            run1_Em = [];
        end
        if ismember(2,runs_ok_Ep_tmp)
            run2_Em = choiceND_percentage_perRun_tmp.Em.run2;
        else
            run2_Em = [];
        end
        choice_hE.Em(iS) = mean([run1_Em, run2_Em],2,'omitnan');
    end % filter case where no mental effort run
    
    %% average across tasks
    if ~isempty(runs.runsToKeep)
        if ismember(1, runs.runsToKeep)
            run1_data = choiceND_percentage_perRun_tmp.run1;
        else
            run1_data = [];
        end
        if ismember(2, runs.runsToKeep)
            run2_data = choiceND_percentage_perRun_tmp.run2;
        else
            run2_data = [];
        end
        if ismember(3, runs.runsToKeep)
            run3_data = choiceND_percentage_perRun_tmp.run3;
        else
            run3_data = [];
        end
        if ismember(4, runs.runsToKeep)
            run4_data = choiceND_percentage_perRun_tmp.run4;
        else
            run4_data = [];
        end
        choice_hE.EpEm(iS) = mean([run1_data,...
            run2_data,...
            run3_data,...
            run4_data],2,'omitnan');
    end % filter case where no run is to be kept (leave NaN for this subject)
end % subject loop

%% average and SEM
[choice_hE.mean.Ep, choice_hE.sem.Ep] = mean_sem_sd(choice_hE.Ep,2);
[choice_hE.mean.Em, choice_hE.sem.Em] = mean_sem_sd(choice_hE.Em,2);
[choice_hE.mean.EpEm, choice_hE.sem.EpEm] = mean_sem_sd(choice_hE.EpEm,2);

%% figure
if fig_disp == 1
    Ep_col = [0 1 0];
    Em_col = [0 1 0];
    EpEm_col = [0 0 1];
    pSize = 30;
    
    fig;
    violinplot([choice_hE.Ep', choice_hE.Em',...
        choice_hE.EpEm'],...
        {'Ep','Em','EpEm'},...
        'ViolinColor',[Ep_col;Em_col;EpEm_col]);
    ylabel('High E choices (%)');
    xlabel('Task');
    legend_size(pSize);
end % figure

end % function