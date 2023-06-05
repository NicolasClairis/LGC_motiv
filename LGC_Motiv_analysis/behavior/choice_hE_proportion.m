function[choice_hE] = choice_hE_proportion(study_nm, condition, fig_disp)
% [choice_hE] = choice_hE_proportion(study_nm, condition, fig_disp)
% choice_hE_proportion will extract the proportion of effortful choices per
%individual and per task (physical/mental) and also across the two tasks.
%
% INPUTS
% study_nm: study name 'study1'/'study2'
%
% condition: condition to use for filtering subjects and runs
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
[subject_id, NS] = LGCM_subject_selection(study_nm, condition, 'all');

%% working directories
which_pc = 'lab';
switch which_pc
    case 'lab'
        rootPath = [fullfile('E:',study_nm),filesep];
    case 'home'
        rootPath = [fullfile('L:','human_data_private','raw_data_subject'),filesep];
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
    cd([rootPath,'CID',sub_nm,filesep,'behavior']);
    % select ok runs
    [runs, n_runs] = runs_definition(study_nm, sub_nm, condition);
    runs_ok_Ep_tmp = run_conversion(runs.Ep.runsToKeep);
    runs_ok_Em_tmp = run_conversion(runs.Em.runsToKeep);
    
    [~,...
        choiceND_perRun_tmp] = choiceNDproportion_perRun(sub_nm, figDispIndiv);
%     switch 
    choice_hE.Ep(iS) = choiceND_perRun_tmp.Ep.run1.hEchosen +...
        choiceND_perRun_tmp.Ep.run2.hEchosen;
    choice_hE.Em(iS) = choiceND_perRun_tmp.Em.run1.hEchosen +...
        choiceND_perRun_tmp.Em.run2.hEchosen;
    choice_hE.EpEm(iS) = choiceND_perRun_tmp.Ep.run1.hEchosen +...
        choiceND_perRun_tmp.Ep.run2.hEchosen +...
        choiceND_perRun_tmp.Em.run1.hEchosen +...
        choiceND_perRun_tmp.Em.run2.hEchosen;
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
    
    fig;
    violinplot([choice_hE.Ep, choice_hE.Em,...
        choice_hE.EpEm],...
        {'Ep','Em','EpEm'},...
        'ViolinColor',[Ep_col;Em_col;EpEm_col]);
    ylabel('Nb of effortful choices');
    xlabel('Task');
end % figure

end % function