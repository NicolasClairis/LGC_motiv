function[choice_hE, r_corr, pval] = choice_hE_proportion_f_GL(study_nm, condition, subject_id, fig_disp)
% [choice_hE, r_corr, pval] = choice_hE_proportion_f_GL(study_nm, condition, subject_id, fig_disp)
% choice_hE_proportion_f_GL will extract the proportion of effortful choices per
%individual and per task (physical/mental) and also across the two tasks
%but splitting the trials between the gain and loss trials.
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
% and across tasks for gains and losses separately
%
% r_corr: correlation coefficient comparing proportion of gain to
% proportion of loss choices
%
% pval: structure with p.value for paired t.test comparing gains to losses
% as well as correlation between proportion of gain and loss high effort 
% choices

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
which_pc = 'lab';
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
GL_names = {'G','L'};
n_GL = length(GL_names);
for iGL = 1:n_GL
    GL_nm = GL_names{iGL};
    [choice_hE.Ep.(GL_nm),...
        choice_hE.Em.(GL_nm),...
        choice_hE.EpEm.(GL_nm)] = deal(NaN(1,NS));
end % gain/loss

for iS = 1:NS
    sub_nm = subject_id{iS};
    sub_folder = [rootPath,'CID',sub_nm,filesep,'behavior'];
    % select ok runs
    runs = runs_definition(study_nm, sub_nm, condition);
    runs_ok_Ep_tmp = run_conversion(runs.Ep.runsToKeep);
    runs_ok_Em_tmp = run_conversion(runs.Em.runsToKeep);
    
    [choiceND_percentage_perRun_tmp,...
        ~] = choiceNDproportion_perRun_f_GL(sub_nm, figDispIndiv, sub_folder);
    
    %% physical task
    if ~isempty(runs_ok_Ep_tmp)
        if ismember(1,runs_ok_Ep_tmp)
            run1_G_Ep = choiceND_percentage_perRun_tmp.Ep.run1.G;
            run1_L_Ep = choiceND_percentage_perRun_tmp.Ep.run1.L;
        else
            run1_G_Ep = [];
            run1_L_Ep = [];
        end
        if ismember(2,runs_ok_Ep_tmp)
            run2_G_Ep = choiceND_percentage_perRun_tmp.Ep.run2.G;
            run2_L_Ep = choiceND_percentage_perRun_tmp.Ep.run2.L;
        else
            run2_G_Ep = [];
            run2_L_Ep = [];
        end
        % pool gain sessions
        choice_hE.Ep.G(iS) = mean([run1_G_Ep, run2_G_Ep],2,'omitnan');
        % pool loss sessions
        choice_hE.Ep.L(iS) = mean([run1_L_Ep, run2_L_Ep],2,'omitnan');
    end % filter case where no physical effort run
    
    %% mental task
    if ~isempty(runs_ok_Ep_tmp)
        if ismember(1,runs_ok_Em_tmp)
            run1_G_Em = choiceND_percentage_perRun_tmp.Em.run1.G;
            run1_L_Em = choiceND_percentage_perRun_tmp.Em.run1.L;
        else
            run1_G_Em = [];
            run1_L_Em = [];
        end
        if ismember(2,runs_ok_Ep_tmp)
            run2_G_Em = choiceND_percentage_perRun_tmp.Em.run2.G;
            run2_L_Em = choiceND_percentage_perRun_tmp.Em.run2.L;
        else
            run2_G_Em = [];
            run2_L_Em = [];
        end
        % pool gain sessions
        choice_hE.Em.G(iS) = mean([run1_G_Em, run2_G_Em],2,'omitnan');
        % pool loss sessions
        choice_hE.Em.L(iS) = mean([run1_L_Em, run2_L_Em],2,'omitnan');
    end % filter case where no mental effort run
    
    %% average across tasks
    if ~isempty(runs.runsToKeep)
        if ismember(1, runs.runsToKeep)
            run1_G_data = choiceND_percentage_perRun_tmp.run1.G;
            run1_L_data = choiceND_percentage_perRun_tmp.run1.L;
        else
            run1_G_data = [];
            run1_L_data = [];
        end
        if ismember(2, runs.runsToKeep)
            run2_G_data = choiceND_percentage_perRun_tmp.run2.G;
            run2_L_data = choiceND_percentage_perRun_tmp.run2.L;
        else
            run2_G_data = [];
            run2_L_data = [];
        end
        if ismember(3, runs.runsToKeep)
            run3_G_data = choiceND_percentage_perRun_tmp.run3.G;
            run3_L_data = choiceND_percentage_perRun_tmp.run3.L;
        else
            run3_G_data = [];
            run3_L_data = [];
        end
        if ismember(4, runs.runsToKeep)
            run4_G_data = choiceND_percentage_perRun_tmp.run4.G;
            run4_L_data = choiceND_percentage_perRun_tmp.run4.L;
        else
            run4_G_data = [];
            run4_L_data = [];
        end
        choice_hE.EpEm.G(iS) = mean([run1_G_data,...
            run2_G_data,...
            run3_G_data,...
            run4_G_data],2,'omitnan');
        choice_hE.EpEm.L(iS) = mean([run1_L_data,...
            run2_L_data,...
            run3_L_data,...
            run4_L_data],2,'omitnan');
    end % filter case where no run is to be kept (leave NaN for this subject)
end % subject loop

%% average and SEM
[choice_hE.mean.Ep.G, choice_hE.sem.Ep.G] = mean_sem_sd(choice_hE.Ep.G,2);
[choice_hE.mean.Ep.L, choice_hE.sem.Ep.L] = mean_sem_sd(choice_hE.Ep.L,2);
[choice_hE.mean.Em.G, choice_hE.sem.Em.G] = mean_sem_sd(choice_hE.Em.G,2);
[choice_hE.mean.Em.L, choice_hE.sem.Em.L] = mean_sem_sd(choice_hE.Em.L,2);
[choice_hE.mean.EpEm.G, choice_hE.sem.EpEm.G] = mean_sem_sd(choice_hE.EpEm.G,2);
[choice_hE.mean.EpEm.L, choice_hE.sem.EpEm.L] = mean_sem_sd(choice_hE.EpEm.L,2);

%% compare proportion of gain and loss choices
[~,pval.ttest_G_vs_L.EpEm] = ttest(choice_hE.EpEm.G, choice_hE.EpEm.L);
[~,pval.ttest_G_vs_L.Ep] = ttest(choice_hE.Ep.G, choice_hE.Ep.L);
[~,pval.ttest_G_vs_L.Em] = ttest(choice_hE.Em.G, choice_hE.Em.L);

%% test correlations between proportion of gain and loss choices
[r_corr.EpEm, pval.r_corr.EpEm] = corr(choice_hE.EpEm.G', choice_hE.EpEm.L');
[r_corr.Ep, pval.r_corr.Ep] = corr(choice_hE.Ep.G', choice_hE.Ep.L');
[r_corr.Em, pval.r_corr.Em] = corr(choice_hE.Em.G', choice_hE.Em.L');

%% figure
if fig_disp == 1
    %     Ep_col = [0 1 0];
    %     Em_col = [0 1 0];
    %     EpEm_col = [0 0 1];
    G_col = [241 163 64]./255;
    L_col = [153 142 195]./255;
    pSize = 30;
    
    fig;
    violinplot([choice_hE.Ep.G',...
        choice_hE.Ep.L',...
        choice_hE.Em.G',...
        choice_hE.Em.L',...
        choice_hE.EpEm.G',...
        choice_hE.EpEm.L'],...
        {'G','L','G','L','G','L'},...
        'ViolinColor',...
        [G_col; L_col;...
        G_col; L_col;...
        G_col; L_col]);
    ylabel('High Effort Choices (%)');
    xlabel('Task');
    legend_size(pSize);
    xlim([0 7]);
    ylim([0 102]);
end % figure

end % function