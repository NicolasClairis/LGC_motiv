function [betas, pval] = metabolite_f_choice_hE(study_nm)
% [betas, pval] = metabolite_f_choice_hE(study_nm)
% metabolite_f_choice_hE will look whether metabolite levels are correlated
% to the proportion of high effort choices.
%
% INPUTS
% study_nm: study name ('study1'/'study2') will pick study 1 by default if
% not defined
%
% OUTPUTS
% betas: structure with beta for the slope metabolite=f(choice %) across
% efforts and split per high effort level proposed
%
% pval: structure with p.value for the test corresponding to each test
%

%% subject selection
if ~exist('study_nm','var') || isempty(study_nm)
    study_nm = 'study1';
end
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% working directory
computerRoot = LGCM_root_paths;
studyBehaviorFolder = [computerRoot, filesep, study_nm, filesep];

%% main parameters
tasks = {'Ep','Em'};
nTasks = length(tasks);
nRunsPerTask = 2;
n_hE_levels = 3;

%% initialize figures
% figure metabolite = f(all effortful choices %)
mb_f_E_hdl = fig;
% figure metabolite = f(effortful choices % split per effort level)
mb_f_E_perE_hdl = fig;
% color gradient for effort
col.E1 = [229 245 249]./255;
col.E2 = [153 216 201]./255;
col.E3 = [44 162 95]./255;

%% load metabolites
[metabolite_allSubs,...
    MRS_ROI_nm, metabolite_nm] = metabolite_extraction(study_nm, subject_id);
goodSubs = ~isnan(metabolite_allSubs);

for iT = 1:nTasks
     task_nm = tasks{iT};
     
     switch task_nm
         case 'Ep'
             task_fullName = 'physical';
         case 'Em'
             task_fullName = 'mental';
         otherwise
             error('error with task name');
     end
     
     [choice.(task_nm).E,...
         choice.(task_nm).E1,...
         choice.(task_nm).E2,...
         choice.(task_nm).E3] = deal(zeros(1,NS));
     
     [choice_perSub.(task_nm).E,...
         choice_perSub.(task_nm).E1,...
         choice_perSub.(task_nm).E2,...
         choice_perSub.(task_nm).E3] = deal(NaN(nRunsPerTask,NS));
     
     for iS = 1:NS
         sub_nm = subject_id{iS};
         subBehaviorFolder = [studyBehaviorFolder, 'CID',sub_nm, filesep, 'behavior',filesep];
         runsStruct = runs_definition(study_nm, sub_nm, condition);
         okRuns = runsStruct.runsToKeep;
         subRunTaskNames = runsStruct.tasks;
         
         jRun = 0;
         for iRun = 1:length(okRuns)
             kRun = okRuns(iRun);
             run_nm = num2str(kRun);
             if strcmp(subRunTaskNames{iRun},task_nm)
                 jRun = jRun + 1;
                 [choice_highE_tmp] = extract_choice_hE(subBehaviorFolder, sub_nm, run_nm, task_fullName);
                 [hE_level_tmp] = extract_hE_level(subBehaviorFolder, sub_nm, run_nm, task_fullName);
                 
                 % extract relevant info
                 choice_perSub.(task_nm).E(jRun,iS) = sum(choice_highE_tmp == 1)./sum(ismember(choice_highE_tmp,[0,1]));
                 for iE = 1:n_hE_levels
                     E_nm = ['E',num2str(iE)];
                     E_idx = hE_level_tmp == iE;
                     choice_perSub.(task_nm).(E_nm)(jRun,iS) = sum(choice_highE_tmp(E_idx) == 1)./sum(ismember(choice_highE_tmp(E_idx),[0,1]));
                 end % loop on high effort level
             end % filter task type
         end % run loop
         
         % average % of choices per subject
         choice.(task_nm).E(1,iS) = mean(choice_perSub.(task_nm).E(:,iS),1,'omitnan');
         choice.(task_nm).E1(1,iS) = mean(choice_perSub.(task_nm).E1(:,iS),1,'omitnan');
         choice.(task_nm).E2(1,iS) = mean(choice_perSub.(task_nm).E2(:,iS),1,'omitnan');
         choice.(task_nm).E3(1,iS) = mean(choice_perSub.(task_nm).E3(:,iS),1,'omitnan');
     end % subject loop
     
     %% test across subjects
     % metabolites = f(choice all E)
     [betas.(task_nm).mb_f_c_E, ~, stats_tmp] = glmfit(choice.(task_nm).E(goodSubs), metabolite_allSubs(goodSubs),'normal');
     pval.(task_nm).mb_f_c_E = stats_tmp.p;
     choice_sorted_tmp.E = sort(choice.(task_nm).E(goodSubs));
     mb_f_choice_fit.E = glmval(betas.(task_nm).mb_f_c_E, choice_sorted_tmp.E, 'identity');
     
     % metabolites = f(choice per E)
     for iE = 1:n_hE_levels
         E_nm = ['E',num2str(iE)];
         [betas.(task_nm).(['mb_f_c_',E_nm]), ~, stats_tmp] = glmfit(choice.(task_nm).(E_nm)(goodSubs), metabolite_allSubs(goodSubs),'normal');
         pval.(task_nm).(['mb_f_c_',E_nm]) = stats_tmp.p;
         choice_sorted_tmp.(E_nm) = sort(choice.(task_nm).(E_nm)(goodSubs));
         mb_f_choice_fit.(E_nm) = glmval(betas.(task_nm).(['mb_f_c_',E_nm]), choice_sorted_tmp.(E_nm), 'identity');
     end % effort loop
     
     %% figure
     pSize = 30;
     lWidth = 3;
     grey = [143 143 143]./255;
     %% figure metabolite = f(all effortful choices %)
     figure(mb_f_E_hdl);
     % metabolite =  f(all choices %)
     subplot(1,2,iT);
     hold on;
     scat_hdl = scatter(choice.(task_nm).E, metabolite_allSubs);
     fit_hdl = plot(choice_sorted_tmp.E, mb_f_choice_fit.E);
     scat_hdl.LineWidth = lWidth;
     scat_hdl.MarkerEdgeColor = 'k';
     fit_hdl.LineWidth = lWidth;
     fit_hdl.LineStyle = '--';
     fit_hdl.Color = grey;
     % mark 50% choices
     line([0.5 0.5],ylim(),...
         'Color','k','LineStyle','-',...
         'LineWidth',lWidth);
     xlabel(['Choices (%) - ',task_nm]);
     ylabel([MRS_ROI_nm,' ', metabolite_nm]);
     legend_size(pSize);
     
     
     %% figure metabolite = f(effortful choices % split per effort level)
     figure(mb_f_E_perE_hdl);
     
     % metabolite =  f(effortful choices % split per effort level)
     for iE = 1:n_hE_levels
         E_nm = ['E',num2str(iE)];
         subplot(3,2,(iT + nTasks*(iE - 1)));
         hold on;
         scat_perE_hdl.(E_nm) = scatter(choice.(task_nm).(E_nm), metabolite_allSubs);
         fit_perE_hdl.(E_nm) = plot(choice_sorted_tmp.(E_nm), mb_f_choice_fit.(E_nm));
         scat_perE_hdl.(E_nm).MarkerEdgeColor = col.(E_nm);
         scat_perE_hdl.(E_nm).LineWidth = lWidth;
         fit_perE_hdl.(E_nm).LineWidth = lWidth;
         fit_perE_hdl.(E_nm).LineStyle = '--';
         fit_perE_hdl.(E_nm).Color = col.(E_nm);
         % mark 50% choices
         line([0.5 0.5],ylim(),...
             'Color','k','LineStyle','-',...
             'LineWidth',lWidth);
         xlabel([E_nm,' choices (%) - ',task_nm]);
         ylabel([MRS_ROI_nm,' ', metabolite_nm]);
         legend_size(pSize);
     end % effort level loop

end % task loop

end % function