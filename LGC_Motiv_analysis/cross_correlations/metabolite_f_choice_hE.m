function [] = metabolite_f_choice_hE(study_nm)
% [] = metabolite_f_choice_hE(study_nm)
% metabolite_f_choice_hE will look whether metabolite levels are correlated
% to the proportion of high effort choices.

if ~exist('study_nm','var') || isempty(study_nm)
    study_nm = 'study1';
end

tasks = {'Ep','Em'};
nTasks = length(tasks);
% figure metabolite = f(all effortful choices %)
mb_f_E_hdl = fig;
% figure metabolite = f(effortful choices % split per effort level)
mb_f_E_perE_hdl = fig;

%% load metabolites
[metabolite_allSubs, MRS_ROI_nm, metabolite_nm] = metabolite_extraction(study_nm, subject_id);
for iT = 1:nTasks
     task_nm = tasks{iT};
     
     for iS = 1:NS
         
     end % subject loop
     
     %% average 
     
     %% figure
     pSize = 30;
     lWidth = 3;
     grey = [143 143 143]./255;
     %% figure metabolite = f(all effortful choices %)
     figure(mb_f_E_hdl);
     % metabolite =  f(all choices %)
     subplot(1,2,iT);
     hold on;
     scat_hdl = scatter(, metabolite_allSubs);
     fit_hdl = plot();
     fit_hdl.LineWidth = lWidth;
     fit_hdl.LineStyle = '-';
     fit_hdl.Color = grey;
     xlabel('Choices (%)');
     ylabel([MRS_ROI_nm,' ', metabolite_nm]);
     legend_size(pSize);
     
     
     %% figure metabolite = f(effortful choices % split per effort level)
     figure(mb_f_E_perE_hdl);
     
     % metabolite =  f(effortful choices % split per effort level)
     subplot(1,2,iT);
     hold on;
     for iE = 1:n_E_levels
         scat_perE_hdl.(['E',num2str(iE)]) = scatter(, metabolite_allSubs);
         fit_perE_hdl.(['E',num2str(iE)]) = plot();
         fit_perE_hdl.(['E',num2str(iE)]).LineWidth = lWidth;
         fit_perE_hdl.(['E',num2str(iE)]).LineStyle = '-';
         fit_perE_hdl.(['E',num2str(iE)]).Color = grey;
     end % effort level loop
     xlabel('Choices (%)');
     ylabel([MRS_ROI_nm,' ', metabolite_nm]);
     legend_size(pSize);

end % task loop

end % function