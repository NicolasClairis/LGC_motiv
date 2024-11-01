% function[condition, NS] = model_comparison()
% script to compare the different computational models. The script will
% allow to select the subjects and models to consider. It will subsequently
% compare these models using VBA model comparison VBA_groupBMC.m tool.
% The final output will show how the different models compare across
% different model comparison metrics including:
% exceedance probability and estimated model frequency derived from
% VBA_groupBMC, R², AIC, MAE and BIC.

%% define subjects to consider
[study_nm, condition, ~, subject_id, NS] = sub_id;

%% working directories
% root_path = fullfile('C:','Users','clairis','Desktop');
root_path = fullfile('C:','Users','Nicolas Clairis','Documents');
mdl_results_folder = fullfile(root_path,...
    'GitHub','LGC_motiv','LGC_Motiv_results',...
    study_nm,'bayesian_modeling');
%% select models to compare
model_names = {'1','2','3','4','5','6','7'};
model_select_idx = listdlg('ListString',model_names);
models_selected_names = model_names(model_select_idx);
n_models = length(models_selected_names);

%% prepare variables of interest
[free_energy, R2, AIC, BIC, RMSE, MAE] = deal(NaN(n_models, NS));

%% load model quality data
for iMdl = 1:n_models
    mdl_nm = models_selected_names{iMdl};
    mdl_fullnm = ['mdl_',mdl_nm];
    data_struct_tmp = load([mdl_results_folder, filesep,...
        'bayesian_model_',mdl_nm,'_results.mat']);
    mdl_quality.(mdl_fullnm) = data_struct_tmp.mdl_quality;
    subject_id_tmp = data_struct_tmp.subject_id;
    for iS = 1:NS
        sub_nm = subject_id{iS};
        sub_idx = strcmp(sub_nm, subject_id_tmp);
        if sum(sub_idx) > 0
            free_energy(iMdl,iS) = mdl_quality.(mdl_fullnm).free_energy(sub_idx);
            R2(iMdl,iS) = mdl_quality.(mdl_fullnm).R2(sub_idx);
            AIC(iMdl,iS) = mdl_quality.(mdl_fullnm).AIC(sub_idx);
            BIC(iMdl,iS) = mdl_quality.(mdl_fullnm).BIC(sub_idx);
            RMSE(iMdl,iS) = mdl_quality.(mdl_fullnm).RMSE(sub_idx);
            MAE(iMdl,iS) = mdl_quality.(mdl_fullnm).MAE(sub_idx);
        end
    end % subject loop
end % model loop

%% perform the model comparison
options.modelNames = models_selected_names;
options.DisplayWin = 1;
[posterior_BMC, out_BMC] = VBA_groupBMC(free_energy, options);

% display results in one figure
plot_model_comparison(out_BMC,AIC,MAE,BIC,R2, models_selected_names);
plot_model_comparison2(out_BMC,AIC,MAE,BIC,R2, RMSE, models_selected_names);

% end % function