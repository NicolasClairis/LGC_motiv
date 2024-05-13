function[] = reorganize_bayesian_mdl_data()
% reorganize_bayesian_mdl_data will reorganize the data in files similar to
% what was done previously in order to avoid having to redo all the scripts
% which are based on the model.
% In particular, reorganize_bayesian_mdl_data will create a
% bayesian_prm.mat, a bayesian_pChoice_data.mat, and a
% bayesian_deltaNV_data.mat file which can further be used by the other
% scripts.

%% define subjects
study_nm = 'study1';
condition = 'behavior_noSatTaskSub_noSatRun';
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% working directory
root_folder = fullfile('C:','Users','clairis','Desktop');
data_folder = [fullfile(root_folder,'GitHub','LGC_motiv',...
    'LGC_Motiv_results',study_nm,'bayesian_modeling'),filesep];

%% identify the data already produced
result_file_names = ls([data_folder,'bayesian_model_*_results.mat']);
n_models = size(result_file_names,1);

%% loop through models and extract the data to produce a summary function
for iMdl = 1:n_models
    result_file_nm = result_file_names(iMdl,:);
    
    %% extract actual model number
    mdl_nm = strrep(result_file_nm,'bayesian_model_','');
    mdl_nm = strrep(mdl_nm,'_results.mat','');
    mdl_fullnm = ['mdl_',mdl_nm];
    
    %% extract parameters
    data_tmp = load([data_folder,result_file_nm]);
    parameter_names = fieldnames(data_tmp.prm);
    nPrm = length(parameter_names);
    bayesian_mdl.(mdl_fullnm) = 
end % model loop

%% save resulting files
save([data_folder,'behavioral_prm.mat'],'bayesian_mdl');
save([data_folder,'bayesian_pChoice_data.mat'],'bayesian_pChoice');
save([data_folder,'bayesian_deltaNV_data.mat'],'bayesian_deltaNV');

end % function