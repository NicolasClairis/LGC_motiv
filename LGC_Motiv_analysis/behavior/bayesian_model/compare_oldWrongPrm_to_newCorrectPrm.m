% script to compare old parameters that were wrong (because positivity
% constraint was not properly considered to estimate the correct
% posteriors: instead of looking at mean(log(1+exp(X)) and
% var(log(1+exp(X)), we looked at log(1+exp(mean(X))) and
% log(1+exp(var(X))) which was wrong

%% folder where data is stored
data_dir = fullfile('C:','Users','Nicolas Clairis',...
    'Documents','GitHub','LGC_motiv','LGC_Motiv_results',...
    'study1','bayesian_modeling');

%% define model to focus on
mdl_nm = '5';

%% load the data
old_prm = getfield(load(fullfile(data_dir,['bayesian_model_',mdl_nm,'_results_old.mat']),'prm'),'prm'); % previous parameters (wrong)
prm = getfield(load(fullfile(data_dir,['bayesian_model_',mdl_nm,'_results.mat']),'prm'),'prm'); % current parameters (correct)

prm_names = fieldnames(prm);
for iP = 1:length(prm_names)
    prm_nm = prm_names{iP};
    
    % correlation between old and new parameter
    [r_corr.(prm_nm), pval.(prm_nm)] = corr(old_prm.(prm_nm)', prm.(prm_nm)');

    fig;
    scatter(old_prm.(prm_nm), prm.(prm_nm));
    xlabel(['wrong ',prm_nm]);
    ylabel(['true ',prm_nm]);
    hold on;
    ref_hdl = refline(1,0); % add line to show perfect ideal correlation
    ref_hdl.Color = [0.5 0.5 0.5];
    place_r_and_pval(r_corr.(prm_nm), pval.(prm_nm)); % add info about p.value and r on the graph
end