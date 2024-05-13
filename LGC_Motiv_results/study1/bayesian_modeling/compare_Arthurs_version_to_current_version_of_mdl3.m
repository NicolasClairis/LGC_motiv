% script to compare former parameters extracted by Arthur with current version
%
%
% Script changes comprise:
%
% - a higher maximal number of VBA iterations (Arthur fixed it at 10, I
% brought back default value of 32 iterations)
%
% - different runs excluded from saturation => I use a more stringent
% threshold (94%) while he only used 100% as a threshold.
%
% - he used current efficiency while removing first trial of each run since
% he was planning to remove previous efficiency

%% working dir
root_dir = fullfile('C:','Users','clairis','Desktop');
working_dir = fullfile(root_dir,'GitHub','LGC_motiv',...
    'LGC_Motiv_results','study1','bayesian_modeling');
%% load Arthur data
bhv_prm_v1 = getfield(getfield(load([working_dir, filesep,'archives',filesep,'behavioral_prm.mat']),'bayesian_mdl'),'mdl_3');
subject_id_v1 = bhv_prm_v1.subject_id;
NS_v1 = length(subject_id_v1);
% replace kBiasM by kBias to make data fit between v1 and v2
bhv_prm_v1.kBias = bhv_prm_v1.kBiasM;
bhv_prm_v1 = rmfield(bhv_prm_v1, 'kBiasM');
bhv_prm_v1 = rmfield(bhv_prm_v1, 'subject_id');

%% load my version
bhv_prm_v2_struct = load('bayesian_model_3_results.mat');
bhv_prm_v2 = bhv_prm_v2_struct.prm;
subject_id_v2 = bhv_prm_v2_struct.subject_id;
NS_v2 = length(subject_id_v2);

%% compare data
[NS,idx] = max([NS_v1, NS_v2]);
switch idx
    case 1
        subject_id = subject_id_v1;
    case 2
        subject_id = subject_id_v2;
end
prm_names = fieldnames(bhv_prm_v2);
nPrm = length(prm_names);
fig;
for iP = 1:nPrm
    prm_nm = prm_names{iP};
    [prm_v1.(prm_nm), prm_v2.(prm_nm)] = deal(NaN(1, NS));
    
    for iS = 1:NS
        sub_nm = subject_id{iS};
        % identify subject for each vector of parameters
        sub_idx_v1 = strcmp(sub_nm,subject_id_v1);
        sub_idx_v2 = strcmp(sub_nm,subject_id_v2);
        
        % extract corresponding parameter
        if sum(sub_idx_v1) > 0
            prm_v1.(prm_nm)(iS) = bhv_prm_v1.(prm_nm)(sub_idx_v1);
        end
        if sum(sub_idx_v2) > 0
            prm_v2.(prm_nm)(iS) = bhv_prm_v2.(prm_nm)(sub_idx_v2);
        end
    end % subject loop
    
    % test correlation
    ok_subs = ~isnan(prm_v1.(prm_nm).*prm_v2.(prm_nm));
    [rho.(prm_nm), pval.(prm_nm)] = corr(prm_v1.(prm_nm)(ok_subs)', prm_v2.(prm_nm)(ok_subs)');
    % display data
    subplot(3,3,iP);
    scat_hdl = scatter(prm_v1.(prm_nm), prm_v2.(prm_nm));
    scat_hdl_upgrade(scat_hdl);
    xlabel([prm_nm,' - Arthur version']);
    ylabel([prm_nm,' - New version']);
    place_r_and_pval(rho.(prm_nm), pval.(prm_nm));
end % parameter loop


% summary: most parameters are still relatively close (r > 0.8 for all)
% but the changes still induce some non-negligible changes that may impact
% the results