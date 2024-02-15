%% script to correlate brain metabolites with behavioral parameters
%
% designed by N.Clairis - 2022

%% define all subjects
study_nm = 'study1';
condition = subject_condition();
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% define metabolite and ROI you want to focus on
[metabolite_allSubs_v0, MRS_ROI_nm, mb_nm] = metabolite_extraction(study_nm, subject_id);

%% extract behavioral parameters
prm = prm_extraction(study_nm, subject_id);
parameters = fieldnames(prm);

%% try uncorrected (all subjects) and corrected for outliers (mean+/-3*SD) analysis
subsIncluded = {'allSubs','withoutOutliers'};

for iUncorrCorr = 1:length(subsIncluded)
    uncCorr_nm = subsIncluded{iUncorrCorr};
    
    %% perform GLM and correlation between behavioral parameters and spectroscopy
    n_prm = length(parameters);
    % remove NaN subjects
    for iPrm = 1:n_prm
        prm_nm = parameters{iPrm};
        prm_tmp = prm.(prm_nm);
        if ~strcmp(prm_nm,'CID') % no sense in this case
            % filter subjects depending on situation
            switch iUncorrCorr
                case 1 % all subjects included as long as no-NaN values
                    metabolite_allSubs = metabolite_allSubs_v0;
                    goodSubs.(uncCorr_nm).(prm_nm) = ~isnan(metabolite_allSubs_v0.*prm_tmp);
                case 2 % remove any outlier in any measure
                    [~,~,metabolite_allSubs] = rmv_outliers_3sd(metabolite_allSubs_v0);
                    [~,~,prm_tmp] = rmv_outliers_3sd(prm_tmp);
                    goodSubs.(uncCorr_nm).(prm_nm) = ~isnan(metabolite_allSubs.*prm_tmp);
            end
            
            NS_goodS.(uncCorr_nm).(prm_nm) = sum(goodSubs.(uncCorr_nm).(prm_nm));
            [r_corr.(uncCorr_nm).(prm_nm),pval_corr.(uncCorr_nm).(prm_nm)] = corr(prm_tmp(goodSubs.(uncCorr_nm).(prm_nm))',...
                metabolite_allSubs_v0(goodSubs.(uncCorr_nm).(prm_nm))');
            [~, betas.(uncCorr_nm).(prm_nm),...
                pval.(uncCorr_nm).(prm_nm), ~,...
                mb_sorted.(uncCorr_nm).(prm_nm),...
                prm_fit_mbSorted.(uncCorr_nm).(prm_nm)] = glm_package(metabolite_allSubs_v0(goodSubs.(uncCorr_nm).(prm_nm))',...
                prm_tmp(goodSubs.(uncCorr_nm).(prm_nm))', 'normal', 'on');
        end
    end % parameters loop
    
    %% display results
    pSize = 50;
    for iPrm = 1:n_prm
        prm_nm = parameters{iPrm};
        if ~strcmp(prm_nm,'CID') % no sense in this case
            fig;
            scat_hdl = scatter(metabolite_allSubs_v0(goodSubs.(uncCorr_nm).(prm_nm))',...
                prm.(prm_nm)(goodSubs.(uncCorr_nm).(prm_nm))');
            fit_hdl = plot(mb_sorted.(uncCorr_nm).(prm_nm),...
                prm_fit_mbSorted.(uncCorr_nm).(prm_nm));
            scat_hdl_upgrade(scat_hdl);
            fit_hdl_upgrade(fit_hdl);
            xlabel([MRS_ROI_nm,' - ',mb_nm,' (mM)']);
            ylabel(prm_nm);
            place_r_and_pval(r_corr.(uncCorr_nm).(prm_nm), pval_corr.(uncCorr_nm).(prm_nm));
            legend_size(pSize);
        end
    end % parameter loop
    
end % raw/outlier cleaned