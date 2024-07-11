%% script to correlate brain metabolites with behavioral parameters
%
% designed by N.Clairis - 2022

%% display correlations
fig_disp = 1; % 1=yes, 0=no

%% define all subjects
study_nm = 'study1';
condition = subject_condition();
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% define metabolite and ROI you want to focus on
[metabolite_allSubs_v0, MRS_ROI_nm, mb_nm] = metabolite_extraction(study_nm, subject_id);

%% extract behavioral parameters
prm = prm_extraction(study_nm, subject_id);
parameters = fieldnames(prm);

%% load age and sex of the participants
[questionnaires, categ_quests, n_quest_categ, subject_id] = extract_questionnaires(study_nm, subject_id, NS);
age = questionnaires.general.age;
sex = questionnaires.general.sex;

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
            goodS_tmp = goodSubs.(uncCorr_nm).(prm_nm);
            
            NS_goodS.(uncCorr_nm).(prm_nm) = sum(goodS_tmp);
            % perform GLM including age and sex co-regressors
            x_vars = [age(goodS_tmp)', sex(goodS_tmp)', metabolite_allSubs_v0(goodS_tmp)'];
            [~, betas_tmp,...
                pval_tmp, ~,...
                mb_sorted.(uncCorr_nm).(prm_nm),...
                prm_fit_mbSorted.(uncCorr_nm).(prm_nm)] = glm_package(x_vars,...
                prm_tmp(goodS_tmp)', 'normal', 'on');
            betas.(uncCorr_nm).(prm_nm).cstt = betas_tmp(1);
            betas.(uncCorr_nm).(prm_nm).age = betas_tmp(2);
            betas.(uncCorr_nm).(prm_nm).sex = betas_tmp(3);
            betas.(uncCorr_nm).(prm_nm).(prm_nm) = betas_tmp(4);
            pval.(uncCorr_nm).(prm_nm).cstt         = pval_tmp(1);
            pval.(uncCorr_nm).(prm_nm).age          = pval_tmp(2);
            pval.(uncCorr_nm).(prm_nm).sex          = pval_tmp(3);
            pval.(uncCorr_nm).(prm_nm).(prm_nm)     = pval_tmp(4);
            
            % orthogonalize variable of interest to age and sex and perform
            % correlation on the orthogonalized variable
            % 1) orthogonalization
            prm_tmp_orthogonalized_age_sex = prm_tmp(goodS_tmp)' -...
                betas.(uncCorr_nm).(prm_nm).age.*age(goodS_tmp)' -...
                betas.(uncCorr_nm).(prm_nm).sex.*sex(goodS_tmp)';
            % 2) correlation with orthogonalized variables
            [r_corr.(uncCorr_nm).(prm_nm),pval_corr.(uncCorr_nm).(prm_nm)] = corr(prm_tmp_orthogonalized_age_sex,...
                metabolite_allSubs_v0(goodS_tmp)');
            
        end
    end % parameters loop
    
    %% display results
    if fig_disp == 1
        pSize = 50;
        for iPrm = 1:n_prm
            prm_nm = parameters{iPrm};
            if ~strcmp(prm_nm,'CID') % no sense in this case
                goodS_tmp = goodSubs.(uncCorr_nm).(prm_nm);
                fig;
                scat_hdl = scatter(metabolite_allSubs_v0(goodS_tmp)',...
                    prm.(prm_nm)(goodS_tmp)');
%                 fit_hdl = plot(mb_sorted.(uncCorr_nm).(prm_nm),...
%                     prm_fit_mbSorted.(uncCorr_nm).(prm_nm));
                scat_hdl_upgrade(scat_hdl);
%                 fit_hdl_upgrade(fit_hdl);
                xlabel([MRS_ROI_nm,' - ',mb_nm,' (mM)']);
                ylabel(prm_nm);
                place_r_and_pval(r_corr.(uncCorr_nm).(prm_nm),...
                    pval_corr.(uncCorr_nm).(prm_nm));
                legend_size(pSize);
            end
        end % parameter loop
    end % figure display
end % raw/outlier cleaned