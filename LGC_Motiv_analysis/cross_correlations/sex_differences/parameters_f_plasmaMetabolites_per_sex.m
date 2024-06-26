%% script to correlate plasma metabolites with behavioral parameters
%
% designed by N.Clairis - 2024

%% display correlations
fig_disp = 1; % 1=yes, 0=no

%% define all subjects
study_nm = 'study1';
condition = subject_condition();
[male_CIDS, female_CIDS, male_NS, female_NS, condition] = subject_selection_per_sex();

%% define metabolite and ROI you want to focus on
[plasmaM_males, mb_names, n_mb] = load_plasma_metabolites(male_CIDS);
[plasmaM_females, mb_names, n_mb] = load_plasma_metabolites(female_CIDS);
% metabolite selection
plasma_mb_males = plasmaM_males.Gln./1000; % convert from μM to mM
plasma_mb_females = plasmaM_females.Gln./1000; % convert from μM to mM
mb_nm = 'glutamine';

metabolite_males = plasma_mb_males;
metabolite_females = plasma_mb_females;

%% extract behavioral parameters
prm_males = prm_extraction(study_nm, male_CIDS);
prm_females = prm_extraction(study_nm, female_CIDS);
parameters = fieldnames(prm_males);

%% extract color for each sex
[col_m1, col_f1, col_m2, col_f2] = col_per_sex;

%% try uncorrected (all subjects) and corrected for outliers (mean+/-3*SD) analysis
subsIncluded = {'allSubs','withoutOutliers'};

for iUncorrCorr = 1:length(subsIncluded)
    uncCorr_nm = subsIncluded{iUncorrCorr};
    
    %% perform GLM and correlation between behavioral parameters and spectroscopy
    n_prm = length(parameters);
    % remove NaN subjects
    for iPrm = 1:n_prm
        prm_nm = parameters{iPrm};
        prm_males_tmp = prm_males.(prm_nm);
        prm_females_tmp = prm_females.(prm_nm);
        if ~strcmp(prm_nm,'CID') % no sense in this case
            % filter subjects depending on situation
            switch iUncorrCorr
                case 1 % all subjects included as long as no-NaN values
                    goodSubs.males.(uncCorr_nm).(prm_nm) = ~isnan(metabolite_males.*prm_males_tmp);
                    goodSubs.females.(uncCorr_nm).(prm_nm) = ~isnan(metabolite_females.*prm_females_tmp);
                case 2 % remove any outlier in any measure
                    [~,~,metabolite_males] = rmv_outliers_3sd(metabolite_males);
                    [~,~,prm_males_tmp] = rmv_outliers_3sd(prm_males_tmp);
                    [~,~,prm_females_tmp] = rmv_outliers_3sd(prm_females_tmp);
                    goodSubs.males.(uncCorr_nm).(prm_nm) = ~isnan(metabolite_males.*prm_males_tmp);
                    goodSubs.females.(uncCorr_nm).(prm_nm) = ~isnan(metabolite_females.*prm_females_tmp);
            end
            
            % extract nb of correct subjects
            NS_goodS.males.(uncCorr_nm).(prm_nm) = sum(goodSubs.males.(uncCorr_nm).(prm_nm));
            NS_goodS.females.(uncCorr_nm).(prm_nm) = sum(goodSubs.females.(uncCorr_nm).(prm_nm));
            % perform correlation
            [r_corr.males.(uncCorr_nm).(prm_nm), pval_corr.males.(uncCorr_nm).(prm_nm)] = corr(prm_males_tmp(goodSubs.males.(uncCorr_nm).(prm_nm))',...
                metabolite_males(goodSubs.males.(uncCorr_nm).(prm_nm))');
            [r_corr.females.(uncCorr_nm).(prm_nm), pval_corr.females.(uncCorr_nm).(prm_nm)] = corr(prm_females_tmp(goodSubs.females.(uncCorr_nm).(prm_nm))',...
                metabolite_females(goodSubs.females.(uncCorr_nm).(prm_nm))');
            % GLM to extract fit
            [~, betas.males.(uncCorr_nm).(prm_nm),...
                pval.males.(uncCorr_nm).(prm_nm), ~,...
                mb_sorted.males.(uncCorr_nm).(prm_nm),...
                prm_fit_mbSorted.males.(uncCorr_nm).(prm_nm)] = glm_package(metabolite_males(goodSubs.males.(uncCorr_nm).(prm_nm))',...
                prm_males_tmp(goodSubs.males.(uncCorr_nm).(prm_nm))', 'normal', 'on');
            [~, betas.females.(uncCorr_nm).(prm_nm),...
                pval.females.(uncCorr_nm).(prm_nm), ~,...
                mb_sorted.females.(uncCorr_nm).(prm_nm),...
                prm_fit_mbSorted.females.(uncCorr_nm).(prm_nm)] = glm_package(metabolite_females(goodSubs.females.(uncCorr_nm).(prm_nm))',...
                prm_females_tmp(goodSubs.females.(uncCorr_nm).(prm_nm))', 'normal', 'on');
        end
    end % parameters loop
    
    %% display results
    if fig_disp == 1
        pSize = 50;
        for iPrm = 1:n_prm
            prm_nm = parameters{iPrm};
            if ~strcmp(prm_nm,'CID') % no sense in this case
                fig;
                % correlation in males
                scat_male_hdl = scatter(metabolite_males(goodSubs.males.(uncCorr_nm).(prm_nm))',...
                    prm_males.(prm_nm)(goodSubs.males.(uncCorr_nm).(prm_nm))');
                fit_male_hdl = plot(mb_sorted.males.(uncCorr_nm).(prm_nm),...
                    prm_fit_mbSorted.males.(uncCorr_nm).(prm_nm));
                % correlation in females
                scat_female_hdl = scatter(metabolite_females(goodSubs.females.(uncCorr_nm).(prm_nm))',...
                    prm_females.(prm_nm)(goodSubs.females.(uncCorr_nm).(prm_nm))');
                fit_female_hdl = plot(mb_sorted.females.(uncCorr_nm).(prm_nm),...
                    prm_fit_mbSorted.females.(uncCorr_nm).(prm_nm));

                scat_hdl_upgrade(scat_male_hdl);
                scat_hdl_upgrade(scat_female_hdl);
                fit_hdl_upgrade(fit_male_hdl);
                fit_hdl_upgrade(fit_female_hdl);
                % change color for each sex
                scat_male_hdl.MarkerFaceColor = col_m1;
                scat_female_hdl.MarkerFaceColor = col_f1;
                fit_male_hdl.Color = col_m2;
                fit_female_hdl.Color = col_f2;
                xlabel(['plasma ',mb_nm,' (mM)']);
                ylabel(prm_nm);
                % add correlation coeffs and p.values
                [txt1_hdl, txt2_hdl, txtSize] = place_r_and_pval_2(r_corr.males.(uncCorr_nm).(prm_nm),...
                    pval_corr.males.(uncCorr_nm).(prm_nm),...
                    fit_male_hdl.Color,...
                    r_corr.females.(uncCorr_nm).(prm_nm),...
                    pval_corr.females.(uncCorr_nm).(prm_nm),...
                    fit_female_hdl.Color);
                legend_size(pSize);
            end
        end % parameter loop
    end % figure display
end % raw/outlier cleaned