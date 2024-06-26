%% script to correlate brain metabolites with behavioral parameters
%
% designed by N.Clairis - 2022

%% display correlations
fig_disp = 1; % 1=yes, 0=no

%% define all subjects
study_nm = 'study1';
[male_CIDS, female_CIDS, male_NS, female_NS, condition] = subject_selection_per_sex();

%% define metabolite and ROI you want to focus on
[metabolite_males_v0, MRS_ROI_nm, mb_nm] = metabolite_extraction(study_nm, male_CIDS);
[metabolite_females_v0, MRS_ROI_nm, mb_nm] = metabolite_extraction(study_nm, female_CIDS);

%% extract behavioral parameters
prm_males = prm_extraction(study_nm, male_CIDS);
prm_females = prm_extraction(study_nm, female_CIDS);
parameters = fieldnames(prm_males);

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
                    metabolite_males = metabolite_males_v0;
                    metabolite_females = metabolite_females_v0;
                    goodSubs.males.(uncCorr_nm).(prm_nm) = ~isnan(metabolite_males_v0.*prm_males_tmp);
                    goodSubs.females.(uncCorr_nm).(prm_nm) = ~isnan(metabolite_females_v0.*prm_females_tmp);
                case 2 % remove any outlier in any measure
                    [~,~,metabolite_males] = rmv_outliers_3sd(metabolite_males_v0);
                    [~,~,metabolite_females] = rmv_outliers_3sd(metabolite_females_v0);
                    [~,~,prm_males_tmp] = rmv_outliers_3sd(prm_males_tmp);
                    [~,~,prm_females_tmp] = rmv_outliers_3sd(prm_females_tmp);
                    goodSubs.males.(uncCorr_nm).(prm_nm) = ~isnan(metabolite_males.*prm_males_tmp);
                    goodSubs.females.(uncCorr_nm).(prm_nm) = ~isnan(metabolite_females.*prm_females_tmp);
            end
            
            NS_goodS.males.(uncCorr_nm).(prm_nm) = sum(goodSubs.males.(uncCorr_nm).(prm_nm));
            NS_goodS.females.(uncCorr_nm).(prm_nm) = sum(goodSubs.females.(uncCorr_nm).(prm_nm));
            [r_corr.(uncCorr_nm).(prm_nm).males, pval_corr.(uncCorr_nm).(prm_nm).males] = corr(prm_males_tmp(goodSubs.males.(uncCorr_nm).(prm_nm))',...
                metabolite_males_v0(goodSubs.males.(uncCorr_nm).(prm_nm))');
            [r_corr.(uncCorr_nm).(prm_nm).females, pval_corr.(uncCorr_nm).(prm_nm).females] = corr(prm_females_tmp(goodSubs.females.(uncCorr_nm).(prm_nm))',...
                metabolite_females_v0(goodSubs.females.(uncCorr_nm).(prm_nm))');
            [~, betas.(uncCorr_nm).(prm_nm).males,...
                pval.(uncCorr_nm).(prm_nm).males, ~,...
                mb_sorted.(uncCorr_nm).(prm_nm).males,...
                prm_fit_mbSorted.(uncCorr_nm).(prm_nm).males] = glm_package(metabolite_males_v0(goodSubs.males.(uncCorr_nm).(prm_nm))',...
                prm_males_tmp(goodSubs.males.(uncCorr_nm).(prm_nm))', 'normal', 'on');
            [~, betas.(uncCorr_nm).(prm_nm).females,...
                pval.(uncCorr_nm).(prm_nm).females, ~,...
                mb_sorted.(uncCorr_nm).(prm_nm).females,...
                prm_fit_mbSorted.(uncCorr_nm).(prm_nm).females] = glm_package(metabolite_females_v0(goodSubs.females.(uncCorr_nm).(prm_nm))',...
                prm_females_tmp(goodSubs.females.(uncCorr_nm).(prm_nm))', 'normal', 'on');
        end % filter to only consider parameters and no other subfield
    end % parameters loop
    
    %% display results
    if fig_disp == 1
        pSize = 50;
        % load color for each sex
        [col_m1, col_f1, col_m2, col_f2] = col_per_sex;

        % loop through parameters
        for iPrm = 1:n_prm
            prm_nm = parameters{iPrm};
            if ~strcmp(prm_nm,'CID') % no sense in this case
                fig;
                scat_males_hdl = scatter(metabolite_males_v0(goodSubs.males.(uncCorr_nm).(prm_nm))',...
                    prm_males.(prm_nm)(goodSubs.males.(uncCorr_nm).(prm_nm))');
                scat_females_hdl = scatter(metabolite_females_v0(goodSubs.females.(uncCorr_nm).(prm_nm))',...
                    prm_females.(prm_nm)(goodSubs.females.(uncCorr_nm).(prm_nm))');
                fit_males_hdl = plot(mb_sorted.(uncCorr_nm).(prm_nm).males,...
                    prm_fit_mbSorted.(uncCorr_nm).(prm_nm).males);
                fit_females_hdl = plot(mb_sorted.(uncCorr_nm).(prm_nm).females,...
                    prm_fit_mbSorted.(uncCorr_nm).(prm_nm).females);
                scat_hdl_upgrade(scat_males_hdl);
                scat_hdl_upgrade(scat_females_hdl);
                fit_hdl_upgrade(fit_males_hdl);
                fit_hdl_upgrade(fit_females_hdl);
                scat_males_hdl.MarkerFaceColor = col_m1;
                scat_females_hdl.MarkerFaceColor = col_f1;
                fit_males_hdl.Color = col_m2;
                fit_females_hdl.Color = col_f2;
                xlabel([MRS_ROI_nm,' - ',mb_nm,' (mM)']);
                ylabel(prm_nm);
                place_r_and_pval_2(r_corr.(uncCorr_nm).(prm_nm).males,...
                    pval_corr.(uncCorr_nm).(prm_nm).males,...
                    col_m2,...
                    r_corr.(uncCorr_nm).(prm_nm).females,...
                    pval_corr.(uncCorr_nm).(prm_nm).females,...
                    col_f2);
                legend_size(pSize);
            end
        end % parameter loop
    end % figure display
end % raw/outlier cleaned