function[correl] = questionnaire_IPAQ_f_Lac(disp_fig)

%% display figure by default
if ~exist('disp_fig','var') || isempty(disp_fig) || ~ismember(disp_fig,[0,1])
    disp_fig = 1;
end

%% subject selection
[study_nm, condition, subject_id, NS] = subject_selection;

%% load IPAQ values
[IPAQ_score, IPAQ_struct, IPAQinactivity] = IPAQ_rescoring(study_nm, subject_id, NS);
IPAQ.total = IPAQ_score';
IPAQ.high = IPAQ_struct.MET_high';
IPAQ.mod = IPAQ_struct.MET_mod';
IPAQ.low = IPAQ_struct.MET_low';
IPAQ.sit = IPAQinactivity';
IPAQ_fields = fieldnames(IPAQ);
n_IPAQ = length(IPAQ_fields);

%% load lactate values
% plasma
[plasma_Lac_struct] = load_plasma_Lac(subject_id);
Lac.plasma = plasma_Lac_struct.Lac'./1000;
% brain
[brainMetabolites] = metabolite_load(subject_id);
Lac.dmPFC = brainMetabolites.dmPFC.Lac';
Lac.aIns = brainMetabolites.aIns.Lac';

Lac_fields = fieldnames(Lac);
n_Lac = length(Lac_fields);

%% perform correlation tests between IPAQ measures and lactate with linear and U-shape tests
for iIPAQ = 1:n_IPAQ
    IPAQ_field_nm = IPAQ_fields{iIPAQ};
    for iLac = 1:n_Lac
        lac_nm = Lac_fields{iLac};
        corr_nm = [IPAQ_field_nm,'_f_',lac_nm];
        %% linear tests
        % raw values
        goodS.(corr_nm).raw = ~isnan(IPAQ.(IPAQ_field_nm).*Lac.(lac_nm));
        NS_goodS.(corr_nm).noOutliers = sum(goodS.(corr_nm).raw );
        [correl.(corr_nm).raw.linear.r_corr, correl.(corr_nm).raw.linear.pval] = corr(Lac.(lac_nm)(goodS.(corr_nm).raw), IPAQ.(IPAQ_field_nm)(goodS.(corr_nm).raw));
        [~, ~, ~, ~, Lac_sorted.(corr_nm).raw.linear, IPAQ_sorted.(corr_nm).raw.linear] = glm_package(Lac.(lac_nm), IPAQ.(IPAQ_field_nm),'normal'); % extract corresponding fit
        
        % remove outliers
        [~,~,~,idx_goodS.(IPAQ_field_nm)] = rmv_outliers_3sd(IPAQ.(IPAQ_field_nm));
        [~,~,~,idx_goodS.(lac_nm)] = rmv_outliers_3sd(Lac.(lac_nm));
        goodS.(corr_nm).noOutliers = (idx_goodS.(IPAQ_field_nm)'.*idx_goodS.(lac_nm)') == 1;
        NS_goodS.(corr_nm).noOutliers = sum(goodS.(corr_nm).noOutliers);
        [correl.(corr_nm).noOutliers.linear.r_corr, correl.(corr_nm).noOutliers.linear.pval] = corr(Lac.(lac_nm)(goodS.(corr_nm).noOutliers), IPAQ.(IPAQ_field_nm)(goodS.(corr_nm).noOutliers));
        [~, ~, ~, ~, Lac_sorted.(corr_nm).noOutliers.linear, IPAQ_sorted.(corr_nm).noOutliers.linear] = glm_package(Lac.(lac_nm)(goodS.(corr_nm).noOutliers), IPAQ.(IPAQ_field_nm)(goodS.(corr_nm).noOutliers),'normal'); % extract corresponding fit
        
        %% U-shape (mean-centered + squared)
        Ushaped_lac = (Lac.(lac_nm) - mean(Lac.(lac_nm),'omitnan')).^2;
        % raw values
        [correl.(corr_nm).raw.Ushape.r_corr, correl.(corr_nm).raw.Ushape.pval] = corr(Ushaped_lac(goodS.(corr_nm).raw), IPAQ.(IPAQ_field_nm)(goodS.(corr_nm).raw));
        [~, ~, ~, IPAQ_fit.(corr_nm).raw.Ushape] = glm_package(Ushaped_lac, IPAQ.(IPAQ_field_nm),'normal'); % extract corresponding fit
        [Lac_sorted.(corr_nm).raw.Ushape, sort_idx] = sort(Lac.(lac_nm));
        IPAQ_sorted.(corr_nm).raw.Ushape = IPAQ_fit.(corr_nm).raw.Ushape(sort_idx);
        % remove outliers (identified before transformation)
        [correl.(corr_nm).noOutliers.Ushape.r_corr, correl.(corr_nm).noOutliers.Ushape.pval] = corr(Ushaped_lac(goodS.(corr_nm).noOutliers), IPAQ.(IPAQ_field_nm)(goodS.(corr_nm).noOutliers));
        [~, ~, ~, IPAQ_fit.(corr_nm).noOutliers.Ushape] = glm_package(Ushaped_lac(goodS.(corr_nm).noOutliers), IPAQ.(IPAQ_field_nm)(goodS.(corr_nm).noOutliers),'normal'); % extract corresponding fit
        [Lac_sorted.(corr_nm).noOutliers.Ushape, sort_idx2] = sort(Lac.(lac_nm)(goodS.(corr_nm).noOutliers));
        IPAQ_sorted.(corr_nm).noOutliers.Ushape = IPAQ_fit.(corr_nm).noOutliers.Ushape(sort_idx2);
    end % plasma/dmPFC/aIns lactate
end % IPAQ measures

%% display results
if disp_fig == 1
    for iIPAQ = 1:n_IPAQ
        IPAQ_field_nm = IPAQ_fields{iIPAQ};
        
        %% linear test
        fig;
        
        for iLac = 1:n_Lac
            lac_nm = Lac_fields{iLac};
            corr_nm = [IPAQ_field_nm,'_f_',lac_nm];
            
            % raw tests
            subplot(2,3,iLac); hold on;
            scat_hdl = scatter(Lac.(lac_nm), IPAQ.(IPAQ_field_nm));
            scat_hdl_upgrade(scat_hdl);
            fit_hdl = plot(Lac_sorted.(corr_nm).raw.linear, IPAQ_sorted.(corr_nm).raw.linear);
            fit_hdl_upgrade(fit_hdl);
            place_r_and_pval(correl.(corr_nm).raw.linear.r_corr,...
                correl.(corr_nm).raw.linear.pval);
            xlabel([lac_nm,' lactate (mM)']);
            ylabel(['IPAQ ',IPAQ_field_nm]);
            
            % test with no outliers
            subplot(2,3,iLac + 3); hold on;
            scat_hdl = scatter(Lac.(lac_nm)(goodS.(corr_nm).noOutliers), IPAQ.(IPAQ_field_nm)(goodS.(corr_nm).noOutliers));
            scat_hdl_upgrade(scat_hdl);
            fit_hdl = plot(Lac_sorted.(corr_nm).noOutliers.linear, IPAQ_sorted.(corr_nm).noOutliers.linear);
            fit_hdl_upgrade(fit_hdl);
            place_r_and_pval(correl.(corr_nm).noOutliers.linear.r_corr,...
                correl.(corr_nm).noOutliers.linear.pval);
            xlabel([lac_nm,' lactate (mM)']);
            ylabel(['IPAQ ',IPAQ_field_nm]);
        end % lactate
        
        %% U-shaped test
        fig;
        
        for iLac = 1:n_Lac
            lac_nm = Lac_fields{iLac};
            corr_nm = [IPAQ_field_nm,'_f_',lac_nm];
            
            % raw tests
            subplot(2,3,iLac); hold on;
            scat_hdl = scatter(Lac.(lac_nm), IPAQ.(IPAQ_field_nm));
            scat_hdl_upgrade(scat_hdl);
            fit_hdl = plot(Lac_sorted.(corr_nm).raw.Ushape, IPAQ_sorted.(corr_nm).raw.Ushape);
            fit_hdl_upgrade(fit_hdl);
            place_r_and_pval(correl.(corr_nm).raw.Ushape.r_corr,...
                correl.(corr_nm).raw.Ushape.pval);
            xlabel([lac_nm,' lactate (mM)']);
            ylabel(['IPAQ ',IPAQ_field_nm]);
            
            % test with no outliers
            subplot(2,3,iLac + 3); hold on;
            scat_hdl = scatter(Lac.(lac_nm)(goodS.(corr_nm).noOutliers), IPAQ.(IPAQ_field_nm)(goodS.(corr_nm).noOutliers));
            scat_hdl_upgrade(scat_hdl);
            fit_hdl = plot(Lac_sorted.(corr_nm).noOutliers.Ushape, IPAQ_sorted.(corr_nm).noOutliers.Ushape);
            fit_hdl_upgrade(fit_hdl);
            place_r_and_pval(correl.(corr_nm).noOutliers.Ushape.r_corr,...
                correl.(corr_nm).noOutliers.Ushape.pval);
            xlabel([lac_nm,' lactate (mM)']);
            ylabel(['IPAQ ',IPAQ_field_nm]);
        end % lactate
        
    end % IPAQ loop
end % display figures

end % function