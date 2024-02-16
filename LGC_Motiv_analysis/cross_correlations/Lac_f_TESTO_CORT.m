function[r_corr, pval] = Lac_f_TESTO_CORT()
% [r_corr, pval] = Lac_f_TESTO_CORT()
% Lac_f_CORT will look at the correlation between plasma,
% dmPFC/dACC and anterior insula measures of lactate with all salivary 
% measures of cortisol, testosterone and testosterone/cortisol.
%
% OUTPUTS
% r_corr: correlation coefficient between plasma, dmPFC/dACC and aIns
% lactate and subjective stress ratings
%
% pval: structure with corresponding p.value

%% subject selection
[study_nm, condition, ~, subject_id, NS] = sub_id;

%% load lactate
% load plasma
[Lac_struct] = load_plasma_Lac(subject_id);
plasma_Lac = Lac_struct.Lac;
% load brain lactate
[metabolites] = metabolite_load(subject_id);
dmPFC_Lac = metabolites.dmPFC.Lac;
aIns_Lac = metabolites.aIns.Lac;

%% load hormonal concentrations
n_timepoints = 4;
timepoint_names = {'preMRS','postMRS',...
    'prefMRI','postfMRI'};
% load cortisol
[CORT_data] = load_CORT(study_nm, subject_id);
cortisol = CORT_data.CORT;
% load testosterone
[TESTO_data] = load_TESTO(study_nm, subject_id);
testosterone = TESTO_data.TESTO;
% compute testosterone/cortisol ratio
testoCort_ratio = testosterone./cortisol;

%% perform each correlation (with and without outliers)
for iRawCorr = 1:2
    switch iRawCorr
        case 1
            raw_or_corr_nm = 'raw';
        case 2
            raw_or_corr_nm = 'filtered';
    end
    
    for iTP = 1:n_timepoints
        timepoint_nm = timepoint_names{iTP};
        
        %% plasma lactate
        % cortisol
        plasma_cort_Lac_nm = [timepoint_nm,'_cort_f_plasma_Lac'];
        [goodS.(raw_or_corr_nm).(plasma_cort_Lac_nm)] = filter_fn(raw_or_corr_nm, plasma_Lac, cortisol(iTP, :));
        [r_corr.(raw_or_corr_nm).(plasma_cort_Lac_nm),...
            pval.(raw_or_corr_nm).(plasma_cort_Lac_nm)] = corr(plasma_Lac(goodS.(raw_or_corr_nm).(plasma_cort_Lac_nm))',...
            cortisol(iTP, goodS.(raw_or_corr_nm).(plasma_cort_Lac_nm))');
        % testosterone
        plasma_testo_Lac_nm = [timepoint_nm,'_testo_f_plasma_Lac'];
        [goodS.(raw_or_corr_nm).(plasma_testo_Lac_nm)] = filter_fn(raw_or_corr_nm, plasma_Lac, testosterone(iTP, :));
        [r_corr.(raw_or_corr_nm).(plasma_testo_Lac_nm),...
            pval.(raw_or_corr_nm).(plasma_testo_Lac_nm)] = corr(plasma_Lac(goodS.(raw_or_corr_nm).(plasma_testo_Lac_nm))',...
            testosterone(iTP, goodS.(raw_or_corr_nm).(plasma_testo_Lac_nm))');
        % testosterone/cortisol
        plasma_testocort_Lac_nm = [timepoint_nm,'_testocort_f_plasma_Lac'];
        [goodS.(raw_or_corr_nm).(plasma_testocort_Lac_nm)] = filter_fn(raw_or_corr_nm, plasma_Lac, testoCort_ratio(iTP, :));
        [r_corr.(raw_or_corr_nm).(plasma_testocort_Lac_nm),...
            pval.(raw_or_corr_nm).(plasma_testocort_Lac_nm)] = corr(plasma_Lac(goodS.(raw_or_corr_nm).(plasma_testocort_Lac_nm))',...
            testoCort_ratio(iTP, goodS.(raw_or_corr_nm).(plasma_testocort_Lac_nm))');
        
        %% dmPFC/dACC Lactate
        % cortisol
        dmPFC_Lac_cort_nm = [timepoint_nm,'_cort_f_dmPFC_Lac'];
        [goodS.(raw_or_corr_nm).(dmPFC_Lac_cort_nm)] = filter_fn(raw_or_corr_nm, dmPFC_Lac, cortisol(iTP, :));
        [r_corr.(raw_or_corr_nm).(dmPFC_Lac_cort_nm),...
            pval.(raw_or_corr_nm).(dmPFC_Lac_cort_nm)] = corr(dmPFC_Lac(goodS.(raw_or_corr_nm).(dmPFC_Lac_cort_nm))',...
            cortisol(iTP, goodS.(raw_or_corr_nm).(dmPFC_Lac_cort_nm))');
        % testosterone
        dmPFC_Lac_testo_nm = [timepoint_nm,'_testo_f_dmPFC_Lac'];
        [goodS.(raw_or_corr_nm).(dmPFC_Lac_testo_nm)] = filter_fn(raw_or_corr_nm, dmPFC_Lac, testosterone(iTP, :));
        [r_corr.(raw_or_corr_nm).(dmPFC_Lac_testo_nm),...
            pval.(raw_or_corr_nm).(dmPFC_Lac_testo_nm)] = corr(dmPFC_Lac(goodS.(raw_or_corr_nm).(dmPFC_Lac_testo_nm))',...
            testosterone(iTP, goodS.(raw_or_corr_nm).(dmPFC_Lac_testo_nm))');
        % testosterone/cortisol
        dmPFC_Lac_testocort_nm = [timepoint_nm,'_testocort_f_dmPFC_Lac'];
        [goodS.(raw_or_corr_nm).(dmPFC_Lac_testocort_nm)] = filter_fn(raw_or_corr_nm, dmPFC_Lac, testoCort_ratio(iTP, :));
        [r_corr.(raw_or_corr_nm).(dmPFC_Lac_testocort_nm),...
            pval.(raw_or_corr_nm).(dmPFC_Lac_testocort_nm)] = corr(dmPFC_Lac(goodS.(raw_or_corr_nm).(dmPFC_Lac_testocort_nm))',...
            testoCort_ratio(iTP, goodS.(raw_or_corr_nm).(dmPFC_Lac_testocort_nm))');
        
        %% aIns Lactate
        % cortisol
        aIns_Lac_cort_nm = [timepoint_nm,'_cort_f_aIns_Lac'];
        [goodS.(raw_or_corr_nm).(aIns_Lac_cort_nm)] = filter_fn(raw_or_corr_nm, aIns_Lac, cortisol(iTP, :));
        [r_corr.(raw_or_corr_nm).(aIns_Lac_cort_nm),...
            pval.(raw_or_corr_nm).(aIns_Lac_cort_nm)] = corr(aIns_Lac(goodS.(raw_or_corr_nm).(aIns_Lac_cort_nm))',...
            cortisol(iTP, goodS.(raw_or_corr_nm).(aIns_Lac_cort_nm))');
        % testosterone
        aIns_Lac_testo_nm = [timepoint_nm,'_testo_f_aIns_Lac'];
        [goodS.(raw_or_corr_nm).(aIns_Lac_testo_nm)] = filter_fn(raw_or_corr_nm, aIns_Lac, testosterone(iTP, :));
        [r_corr.(raw_or_corr_nm).(aIns_Lac_testo_nm),...
            pval.(raw_or_corr_nm).(aIns_Lac_testo_nm)] = corr(aIns_Lac(goodS.(raw_or_corr_nm).(aIns_Lac_testo_nm))',...
            testosterone(iTP, goodS.(raw_or_corr_nm).(aIns_Lac_testo_nm))');
        % testosterone/cortisol
        aIns_Lac_testocort_nm = [timepoint_nm,'_testocort_f_aIns_Lac'];
        [goodS.(raw_or_corr_nm).(aIns_Lac_testocort_nm)] = filter_fn(raw_or_corr_nm, aIns_Lac, testoCort_ratio(iTP, :));
        [r_corr.(raw_or_corr_nm).(aIns_Lac_testocort_nm),...
            pval.(raw_or_corr_nm).(aIns_Lac_testocort_nm)] = corr(aIns_Lac(goodS.(raw_or_corr_nm).(aIns_Lac_testocort_nm))',...
            testoCort_ratio(iTP, goodS.(raw_or_corr_nm).(aIns_Lac_testocort_nm))');
    end % loop over timepoints
    
    %% extract significant correlation if any
    correl_fields = fieldnames(pval.(raw_or_corr_nm));
    for iC = 1:length(correl_fields)
        correl_nm = correl_fields{iC};
        if pval.(raw_or_corr_nm).(correl_nm) < 0.05
            pval.signif.(raw_or_corr_nm).(correl_nm) = pval.(raw_or_corr_nm).(correl_nm);
        end % p.value filter
    end % loop over correlation
end % raw/filtered

end % function