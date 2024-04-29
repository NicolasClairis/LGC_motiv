function[r_corr, pval] = Lac_f_subjective_stress_rtg()
% [r_corr, pval] = Lac_f_subjective_stress_rtg()
% Lac_f_subjective_stress_rtg will look at the correlation between plasma,
% dmPFC/dACC and anterior insula measures of lactate with all subjective
% stress ratings
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

%% load subjective stress ratings
n_stress_rtgs = 4;
stress_rtg = NaN(n_stress_rtgs, NS);
[~,...
    stress_rtg(1,:), stress_rtg(2,:),...
    stress_rtg(3,:), stress_rtg(4,:)] = extract_subjective_stress_ratings(study_nm, subject_id, NS);
stress_rtg_names = {'preMRS_stress','postMRS_stress',...
    'prefMRI_stress','postfMRI_stress'};

%% perform each correlation (with and without outliers)
for iRawCorr = 1:2
    switch iRawCorr
        case 1
            raw_or_corr_nm = 'raw';
        case 2
            raw_or_corr_nm = 'filtered';
    end
    
    for iS_rtg = 1:n_stress_rtgs
        stress_rtg_nm = stress_rtg_names{iS_rtg};
        
        % plasma lactate
        plasma_Lac_nm = [stress_rtg_nm,'_f_plasma_Lac'];
        [goodS.(raw_or_corr_nm).(plasma_Lac_nm)] = filter_fn(raw_or_corr_nm, plasma_Lac, stress_rtg(iS_rtg, :));
        [r_corr.(raw_or_corr_nm).(plasma_Lac_nm),...
            pval.(raw_or_corr_nm).(plasma_Lac_nm)] = corr(plasma_Lac(goodS.(raw_or_corr_nm).(plasma_Lac_nm))', stress_rtg(iS_rtg, goodS.(raw_or_corr_nm).(plasma_Lac_nm))');
        
        % dmPFC/dACC Lactate
        dmPFC_Lac_nm = [stress_rtg_nm,'_f_dmPFC_Lac'];
        [goodS.(raw_or_corr_nm).(dmPFC_Lac_nm)] = filter_fn(raw_or_corr_nm, dmPFC_Lac, stress_rtg(iS_rtg, :));
        [r_corr.(raw_or_corr_nm).(dmPFC_Lac_nm),...
            pval.(raw_or_corr_nm).(dmPFC_Lac_nm)] = corr(dmPFC_Lac(goodS.(raw_or_corr_nm).(dmPFC_Lac_nm))', stress_rtg(iS_rtg, goodS.(raw_or_corr_nm).(dmPFC_Lac_nm))');
                
        % aIns Lactate
        aIns_Lac_nm = [stress_rtg_nm,'_f_aIns_Lac'];
        [goodS.(raw_or_corr_nm).(aIns_Lac_nm)] = filter_fn(raw_or_corr_nm, aIns_Lac, stress_rtg(iS_rtg, :));
        [r_corr.(raw_or_corr_nm).(aIns_Lac_nm),...
            pval.(raw_or_corr_nm).(aIns_Lac_nm)] = corr(aIns_Lac(goodS.(raw_or_corr_nm).(aIns_Lac_nm))', stress_rtg(iS_rtg, goodS.(raw_or_corr_nm).(aIns_Lac_nm))');
        
    end % loop over stress rating
    
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