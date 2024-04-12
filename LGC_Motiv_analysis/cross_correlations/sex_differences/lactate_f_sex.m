function[plasma_Lac, dmPFC_Lac, aIns_Lac] = lactate_f_sex(fig_disp, rmv_outliers_yn)
% [plasma_Lac, dmPFC_Lac, aIns_Lac] = lactate_f_sex(fig_disp, rmv_outliers_yn)
% lactate_f_sex compares plasma, dmPFC/dACC and aIns lactate levels in
% males vs females.
% INPUTS
% fig_disp: display figures? (1) yes (0) no
%
% rmv_outliers_yn: remove median +/- 3*SD outliers yes (1) or no (0)? Yes by
% default
%
% OUTPUTS
% plasma_Lac, dmPFC_Lac, aIns_Lac: structure with p.value (for unpaired
% t.test), and mean and SEM for each measure and each sex

%% define inputs by default
% display figure by default
if ~exist('fig_disp','var') || isempty(fig_disp) || ~ismember(fig_disp,[0,1])
    fig_disp = 1;
end
% remove outliers by default
if ~exist('rmv_outliers_yn','var') || isempty(rmv_outliers_yn) || ~ismember(rmv_outliers_yn,[0,1])
    rmv_outliers_yn = 1;
end

%% subject selection
study_nm = 'study1';
[male_CIDS, female_CIDS, male_NS, female_NS, condition] = subject_selection_per_sex;

%% extract lactate levels
% plasma
plasmaM_males = load_plasma_metabolites(male_CIDS);
plasmaM_females = load_plasma_metabolites(female_CIDS);
plasmaLac_males = plasmaM_males.Lac./1000;
plasmaLac_females = plasmaM_females.Lac./1000;

% load brain metabolites
metabolites_males = metabolite_load(male_CIDS);
metabolites_females = metabolite_load(female_CIDS);
dmPFC_Lac_males = metabolites_males.dmPFC.Lac;
dmPFC_Lac_females = metabolites_females.dmPFC.Lac;
aIns_Lac_males = metabolites_males.aIns.Lac;
aIns_Lac_females = metabolites_females.aIns.Lac;

%% remove outliers
if rmv_outliers_yn == 1
    % plasma
    [~,~,plasmaLac_males] = rmv_outliers_3sd(plasmaLac_males);
    [~,~,plasmaLac_females] = rmv_outliers_3sd(plasmaLac_females);
    % brain
    [~,~,dmPFC_Lac_males] = rmv_outliers_3sd(dmPFC_Lac_males);
    [~,~,dmPFC_Lac_females] = rmv_outliers_3sd(dmPFC_Lac_females);
    [~,~,aIns_Lac_males] = rmv_outliers_3sd(aIns_Lac_males);
    [~,~,aIns_Lac_females] = rmv_outliers_3sd(aIns_Lac_females);
end % outliers

%% extract mean, SEM per sex and perform the comparison
% plasma
[~,plasma_Lac.pval] = ttest2(plasmaLac_males, plasmaLac_females);
[plasma_Lac.m_males, plasma_Lac.sem_males,...
    plasma_Lac.sd_males, plasma_Lac.med_males] = mean_sem_sd(plasmaLac_males,2);
[plasma_Lac.m_females, plasma_Lac.sem_females,...
    plasma_Lac.sd_females, plasma_Lac.med_females] = mean_sem_sd(plasmaLac_females,2);
% dmPFC
[~,dmPFC_Lac.pval] = ttest2(dmPFC_Lac_males, dmPFC_Lac_females);
[dmPFC_Lac.m_males, dmPFC_Lac.sem_males,...
    dmPFC_Lac.sd_males, dmPFC_Lac.med_males] = mean_sem_sd(dmPFC_Lac_males,2);
[dmPFC_Lac.m_females, dmPFC_Lac.sem_females,...
    dmPFC_Lac.sd_females, dmPFC_Lac.med_females] = mean_sem_sd(dmPFC_Lac_females,2);
% aIns
[~,aIns_Lac.pval] = ttest2(aIns_Lac_males, aIns_Lac_females);
[aIns_Lac.m_males, aIns_Lac.sem_males,...
    aIns_Lac.sd_males, aIns_Lac.med_males] = mean_sem_sd(aIns_Lac_males,2);
[aIns_Lac.m_females, aIns_Lac.sem_females,...
    aIns_Lac.sd_females, aIns_Lac.med_females] = mean_sem_sd(aIns_Lac_females,2);

%% figure display
if fig_disp == 1
    % general figure parameters
    [pSize, lW, col, mSize] = general_fig_prm;
    female_col = col.red;
    male_col = col.blue_dark;
    
    % Lactate
    fig;
    
    % show male vs female data
    % dmPFC
    ok_dmPFC_males = ~isnan(dmPFC_Lac_males);
    Violin({dmPFC_Lac_males(ok_dmPFC_males)},1,...
        'ViolinColor',{male_col});
    ok_dmPFC_females = ~isnan(dmPFC_Lac_females);
    Violin({dmPFC_Lac_females(ok_dmPFC_females)},2,...
        'ViolinColor',{female_col});
    % aIns
    ok_aIns_males = ~isnan(aIns_Lac_males);
    Violin({aIns_Lac_males(ok_aIns_males)},4,...
        'ViolinColor',{male_col});
    ok_aIns_females = ~isnan(aIns_Lac_females);
    Violin({aIns_Lac_females(ok_aIns_females)},5,...
        'ViolinColor',{female_col});
    % plasma
    ok_plasma_males = ~isnan(plasmaLac_males);
    Violin({plasmaLac_males(ok_plasma_males)},7,...
        'ViolinColor',{male_col});
    ok_plasma_females = ~isnan(plasmaLac_females);
    Violin({plasmaLac_females(ok_plasma_females)},8,...
        'ViolinColor',{female_col});
    
    % add p.value indication if difference is significant
    % dmPFC
    [l_hdl, star_hdl] = add_pval_comparison(dmPFC_Lac_males,...
        dmPFC_Lac_females,...
        dmPFC_Lac.pval, 1, 2, 'NS');
    % aIns
    [l_hdl, star_hdl] = add_pval_comparison(aIns_Lac_males,...
        aIns_Lac_females,...
        aIns_Lac.pval, 4, 5, 'NS');
    % plasma
    [l_hdl, star_hdl] = add_pval_comparison(plasmaLac_males,...
        plasmaLac_females,...
        plasma_Lac.pval, 7, 8, 'NS');
    ylabel('Lactate (mM)');
    xticks([1,2, 4,5, 7,8]);
    xticklabels({'M','F','M','F','M','F'});
end % figure display
    
end % function