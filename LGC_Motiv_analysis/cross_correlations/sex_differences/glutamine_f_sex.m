function[plasma_Gln, dmPFC_Gln, aIns_Gln] = glutamine_f_sex(fig_disp, rmv_outliers_yn)
% [plasma_Gln, dmPFC_Gln, aIns_Gln] = glutamine_f_sex(fig_disp, rmv_outliers_yn)
% glutamine_f_sex compares plasma, dmPFC/dACC and aIns glutamine levels in
% males vs females.
% INPUTS
% fig_disp: display figures? (1) yes (0) no
%
% rmv_outliers_yn: remove median +/- 3*SD outliers yes (1) or no (0)? Yes by
% default
%
% OUTPUTS
% plasma_Gln, dmPFC_Gln, aIns_Gln: structure with p.value (for unpaired
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

%% extract glutamine levels
% plasma
plasmaM_males = load_plasma_metabolites(male_CIDS);
plasmaM_females = load_plasma_metabolites(female_CIDS);
plasmaGln_males = plasmaM_males.Gln./1000;
plasmaGln_females = plasmaM_females.Gln./1000;

% load brain metabolites
metabolites_males = metabolite_load(male_CIDS);
metabolites_females = metabolite_load(female_CIDS);
dmPFC_Gln_males = metabolites_males.dmPFC.Gln;
dmPFC_Gln_females = metabolites_females.dmPFC.Gln;
aIns_Gln_males = metabolites_males.aIns.Gln;
aIns_Gln_females = metabolites_females.aIns.Gln;

%% remove outliers
if rmv_outliers_yn == 1
    % plasma
    [~,~,plasmaGln_males] = rmv_outliers_3sd(plasmaGln_males);
    [~,~,plasmaGln_females] = rmv_outliers_3sd(plasmaGln_females);
    % brain
    [~,~,dmPFC_Gln_males] = rmv_outliers_3sd(dmPFC_Gln_males);
    [~,~,dmPFC_Gln_females] = rmv_outliers_3sd(dmPFC_Gln_females);
    [~,~,aIns_Gln_males] = rmv_outliers_3sd(aIns_Gln_males);
    [~,~,aIns_Gln_females] = rmv_outliers_3sd(aIns_Gln_females);
end % outliers

%% extract mean, SEM per sex and perform the comparison
% plasma
[~,plasma_Gln.pval] = ttest2(plasmaGln_males, plasmaGln_females);
[plasma_Gln.m_males, plasma_Gln.sem_males,...
    plasma_Gln.sd_males, plasma_Gln.med_males] = mean_sem_sd(plasmaGln_males,2);
[plasma_Gln.m_females, plasma_Gln.sem_females,...
    plasma_Gln.sd_females, plasma_Gln.med_females] = mean_sem_sd(plasmaGln_females,2);
% dmPFC
[~,dmPFC_Gln.pval] = ttest2(dmPFC_Gln_males, dmPFC_Gln_females);
[dmPFC_Gln.m_males, dmPFC_Gln.sem_males,...
    dmPFC_Gln.sd_males, dmPFC_Gln.med_males] = mean_sem_sd(dmPFC_Gln_males,2);
[dmPFC_Gln.m_females, dmPFC_Gln.sem_females,...
    dmPFC_Gln.sd_females, dmPFC_Gln.med_females] = mean_sem_sd(dmPFC_Gln_females,2);
% aIns
[~,aIns_Gln.pval] = ttest2(aIns_Gln_males, aIns_Gln_females);
[aIns_Gln.m_males, aIns_Gln.sem_males,...
    aIns_Gln.sd_males, aIns_Gln.med_males] = mean_sem_sd(aIns_Gln_males,2);
[aIns_Gln.m_females, aIns_Gln.sem_females,...
    aIns_Gln.sd_females, aIns_Gln.med_females] = mean_sem_sd(aIns_Gln_females,2);

%% figure display
if fig_disp == 1
    % general figure parameters
    [pSize, lW, col, mSize] = general_fig_prm;
    female_col = col.red;
    male_col = col.blue_dark;
    
    % glutamine
    fig;
    
    % show male vs female data
    % dmPFC
    ok_dmPFC_males = ~isnan(dmPFC_Gln_males);
    Violin({dmPFC_Gln_males(ok_dmPFC_males)},1,...
        'ViolinColor',{male_col});
    ok_dmPFC_females = ~isnan(dmPFC_Gln_females);
    Violin({dmPFC_Gln_females(ok_dmPFC_females)},2,...
        'ViolinColor',{female_col});
    % aIns
    ok_aIns_males = ~isnan(aIns_Gln_males);
    Violin({aIns_Gln_males(ok_aIns_males)},4,...
        'ViolinColor',{male_col});
    ok_aIns_females = ~isnan(aIns_Gln_females);
    Violin({aIns_Gln_females(ok_aIns_females)},5,...
        'ViolinColor',{female_col});
    % plasma
    ok_plasma_males = ~isnan(plasmaGln_males);
    Violin({plasmaGln_males(ok_plasma_males)},7,...
        'ViolinColor',{male_col});
    ok_plasma_females = ~isnan(plasmaGln_females);
    Violin({plasmaGln_females(ok_plasma_females)},8,...
        'ViolinColor',{female_col});
    
    % add p.value indication if difference is significant
    % dmPFC
    [l_hdl, star_hdl] = add_pval_comparison(dmPFC_Gln_males,...
        dmPFC_Gln_females,...
        dmPFC_Gln.pval, 1, 2, 'NS');
    % aIns
    [l_hdl, star_hdl] = add_pval_comparison(aIns_Gln_males,...
        aIns_Gln_females,...
        aIns_Gln.pval, 4, 5, 'NS');
    % plasma
    [l_hdl, star_hdl] = add_pval_comparison(plasmaGln_males,...
        plasmaGln_females,...
        plasma_Gln.pval, 7, 8, 'NS');
    ylabel('Glutamine (mM)');
    xticks([1,2, 4,5, 7,8]);
    xticklabels({'M','F','M','F','M','F'});
end % figure display
    
end % function