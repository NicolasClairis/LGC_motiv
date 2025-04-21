function[plasma_aa, plasma_Lac, plasma_FA] = plasma_f_sex(fig_disp, rmv_outliers_yn)
% [plasma_aa, plasma_Lac, plasma_FA] = plasma_f_sex(fig_disp, rmv_outliers_yn)
%
% INPUTS
% fig_disp: display figures? (1) yes (0) no
%
% rmv_outliers_yn: remove median +/- 3*SD outliers yes (1) or no (0)? Yes by
% default
%
% OUTPUTS
% plasma_aa: structure with mean, sem and p.value for comparison
% between males and females for different amino-acids (aa) measured in the
% plasma
%
% plasma_Lac: structure with mean, sem and p.value for comparison
% between males and females for lactate measured in the plasma
%
% plasma_FA: structure with mean, sem and p.value for comparison between 
% males and females for different fatty acids measured in the plasma

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

%% extract plasma metabolites
[plasmaM_males, mb_names, n_mb] = load_plasma_metabolites(male_CIDS);
plasmaM_females = load_plasma_metabolites(female_CIDS);

%% split metabolites per category (aa/fatty acids/Lactate)
% extract amino-acids
aa_names = {'Pro','Gly','Ile','Asn','Tyr','Tau','Ser','His','Asp',...
    'Ala','Arg','Thr','Leu','Lys','Gln','Val','Glu','Phe','Met',...
    'Trans_4_Hydroxyproline',...
    'SymmetricDimethylarginine_SDMA_','AsymmetricDimethylarginine_ADMA_',...
    'L_Citrulline','Sarcosine','x3_MethylHistidine','x1_Methylhistidine'};
n_aa = length(aa_names);
for iAA = 1:n_aa
    aa_nm = aa_names{iAA};
    aa.(aa_nm).males = plasmaM_males.(aa_nm);
    aa.(aa_nm).females = plasmaM_females.(aa_nm);
end % loop over amino-acids

% extract lactate
Lac_males = plasmaM_males.Lac;
Lac_females = plasmaM_females.Lac;

% extract fatty acids
fatty_acids_names = {'Acetoacetate','Hydroxybutyrate','Acetate',...
    'Propionate','Isobutyrate','Butyrate','x2_Methylbutyrate',...
    'Isovalerate','Valerate','HexanoicAcid','OctanoicAcid','DecanoicAcid',...
    'a_AminobutyricA_AABA_'};
n_fattyA = length(fatty_acids_names);
for iFA = 1:n_fattyA
    FA_nm = fatty_acids_names{iFA};
    FA.(FA_nm).males = plasmaM_males.(FA_nm);
    FA.(FA_nm).females = plasmaM_females.(FA_nm);
end % loop over amino-acids

%% remove outliers
if rmv_outliers_yn == 1
    % amino-acids
    for iAA = 1:n_aa
        aa_nm = aa_names{iAA};
        [~,~,aa.(aa_nm).males] = rmv_outliers_3sd(aa.(aa_nm).males);
        [~,~,aa.(aa_nm).females] = rmv_outliers_3sd(aa.(aa_nm).females);
    end % loop over amino-acids
    
    % lactate
    [~,~,Lac_males] = rmv_outliers_3sd(Lac_males);
    [~,~,Lac_females] = rmv_outliers_3sd(Lac_females);
    
    % fatty acids
    for iFA = 1:n_fattyA
        FA_nm = fatty_acids_names{iFA};
        [~,~,FA.(FA_nm).males] = rmv_outliers_3sd(FA.(FA_nm).males);
        [~,~,FA.(FA_nm).females] = rmv_outliers_3sd(FA.(FA_nm).females);
    end % loop over amino-acids
end % remove outliers?

%% extract mean, SEM per sex and perform the comparison
% amino-acids
for iAA = 1:n_aa
    aa_nm = aa_names{iAA};
    [~,plasma_aa.(aa_nm).pval] = ttest2(aa.(aa_nm).males, aa.(aa_nm).females);
    [plasma_aa.(aa_nm).m_males, plasma_aa.(aa_nm).sem_males] = mean_sem_sd(aa.(aa_nm).males,2);
    [plasma_aa.(aa_nm).m_females, plasma_aa.(aa_nm).sem_females] = mean_sem_sd(aa.(aa_nm).females,2);
    % note significant amino-acids
    if plasma_aa.(aa_nm).pval < 0.05
        plasma_aa.signif.(aa_nm).pval = plasma_aa.(aa_nm).pval;
        plasma_aa.signif.(aa_nm).m_males = plasma_aa.(aa_nm).m_males;
        plasma_aa.signif.(aa_nm).m_females = plasma_aa.(aa_nm).m_females;
    end 
end % amino-acids loop

% lactate
[~,plasma_Lac.pval] = ttest2(Lac_males, Lac_females);
[plasma_Lac.m_males, plasma_Lac.sem_males] = mean_sem_sd(Lac_males,2);
[plasma_Lac.m_females, plasma_Lac.sem_females] = mean_sem_sd(Lac_females,2);

% fatty acids
for iFA = 1:n_fattyA
    FA_nm = fatty_acids_names{iFA};
    [~,plasma_FA.(FA_nm).pval] = ttest2(FA.(FA_nm).males, FA.(FA_nm).females);
    [plasma_FA.(FA_nm).m_males, plasma_FA.(FA_nm).sem_males] = mean_sem_sd(FA.(FA_nm).males,2);
    [plasma_FA.(FA_nm).m_females, plasma_FA.(FA_nm).sem_females] = mean_sem_sd(FA.(FA_nm).females,2);
    % note significant amino-acids
    if plasma_FA.(FA_nm).pval < 0.05
        plasma_FA.signif.(FA_nm).pval = plasma_FA.(FA_nm).pval;
        plasma_FA.signif.(FA_nm).m_males = plasma_FA.(FA_nm).m_males;
        plasma_FA.signif.(FA_nm).m_females = plasma_FA.(FA_nm).m_females;
    end 
end % fatty acids loop

%% figure display
if fig_disp == 1
    % general figure parameters
    [pSize, lW, col, mSize] = general_fig_prm;
    female_col = col.red;
    male_col = col.blue_dark;

    % amino-acids
    fig;
    for iAA = 1:n_aa
        aa_nm = aa_names{iAA};
        
        % x-coordinates
        jPos_male = 1 + 2*(iAA - 1);
        jPos_female = 2 + 2*(iAA - 1);
        
        % show male vs female data
        ok_males = ~isnan(aa.(aa_nm).males);
        male_violin = Violin({aa.(aa_nm).males(ok_males)},jPos_male,...
            'ViolinColor',{male_col});
        ok_females = ~isnan(aa.(aa_nm).females);
        female_violin = Violin({aa.(aa_nm).females(ok_females)},jPos_female,...
            'ViolinColor',{female_col});
        
        % add p.value indication if difference is significant
        [l_hdl, star_hdl] = add_pval_comparison(aa.(aa_nm).males,...
            aa.(aa_nm).females,...
            plasma_aa.(aa_nm).pval, jPos_male, jPos_female, 'NS');
    end % amino-acids
    xticks(1.5:2:n_aa*2);
    aa_names_bis = aa_names;
    for iAA = 1:n_aa
        if length(aa_names_bis{iAA}) > 3
            aa_names_bis{iAA} = aa_names_bis{iAA}(1:3);
        end
    end
    xticklabels(aa_names_bis);
    ylabel('Concentration (μM)');
    
    
    % Lactate
    fig;
    % x-coordinates
    jPos_male = 1;
    jPos_female = 2;
    
    % show male vs female data
    ok_males = ~isnan(Lac_males);
    male_violin = Violin({Lac_males(ok_males)},jPos_male,...
        'ViolinColor',{male_col});
    ok_females = ~isnan(Lac_females);
    female_violin = Violin({Lac_females(ok_females)},jPos_female,...
        'ViolinColor',{female_col});
    
    % add p.value indication if difference is significant
    [l_hdl, star_hdl] = add_pval_comparison(Lac_males,...
        Lac_females,...
        plasma_Lac.pval, jPos_male, jPos_female, 'NS');
    ylabel('Lactate (μM)');
    xticks(1:2);
    xticklabels({'Males','Females'});
    
    % fatty acids
    fig;
    for iFA = 1:n_fattyA
        FA_nm = fatty_acids_names{iFA};
        
        % x-coordinates
        jPos_male = 1 + 2*(iFA - 1);
        jPos_female = 2 + 2*(iFA - 1);
        
        % show male vs female data
        ok_males = ~isnan(FA.(FA_nm).males);
        male_violin = Violin({FA.(FA_nm).males(ok_males)},jPos_male,...
            'ViolinColor',{male_col});
        ok_females = ~isnan(FA.(FA_nm).females);
        female_violin = Violin({FA.(FA_nm).females(ok_females)},jPos_female,...
            'ViolinColor',{female_col});
        
        % add p.value indication if difference is significant
        [l_hdl, star_hdl] = add_pval_comparison(FA.(FA_nm).males,...
            FA.(FA_nm).females,...
            plasma_FA.(FA_nm).pval, jPos_male, jPos_female, 'NS');
    end % amino-acids
    xticks(1.5:2:n_fattyA*2);
    FA_names_bis = fatty_acids_names;
    for iFA = 1:n_fattyA
        if length(FA_names_bis{iFA}) > 3
            FA_names_bis{iFA} = FA_names_bis{iFA}(1:3);
        end
    end
    xticklabels(FA_names_bis);
    ylabel('Concentration (μM)');
    
end % figure display

end % function