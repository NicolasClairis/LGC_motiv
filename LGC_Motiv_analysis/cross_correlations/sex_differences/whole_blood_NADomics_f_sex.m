function[whole_blood_NAD] = whole_blood_NADomics_f_sex(fig_disp, rmv_outliers_yn)
% [whole_blood_NAD] = whole_blood_NADomics_f_sex(fig_disp, rmv_outliers_yn)
% whole_blood_NADomics_f_sex will test whether there are significant
% differences between males and females in the whole-blood NAD-related
% metabolomics.
%
% INPUTS
% fig_disp: display figures? (1) yes (0) no
%
% rmv_outliers_yn: remove median +/- 3*SD outliers yes (1) or no (0)? Yes by
% default
%
% OUTPUTS
% whole_blood_NAD: structure with mean, sem and p.value for comparison
% between males and females for the NAD-related molecules in the whole
% blood

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

%% load whole-blood NAD data
blood_males = load_blood_NAD(study_nm, male_CIDS);
blood_females = load_blood_NAD(study_nm, female_CIDS);

% identify blood metabolites
blood_mb_names = fieldnames(blood_males);
n_blood_mb = length(blood_mb_names);

%% remove outliers
if rmv_outliers_yn == 1
    for iB = 1:n_blood_mb
        blood_nm = blood_mb_names{iB};
        [~,~,blood_males.(blood_nm)] = rmv_outliers_3sd(blood_males.(blood_nm));
        [~,~,blood_females.(blood_nm)] = rmv_outliers_3sd(blood_females.(blood_nm));
    end % blood metabolite loop
end % outlier removal

%% extract mean, SEM per sex and perform the comparison
for iB = 1:n_blood_mb
    blood_nm = blood_mb_names{iB};
    [~,whole_blood_NAD.(blood_nm).pval] = ttest2(blood_males.(blood_nm), blood_females.(blood_nm));
    [whole_blood_NAD.(blood_nm).m_males, whole_blood_NAD.(blood_nm).sem_males] = mean_sem_sd(blood_males.(blood_nm),2);
    [whole_blood_NAD.(blood_nm).m_females, whole_blood_NAD.(blood_nm).sem_females] = mean_sem_sd(blood_females.(blood_nm),2);
    
    % extract significant data
    if whole_blood_NAD.(blood_nm).pval < 0.05
        whole_blood_NAD.signif.(blood_nm).pval = whole_blood_NAD.(blood_nm).pval;
        whole_blood_NAD.signif.(blood_nm).m_males = whole_blood_NAD.(blood_nm).m_males;
        whole_blood_NAD.signif.(blood_nm).m_females = whole_blood_NAD.(blood_nm).m_females;
    end % p<0.05 filter
end % blood metabolite loop

%% figure display
if fig_disp == 1
    % general figure parameters
    [pSize, lW, col, mSize] = general_fig_prm;
    female_col = col.red;
    male_col = col.blue_dark;
    
    % rename variables to ease the visual display
    blood_mb_names_bis = blood_mb_names;
    for iB = 1:n_blood_mb
        blood_nm = blood_mb_names{iB};
        switch blood_nm
            case 'NAD_div_NADH'
                blood_mb_names_bis{iB} = 'NAD/NADH';
            case 'NADP_div_NADPH'
                blood_mb_names_bis{iB} = 'NADP/NADPH';
            case 'total_NAD_precursors'
                blood_mb_names_bis{iB} = 'pre-NAD';
            case 'total_NAD'
                blood_mb_names_bis{iB} = 'tNAD';
            case 'total_NAD_with_precursors'
                blood_mb_names_bis{iB} = 'pre-NAD + tNAD';
            case 'total_NAD_with_byproducts'
                blood_mb_names_bis{iB} = 'tNAD + post-NAD';
            case 'total_NAD_byproducts'
                blood_mb_names_bis{iB} = 'post-NAD';
            otherwise
                blood_mb_names_bis{iB} = blood_mb_names{iB};
        end
    end % blood metabolite
    
    % global figure
    fig;
    
    % add violin plots with data
    for iB = 1:n_blood_mb
        blood_nm = blood_mb_names{iB};
        
        % x-coordinates
        jPos_male = 1 + 2*(iB - 1);
        jPos_female = 2 + 2*(iB - 1);
        
        % show male vs female data
        ok_males = ~isnan(blood_males.(blood_nm));
        male_violin = Violin({blood_males.(blood_nm)(ok_males)},jPos_male,...
            'ViolinColor',{male_col});
        ok_females = ~isnan(blood_females.(blood_nm));
        female_violin = Violin({blood_females.(blood_nm)(ok_females)},jPos_female,...
            'ViolinColor',{female_col});
    end % loop over blood metabolites
    
    % add p.values
    for iB = 1:n_blood_mb
        blood_nm = blood_mb_names{iB};
        
        % x-coordinates
        jPos_male = 1 + 2*(iB - 1);
        jPos_female = 2 + 2*(iB - 1);
        
        % add p.value indication if difference is significant
        [l_hdl, star_hdl] = add_pval_comparison(blood_males.(blood_nm),...
            blood_females.(blood_nm),...
            whole_blood_NAD.(blood_nm).pval, jPos_male, jPos_female, 'NS');
    end % loop over blood metabolites
    xticks(1.5:2:n_blood_mb*2);
    xticklabels(blood_mb_names_bis);
    ylabel('Concentration (μM)');
    xlim([0 n_blood_mb*2+1]);
end % figure
end % function