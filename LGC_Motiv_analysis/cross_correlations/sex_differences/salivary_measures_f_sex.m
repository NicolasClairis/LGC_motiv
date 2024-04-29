function[CORT, TESTO, TESTOCORT, IL] = salivary_measures_f_sex(rmv_outliers_yn)
% [CORT, TESTO, TESTOCORT, IL] = salivary_measures_f_sex(rmv_outliers_yn)
% salivary_measures_f_sex will load testosterone, cortisol and inflammatory
% interleukin concentrations measured in the saliva with ELISA kits and
% compare whether there are differences between male and females in these
% measures.
%
% INPUTS
% rmv_outliers_yn: remove outliers (1) or not (0)?
%
% OUTPUTS
% CORT: structure with cortisol mean, sem and p.value for comparison
% between males and females
%
% TESTO: structure with testosterone mean, sem and p.value for comparison
% between males and females
%
% TESTOCORT: structure with testosterone/cortisol ratio mean, sem and p.value for comparison
% between males and females
%
% IL: structure with IL-1b, IL-6 and IL-18 mean, sem and p.value for comparison
% between males and females

%% define default inputs
% remove outliers by default
if ~exist('rmv_outliers_yn','var') || isempty(rmv_outliers_yn) || ~ismember(rmv_outliers_yn,[0,1])
    rmv_outliers_yn = 1;
end

%% subject selection
study_nm = 'study1';
[male_CIDS, female_CIDS, male_NS, female_NS, condition] = subject_selection_per_sex;

%% load all salivary measures
n_TestoCort_samples = 4;

% testosterone
[TESTO_data_males] = load_TESTO(study_nm, male_CIDS);
[TESTO_data_females] = load_TESTO(study_nm, female_CIDS);
TESTO_males = TESTO_data_males.TESTO;
TESTO_females = TESTO_data_females.TESTO;

% cortisol
[CORT_data_males] = load_CORT(study_nm, male_CIDS);
[CORT_data_females] = load_CORT(study_nm, female_CIDS);
CORT_males = CORT_data_males.CORT;
CORT_females = CORT_data_females.CORT;

% extract TESTO/CORT ratio after converting in same unit (μg/dL)
TestoCort_males = (TESTO_males.*0.0001)./CORT_males;
TestoCort_females = (TESTO_females.*0.0001)./CORT_females;

% interleukins
[IL_data_males] = load_IL(study_nm, male_CIDS);
[IL_data_females] = load_IL(study_nm, female_CIDS);
IL_1b_males = IL_data_males.IL1b';
IL_1b_females = IL_data_females.IL1b';
IL_6_males = IL_data_males.IL6';
IL_6_females = IL_data_females.IL6';
IL_18_males = IL_data_males.IL18';
IL_18_females = IL_data_females.IL18';

%% remove outliers
if rmv_outliers_yn == 1
    for iSample = 1:n_TestoCort_samples
        % testosterone
        [~,~,TESTO_males(iSample,:)] = rmv_outliers_3sd(TESTO_males(iSample,:));
        [~,~,TESTO_females(iSample,:)] = rmv_outliers_3sd(TESTO_females(iSample,:));
        
        % cortisol
        [~,~,CORT_males(iSample,:)] = rmv_outliers_3sd(CORT_males(iSample,:));
        [~,~,CORT_females(iSample,:)] = rmv_outliers_3sd(CORT_females(iSample,:));
        
        % testosterone/cortisol
        [~,~,TestoCort_males(iSample,:)] = rmv_outliers_3sd(TestoCort_males(iSample,:));
        [~,~,TestoCort_females(iSample,:)] = rmv_outliers_3sd(TestoCort_females(iSample,:));
    end % loop over samples
    
    % interleukins (only 1 time sample)
    % IL-1b
    [~,~,IL_1b_males] = rmv_outliers_3sd(IL_1b_males);
    [~,~,IL_1b_females] = rmv_outliers_3sd(IL_1b_females);
    % IL-6
    [~,~,IL_6_males] = rmv_outliers_3sd(IL_6_males);
    [~,~,IL_6_females] = rmv_outliers_3sd(IL_6_females);
    % IL-18
    [~,~,IL_18_males] = rmv_outliers_3sd(IL_18_males);
    [~,~,IL_18_females] = rmv_outliers_3sd(IL_18_females);
end % outlier removal

%% extract mean + SEM per sex and perform the comparison
[TESTO.pval, CORT.pval, TESTOCORT.pval,...
    TESTO.m_males, CORT.m_males, TESTOCORT.m_males,...
    TESTO.m_females, CORT.m_females, TESTOCORT.m_females,...
    TESTO.sem_males, CORT.sem_males, TESTOCORT.sem_males,...
    TESTO.sem_females, CORT.sem_females, TESTOCORT.sem_females] = deal(NaN(n_TestoCort_samples, 1));

for iSample = 1:n_TestoCort_samples
    % testosterone
    [~,TESTO.pval(iSample)] = ttest2(TESTO_males(iSample,:), TESTO_females(iSample,:));
    [TESTO.m_males(iSample), TESTO.sem_males(iSample)] = mean_sem_sd(TESTO_males(iSample,:),2);
    [TESTO.m_females(iSample), TESTO.sem_females(iSample)] = mean_sem_sd(TESTO_females(iSample,:),2);
    
    % cortisol
    [~,CORT.pval(iSample)] = ttest2(CORT_males(iSample,:), CORT_females(iSample,:));
    [CORT.m_males(iSample), CORT.sem_males(iSample)] = mean_sem_sd(CORT_males(iSample,:),2);
    [CORT.m_females(iSample), CORT.sem_females(iSample)] = mean_sem_sd(CORT_females(iSample,:),2);
    
    % testosterone/cortisol
    [~,TESTOCORT.pval(iSample)] = ttest2(TestoCort_males(iSample,:), TestoCort_females(iSample,:));
    [TESTOCORT.m_males(iSample), TESTOCORT.sem_males(iSample)] = mean_sem_sd(TestoCort_males(iSample,:),2);
    [TESTOCORT.m_females(iSample), TESTOCORT.sem_females(iSample)] = mean_sem_sd(TestoCort_females(iSample,:),2);
end % loop through testosterone and cortisol samples

% interleukins
% IL 1b
[~,IL.IL_1b.pval] = ttest2(IL_1b_males, IL_1b_females);
[IL.IL_1b.m_males, IL.IL_1b.sem_males] = mean_sem_sd(IL_1b_males,2);
[IL.IL_1b.m_females, IL.IL_1b.sem_females] = mean_sem_sd(IL_1b_females,2);
% IL 6
[~,IL.IL_6.pval] = ttest2(IL_6_males, IL_6_females);
[IL.IL_6.m_males, IL.IL_6.sem_males] = mean_sem_sd(IL_6_males,2);
[IL.IL_6.m_females, IL.IL_6.sem_females] = mean_sem_sd(IL_6_females,2);
% IL 18
[~,IL.IL_18.pval] = ttest2(IL_18_males, IL_18_females);
[IL.IL_18.m_males, IL.IL_18.sem_males] = mean_sem_sd(IL_18_males,2);
[IL.IL_18.m_females, IL.IL_18.sem_females] = mean_sem_sd(IL_18_females,2);

%% figure
% colour code
[pSize, lW, col, mSize] = general_fig_prm;
female_col = col.red;
male_col = col.blue_dark;

% display figure for each hormone
fig;

% testosterone
subplot(1,3,1);
for iSample = 1:n_TestoCort_samples
    jPos_male = 1 + 2*(iSample - 1);
    jPos_female = 2 + 2*(iSample - 1);
    
    % show male vs female data
    ok_males = ~isnan(TESTO_males(iSample,:));
    male_violin = Violin({TESTO_males(iSample,ok_males)},jPos_male,...
        'ViolinColor',{male_col});
    ok_females = ~isnan(TESTO_females(iSample,:));
    female_violin = Violin({TESTO_females(iSample,ok_females)},jPos_female,...
        'ViolinColor',{female_col});
    
    % add p.value indication if difference is significant
    if ismember(iSample,[1,4]) % pointless for samples 2 and 3 not measured
        [l_hdl, star_hdl] = add_pval_comparison(TESTO_males(iSample,:),...
            TESTO_females(iSample,:),...
            TESTO.pval(iSample), jPos_male, jPos_female, 'NS');
    end
end % loop through testosterone and cortisol samples
xticks(1.5:2:n_TestoCort_samples*2);
xticklabels({'1','2','3','4'});
ylabel('Testosterone (pg/mL)');
xlabel('Sample');


% cortisol
subplot(1,3,2);
for iSample = 1:n_TestoCort_samples
    jPos_male = 1 + 2*(iSample - 1);
    jPos_female = 2 + 2*(iSample - 1);
    
    % show male vs female data
    ok_males = ~isnan(CORT_males(iSample,:));
    male_violin = Violin({CORT_males(iSample,ok_males)},jPos_male,...
        'ViolinColor',{male_col});
    ok_females = ~isnan(CORT_females(iSample,:));
    female_violin = Violin({CORT_females(iSample,ok_females)},jPos_female,...
        'ViolinColor',{female_col});
    
    % add p.value indication if difference is significant
    [l_hdl, star_hdl] = add_pval_comparison(CORT_males(iSample,:),...
        CORT_females(iSample,:),...
        CORT.pval(iSample), jPos_male, jPos_female, 'NS');
end % loop through testosterone and cortisol samples
xticks(1.5:2:n_TestoCort_samples*2);
xticklabels({'1','2','3','4'});
ylabel('Cortisol (μg/dL)');
xlabel('Sample');

% testosterone/cortisol
subplot(1,3,3);
for iSample = 1:n_TestoCort_samples
    jPos_male = 1 + 2*(iSample - 1);
    jPos_female = 2 + 2*(iSample - 1);
    
    % show male vs female data
    ok_males = ~isnan(TestoCort_males(iSample,:));
    male_violin = Violin({TestoCort_males(iSample,ok_males)},jPos_male,...
        'ViolinColor',{male_col});
    ok_females = ~isnan(TestoCort_females(iSample,:));
    female_violin = Violin({TestoCort_females(iSample,ok_females)},jPos_female,...
        'ViolinColor',{female_col});
    
    % add p.value indication if difference is significant
    if ismember(iSample,[1,4]) % pointless for samples 2 and 3 not measured
        [l_hdl, star_hdl] = add_pval_comparison(TestoCort_males(iSample,:),...
            TestoCort_females(iSample,:),...
            TESTOCORT.pval(iSample), jPos_male, jPos_female, 'NS');
    end
end % loop through testosterone and cortisol samples
xticks(1.5:2:n_TestoCort_samples*2);
xticklabels({'1','2','3','4'});
ylabel('Testosterone/Cortisol');
xlabel('Sample');


% interleukins
fig;

% IL-1b
subplot(1,3,1);
jPos_male = 1;
jPos_female = 2;

% show male vs female data
ok_males = ~isnan(IL_1b_males);
male_violin = Violin({IL_1b_males(ok_males)},jPos_male,...
    'ViolinColor',{male_col});
ok_females = ~isnan(IL_1b_females);
female_violin = Violin({IL_1b_females(ok_females)},jPos_female,...
    'ViolinColor',{female_col});

% add p.value indication if difference is significant
[l_hdl, star_hdl] = add_pval_comparison(IL_1b_males,...
    IL_1b_females,...
    IL.IL_1b.pval, jPos_male, jPos_female, 'NS');
ylabel('IL-1β (pg/mL)');

% IL-6
subplot(1,3,2);
jPos_male = 1;
jPos_female = 2;

% show male vs female data
ok_males = ~isnan(IL_6_males);
male_violin = Violin({IL_6_males(ok_males)},jPos_male,...
    'ViolinColor',{male_col});
ok_females = ~isnan(IL_6_females);
female_violin = Violin({IL_6_females(ok_females)},jPos_female,...
    'ViolinColor',{female_col});

% add p.value indication if difference is significant
[l_hdl, star_hdl] = add_pval_comparison(IL_6_males,...
    IL_6_females,...
    IL.IL_6.pval, jPos_male, jPos_female, 'NS');
ylabel('IL-6 (pg/mL)');

% IL-18
subplot(1,3,3);
jPos_male = 1;
jPos_female = 2;

% show male vs female data
ok_males = ~isnan(IL_18_males);
male_violin = Violin({IL_18_males(ok_males)},jPos_male,...
    'ViolinColor',{male_col});
ok_females = ~isnan(IL_18_females);
female_violin = Violin({IL_18_females(ok_females)},jPos_female,...
    'ViolinColor',{female_col});

% add p.value indication if difference is significant
[l_hdl, star_hdl] = add_pval_comparison(IL_18_males,...
    IL_18_females,...
    IL.IL_18.pval, jPos_male, jPos_female, 'NS');
ylabel('IL-18 (pg/mL)');


end % function