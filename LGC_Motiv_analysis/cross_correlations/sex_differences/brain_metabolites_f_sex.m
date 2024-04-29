function[brain_mb] = brain_metabolites_f_sex(fig_disp, rmv_outliers_yn)
% [brain_mb] = brain_metabolites_f_sex(fig_disp, rmv_outliers_yn)
% brain_metabolites_f_sex will extract the concentration of brain
% metabolites in the dmPFC/dACC and in the anterior insula measured with
% 1H-MRS and test if the concentrations are significantly different between
% male and female.
%
% INPUTS
% fig_disp: display figures? (1) yes (0) no
%
% rmv_outliers_yn: remove median +/- 3*SD outliers yes (1) or no (0)? Yes by
% default
%
% OUTPUTS
% brain_mb: structure with mean, sem and p.value for comparison
% between males and females for the brain metabolite concentrations in the
% dorsomedial prefrontal cortex/dorsal anterior cingulate cortex
% (dmPFC/dACC) and the anterior insula (aINS)

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

%% load brain metabolites
metabolites_males = metabolite_load(male_CIDS);
metabolites_females = metabolite_load(female_CIDS);
mb_names = fieldnames(metabolites_males.dmPFC);
n_mb = length(mb_names);
brain_areas = {'dmPFC','aIns'};
n_brain_areas = length(brain_areas);

%% remove outliers
if rmv_outliers_yn == 1
    for iBrain = 1:n_brain_areas
        brain_area_nm = brain_areas{iBrain};
        for iMb = 1:n_mb
            mb_nm = mb_names{iMb};
            [~,~,metabolites_males.(brain_area_nm).(mb_nm)] = rmv_outliers_3sd(metabolites_males.(brain_area_nm).(mb_nm));
            [~,~,metabolites_females.(brain_area_nm).(mb_nm)] = rmv_outliers_3sd(metabolites_females.(brain_area_nm).(mb_nm));
        end % metabolite
    end % brain area
end % outlier removal

%% extract mean, SEM per sex and perform the comparison
for iBrain = 1:n_brain_areas
    brain_area_nm = brain_areas{iBrain};
    for iMb = 1:n_mb
        mb_nm = mb_names{iMb};
        [~,brain_mb.(brain_area_nm).(mb_nm).pval] = ttest2(metabolites_males.(brain_area_nm).(mb_nm),...
            metabolites_females.(brain_area_nm).(mb_nm));
        [brain_mb.(brain_area_nm).(mb_nm).m_males,...
            brain_mb.(brain_area_nm).(mb_nm).sem_males] = mean_sem_sd(metabolites_males.(brain_area_nm).(mb_nm),2);
        [brain_mb.(brain_area_nm).(mb_nm).m_females,...
            brain_mb.(brain_area_nm).(mb_nm).sem_females] = mean_sem_sd(metabolites_females.(brain_area_nm).(mb_nm),2);
        
    % extract significant data
    if brain_mb.(brain_area_nm).(mb_nm).pval < 0.05
        brain_mb.(brain_area_nm).signif.(mb_nm).pval = brain_mb.(brain_area_nm).(mb_nm).pval;
        brain_mb.(brain_area_nm).signif.(mb_nm).m_males = brain_mb.(brain_area_nm).(mb_nm).m_males;
        brain_mb.(brain_area_nm).signif.(mb_nm).m_females = brain_mb.(brain_area_nm).(mb_nm).m_females;
    end % p<0.05 filter
    end % metabolite
end % brain area

%% figure display
if fig_disp == 1
    % general figure parameters
    [pSize, lW, col, mSize] = general_fig_prm;
    female_col = col.red;
    male_col = col.blue_dark;
    
    % global figures (split data in 2 because too many metabolites)
    fig1 = fig;
    fig2 = fig;
    n_half_mb = n_mb/2;
    
    % loop over brain areas
    for iBrain = 1:n_brain_areas
        brain_area_nm = brain_areas{iBrain};
        figure(fig1);
        subplot(1,n_brain_areas,iBrain);
        figure(fig2);
        subplot(1,n_brain_areas,iBrain);
        
        % rename variables to ease the visual display
        mb_names_bis = mb_names;
        for iMb = 1:n_mb
            mb_nm = mb_names{iMb};
            switch mb_nm
                case 'GPC_PCho'
                    mb_names_bis{iMb} = 'GPC+PCho';
                case 'Cr_PCr'
                    mb_names_bis{iMb} = 'Cr+Pcr';
                case 'Gln_div_Glu'
                    mb_names_bis{iMb} = 'Gln/Glu';
                case 'z_antiox'
                    mb_names_bis{iMb} = 'z(GSH)+z(Tau)';
                case 'antiox'
                    mb_names_bis{iMb} = 'GSH+Tau';
                case 'Glu_div_GSH'
                    mb_names_bis{iMb} = 'Glu/GSH';
                case 'Glu_div_Tau'
                    mb_names_bis{iMb} = 'Glu/Tau';
                case 'Glu_div_antiox'
                    mb_names_bis{iMb} = 'Glu/(GSH+Tau)';
                case 'Glu_div_z_antiox'
                    mb_names_bis{iMb} = 'Glu/(z(GSH)+z(Tau))';
                case 'Glu_div_GABA_div_GSH'
                    mb_names_bis{iMb} = '(Glu/GABA)/GSH';
                case 'Glu_div_GABA'
                    mb_names_bis{iMb} = 'Glu/GABA';
                case 'Asp_plus_Lac'
                    mb_names_bis{iMb} = 'Asp+Lac';
                case 'zAsp_plus_zLac'
                    mb_names_bis{iMb} = 'z(Asp)+z(Lac)';
                case 'Glu_Ushape'
                    mb_names_bis{iMb} = '(Glu-μ(Glu))²';
                case 'GSH_Ushape'
                    mb_names_bis{iMb} = '(GSH-μ(GSH))²';
                case 'Glu_Ushape_div_GSH'
                    mb_names_bis{iMb} = '(Glu-μ(Glu))²/GSH';
                case 'Glu_div_GSH_Ushape'
                    mb_names_bis{iMb} = 'Glu/(GSH-μ(GSH))²';
                case 'Glu_Ushape_div_GSH_Ushape'
                    mb_names_bis{iMb} = '(Glu-μ(Glu))²/(GSH-μ(GSH))²';
                case 'Glu_div_GSH_ratio_Ushape'
                    mb_names_bis{iMb} = '((Glu/GSH)-μ(Glu/GSH))²';
                case 'Glu_div_Asp'
                    mb_names_bis{iMb} = 'Glu/Asp';
                case 'Lac_div_GSH'
                    mb_names_bis{iMb} = 'Lac/GSH';
                case 'Lac_div_antiox'
                    mb_names_bis{iMb} = 'Lac/(GSH+Tau)';
                case 'Glu_plus_Lac_div_GSH'
                    mb_names_bis{iMb} = '(Glu+Lac)/GSH';
                case 'Glu_plus_Lac_div_antiox'
                    mb_names_bis{iMb} = '(Glu+Lac)/(GSH+Tau)';
                otherwise
                    mb_names_bis{iMb} = mb_names{iMb};
            end
        end % blood metabolite
        
        % split metabolite names for the 2 figures
        mb_names_bis1 = mb_names_bis(1:n_half_mb);
        mb_names_bis2 = mb_names_bis((1:n_half_mb)+n_half_mb);
        
        % add violin plots with data
        for iMb = 1:n_mb
            mb_nm = mb_names{iMb};
            
            if iMb <= n_half_mb
                figure(fig1);
                subplot(1,n_brain_areas,iBrain);
                % x-coordinates
                jPos_male = 1 + 2*(iMb - 1);
                jPos_female = 2 + 2*(iMb - 1);
            else
                figure(fig2);
                subplot(1,n_brain_areas,iBrain);
                % x-coordinates
                % reset x_coordinate
                jPos_male = 1 + 2*(iMb-n_half_mb - 1);
                jPos_female = 2 + 2*(iMb-n_half_mb - 1);
            end
            
            % show male vs female data
            ok_males = ~isnan(metabolites_males.(brain_area_nm).(mb_nm));
            male_violin = Violin({metabolites_males.(brain_area_nm).(mb_nm)(ok_males)},jPos_male,...
                'ViolinColor',{male_col});
            ok_females = ~isnan(metabolites_females.(brain_area_nm).(mb_nm));
            female_violin = Violin({metabolites_females.(brain_area_nm).(mb_nm)(ok_females)},jPos_female,...
                'ViolinColor',{female_col});
        end % loop over metabolites
        
        % add p.values
        for iMb = 1:n_mb
            mb_nm = mb_names{iMb};
            
            if iMb <= n_half_mb
                figure(fig1);
                subplot(1,n_brain_areas,iBrain);
                % x-coordinates
                jPos_male = 1 + 2*(iMb - 1);
                jPos_female = 2 + 2*(iMb - 1);
            else
                figure(fig2);
                subplot(1,n_brain_areas,iBrain);
                % x-coordinates
                % reset x_coordinate
                jPos_male = 1 + 2*(iMb-n_half_mb - 1);
                jPos_female = 2 + 2*(iMb-n_half_mb - 1);
            end
            
            
            % add p.value indication if difference is significant
            [l_hdl, star_hdl] = add_pval_comparison(metabolites_males.(brain_area_nm).(mb_nm),...
                metabolites_females.(brain_area_nm).(mb_nm),...
                brain_mb.(brain_area_nm).(mb_nm).pval, jPos_male, jPos_female, '');
        end % loop over metabolites
        
        for iFig = 1:2
            switch iFig
                case 1
                    figure(fig1);
                case 2
                    figure(fig2);
            end
            subplot(1,n_brain_areas,iBrain);
            xticks(1.5:2:n_half_mb*2);
            switch iFig
                case 1
                    xticklabels(mb_names_bis1);
                case 2
                    xticklabels(mb_names_bis2);
            end
            ylabel('Concentration (μM)');
            xlim([0 n_half_mb*2+1]);
            title(brain_area_nm);
        end
    end % brain area loop
end % figure
end % function