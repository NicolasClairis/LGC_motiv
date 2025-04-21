% brainMetabolites_f_stress_and_fatigue will check the correlation between
% stress, fatigue and brain metabolites.

%% outlier filtering
if ~exist('outlierF','var') || isempty(outlierF)
    outlierF_nm = questdlg('Outlier filtering?','Outlier filtering',...
        'No','Yes','Yes');
    switch outlierF_nm
        case 'Yes'
            outlierF = 1;
        case 'No'
            outlierF = 0;
    end
end % outlier filtering

%% subject selection
[study_nm, condition, subject_id, NS, genderFilter] = subject_selection;

%% load fatigue
[fatigue_measures] = fatigue_pool(study_nm, condition, subject_id, NS, genderFilter);
fatigue_vars = fieldnames(fatigue_measures);
fatigue_vars(strcmp(fatigue_vars,'sub_selection')) = [];
n_F_vars = length(fatigue_vars);
fatigue_vars_bis = cell(n_F_vars,1);
% rename fatigue variables for the figure
for iF = 1:n_F_vars
    [fatigue_vars_bis{iF}] = fatigue_nm_rename(fatigue_vars{iF});
end % loop over fatigue variables

%% load stress
[stress_measures] = stress_pool(study_nm, condition, subject_id, NS, genderFilter);
stress_vars = fieldnames(stress_measures);
stress_vars(strcmp(stress_vars,'sub_selection')) = [];
n_stress_vars = length(stress_vars);

%% load dmPFC/dACC and aIns Glu & Glu/GABA
[brainMetabolites, CRLB, brainMetabolites_bis] = metabolite_load(subject_id);

%% pool all variables
behavior_vars = [fatigue_measures.F1', fatigue_measures.F2',...
    fatigue_measures.F3', fatigue_measures.F4',...
    fatigue_measures.F4_F1', fatigue_measures.F4_F3',...
    fatigue_measures.MPSTEFS_mental_energy', fatigue_measures.MPSTEFS_physical_energy',...
    fatigue_measures.MPSTEFS_mental_fatigue', fatigue_measures.MPSTEFS_physical_fatigue',...
    stress_measures.STAI_T', stress_measures.SIAS',...
    stress_measures.PSS14',...
    stress_measures.S1', stress_measures.S2',...
    stress_measures.S3',stress_measures.S4',...
    stress_measures.S4_S1', stress_measures.S4_S3',...
    stress_measures.CORT1', stress_measures.CORT2',...
    stress_measures.CORT3', stress_measures.CORT4'];
bhv_var_names = {'F1','F2','F3','F4','F4_min_F1','F4_min_F3',...
    'MPSTEFS_m_energy','MPSTEFS_p_energy',...
    'MPSTEFS_m_fatigue','MPSTEFS_p_fatigue',...
    'STAIT','SIAS','PSS14',...
    'S1','S2','S3','S4','S4_min_S1','S4_min_S3',...
    'CORT1','CORT2','CORT3','CORT4'};
bhv_var_names_bis = strrep(bhv_var_names,'_min_','-');
bhv_var_names_bis = strrep(bhv_var_names_bis,'_',' ');
n_bhv_vars = size(behavior_vars,2);
brainMb_vars = [brainMetabolites.dmPFC.Glu', brainMetabolites.dmPFC.Glu_div_GABA',...
    brainMetabolites.aIns.Glu', brainMetabolites.aIns.Glu_div_GABA'];
brainMb_var_names = {'dmPFC_Glu','dmPFC_Glu_div_GABA',...
    'aIns_Glu','aIns_Glu_div_GABA'};
brainMb_var_names_bis = strrep(brainMb_var_names,'_div_','/');
brainMb_var_names_bis = strrep(brainMb_var_names_bis,'_',' ');
n_brainMb_vars = size(brainMb_vars,2);

%% correlations
[corr_mtrx, pval_mtrx] = deal(NaN(n_brainMb_vars, n_bhv_vars));
for iBhv = 1:n_bhv_vars
    for iBrainMb = 1:n_brainMb_vars
        % remove outliers
        switch outlierF
            case 0
                bhv_var_tmp = behavior_vars(:,iBhv);
                brainMb_var_tmp = brainMb_vars(:,iBrainMb);
            case 1
                [~, ~, bhv_var_tmp] = rmv_outliers_3sd(behavior_vars(:,iBhv));
                [~, ~, brainMb_var_tmp] = rmv_outliers_3sd(brainMb_vars(:,iBrainMb));
        end
        goodS_tmp = ~isnan(bhv_var_tmp.*brainMb_var_tmp);
        % perform correlations
        corr_nm = [bhv_var_names{iBhv},'_f_',brainMb_var_names{iBrainMb}];
        [r_corr.(corr_nm), pval.(corr_nm)] = corr(behavior_vars(goodS_tmp,iBhv), brainMb_vars(goodS_tmp,iBrainMb));
        % store into general matrix for figure
        corr_mtrx(iBrainMb, iBhv) = r_corr.(corr_nm);
        pval_mtrx(iBrainMb, iBhv) = pval.(corr_nm);
        NS_goodS.(corr_nm) = sum(goodS_tmp);
    end % brain metabolites loop
end % behavioral loop

%% figure
% general parameters
corr_range = [-1, 1];
color_range_choices = redblue(45);
[~, ~, col] = general_fig_prm;

% figure
fig;
subplot_hdl = subplot(1,2,1);
imagesc(corr_mtrx, corr_range);
colormap(subplot_hdl, color_range_choices);
cbar = colorbar;
cbar.Label.String = 'r';
xticks(1:n_bhv_vars);
xticklabels(bhv_var_names_bis);
yticks(1:n_brainMb_vars);
yticklabels(brainMb_var_names_bis);
legend_size(15);
% add stars in the graph if some correlations are significant
for iBhv = 1:n_bhv_vars
    for iBrainMb = 1:n_brainMb_vars
        if pval_mtrx(iBrainMb, iBhv) <= 0.05
            if pval_mtrx(iBrainMb, iBhv) > 0.01 && pval_mtrx(iBrainMb, iBhv) <= 0.05
                pval_hdl = text(iBhv, iBrainMb, '*');
            elseif pval_mtrx(iBrainMb, iBhv) > 0.001 && pval_mtrx(iBrainMb, iBhv) <= 0.01
                pval_hdl = text(iBhv, iBrainMb, '**');
            elseif pval_mtrx(iBrainMb, iBhv) <= 0.001
                pval_hdl = text(iBhv, iBrainMb, '***');
            end % p.value
            % adjust p.value parameters
            pval_hdl.Color = col.black;
            pval_hdl.FontSize = 20;
            pval_hdl.FontWeight = 'bold';
            pval_hdl.HorizontalAlignment = 'center'; % center text on x-axis
            pval_hdl.VerticalAlignment = 'middle'; % center text on y-axis
        end % when p.value is significant
    end % loop over Y variables
end % loop over X variables
