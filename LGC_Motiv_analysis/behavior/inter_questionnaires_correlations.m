
%% extract subjects of interest
study_nm = 'study1';
condition = subject_condition;
gender = 'all';
[subject_id, NS] = LGCM_subject_selection(study_nm, condition, gender);
%% load questionnaires
[excelReadQuestionnairesFile, sub_CID_list] = load_questionnaires_data();
quest_names = {'none';'JPI_RScore';...
    'MPSTEFSPhysicalTraitScore';'MPSTEFSMentalTraitScore';...
    'MADRS_SCorrected';'PunishmentScore';'RewardScore';'SHAPScore';...
    'IPAQ';'IPAQInactivity'
    'PRF_DScore';...
    'PSS_14Score';...
    'SIASScore';'STAITraitScore';...
    'CTQEmotionalAbuseScore';'CTQPhysicalAbuseScore';...
    'CTQPhysicalNeglectScore';'CTQEmotionalNeglectScore';...
    'CTQSexualAbuseScore';'CTQMinimizationDenialScore';
    'AI_EP';'AI_I';'IC_I';'IC_N';'IC_M';'IC_S';'ER';'SA';...
    'Non_PlanningImpulsivenessScore';'MotorImpulsivenessScore';'AttentionalImpulsivenessScore';...
    'EnjoymentOfCompetition';'Contentiousness';'Honesty_Humility';...
    'Emotionality';'Extraversion';'Agreeableness';...
    'Conscientiousness';'OpennessToExperience'};

%% load general variables
[excelReadGeneralFile] = load_gal_data_bis(study_nm);
gal_data_nm = {'none';'BMI';...
    'HeuresDeSommeilLaVeilleDeL_exp_rience';...
    'HeuresDeSommeil_enMoyenne_';...
    'SocialComparisonLadder'};

%% ask what to correlate with what?
nVars = 2;
var_names = cell(1,nVars);
for iVar = 1:nVars
    
    quest_idx = listdlg('PromptString',['Variable of interest ',num2str(iVar),'?'],...
        'ListString',quest_names,...
        'SelectionMode','single');
    if ~strcmp(quest_names{quest_idx},'none')
        var_names{iVar} = quest_names{quest_idx};
        switch iVar
            case 1
                var1_tmp = excelReadQuestionnairesFile.(quest_names{quest_idx});
                var1_sub_tmp = convert_sub_id_from_num_to_cell(excelReadQuestionnairesFile.CID');
            case 2
                var2_tmp = excelReadQuestionnairesFile.(quest_names{quest_idx});
                var2_sub_tmp = convert_sub_id_from_num_to_cell(excelReadQuestionnairesFile.CID');
        end
    else
        gal_idx = listdlg('PromptString',['Variable of interest ',num2str(iVar),'?'],...
            'ListString',gal_data_nm,...
            'SelectionMode','single');
        var_names{iVar} = gal_data_nm{gal_idx};
        switch iVar
            case 1
                var1_tmp = excelReadGeneralFile.(gal_data_nm{gal_idx});
                var1_sub_tmp = strrep(excelReadGeneralFile.CID,'CID','');
            case 2
                var2_tmp = excelReadGeneralFile.(gal_data_nm{gal_idx});
                var2_sub_tmp = strrep(excelReadGeneralFile.CID,'CID','');
        end
    end
end

%% extract variables for subjects of interest
[var1, var2] = deal(NaN(1,NS));
for iS = 1:NS
    sub_nm = subject_id{iS};
    sub_idx_var1 = strcmp(var1_sub_tmp, sub_nm);
    sub_idx_var2 = strcmp(var2_sub_tmp, sub_nm);
    if (sum(sub_idx_var1) == 1) && (sum(sub_idx_var2) == 1)
        var1(iS) = var1_tmp(sub_idx_var1);
        var2(iS) = var2_tmp(sub_idx_var2);
    end
end

%% perform correlation
goodSubs = ~isnan(var1.*var2);
[beta,~,stats] = glmfit(var1(goodSubs), var2(goodSubs), 'normal');
var1_sorted = sort(var1(goodSubs));
var2_fit = glmval(beta, var1_sorted,'identity');
pval = stats.p;
disp(['p = ',num2str(pval(2))]);

%% draw corresponding figure
pSize = 30;
lWidth = 3;
grey = [143 143 143]./255;

fig;
hold on;
scat_hdl = scatter(var1(goodSubs), var2(goodSubs));
scat_hdl.MarkerEdgeColor = 'k';
scat_hdl.LineWidth = lWidth;
fit_hdl = plot(var1_sorted, var2_fit);
fit_hdl.LineWidth = lWidth;
fit_hdl.Color = grey;
fit_hdl.LineStyle = '--';
xlabel(var_names{1});
ylabel(var_names{2});
legend_size(pSize);