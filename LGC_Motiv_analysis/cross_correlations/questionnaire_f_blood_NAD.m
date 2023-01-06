% test whether whole-blood NAD metabolites can (partially) predict the
% scores in the motivational questionnaires related to depression, energy,
% etc.

%% by default display the figures but can be set to 0 if you only want the p.value
figDisp = 1;

%% define subjects
study_nm = 'study1';
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% working directories
root = LGCM_root_paths;
% study
study_nm = 'study1';
studyPath = [root, filesep, study_nm, filesep];
%% extract all blood data
[bloodTable, blood_NAD_sub_List] = load_blood_NAD(study_nm);
bloodMb_names = fieldnames(bloodTable);
n_Mb = length(bloodMb_names);
for iMb = 1:n_Mb
    bloodMb.(bloodMb_names{iMb}) = NaN(1,NS);
end % metabolite loop

%% load questionnaire
potentialQuestionnaires = {'JPI_RScore', 'MADRS_SCorrected',...
    'MPSTEFSPhysicalTraitScore','MPSTEFSMentalTraitScore',...
    'PunishmentScore','RewardScore',...
    'IPAQ','IPAQInactivity'};
selectedCon = listdlg('PromptString','Which questionnaires?',...
    'ListString',potentialQuestionnaires);
nQuest = length(selectedCon);
questToCheck = potentialQuestionnaires(selectedCon);
[excelReadQuestionnairesFile, sub_CID_list] = load_questionnaires_data();
for iQuest = 1:nQuest
    quest_nm = questToCheck{iQuest};
    quest_data.(quest_nm) = NaN(1,NS);
end

%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    sub_nm_bis = ['CID',sub_nm];
    
    %% load questionnaire
    sub_quest_idx = strcmp(sub_CID_list, sub_nm);
    for iQuest = 1:nQuest
        quest_nm = questToCheck{iQuest};
        quest_data.(quest_nm)(iS) = excelReadQuestionnairesFile.(quest_nm)(sub_quest_idx);
    end % questionnaire loop
    %% load blood
    sub_blood_idx = strcmp(sub_nm_bis, blood_NAD_sub_List);
    % extract blood info
    for iMb = 1:n_Mb
        bloodMb_nm = bloodMb_names{iMb};
        bloodMb.(bloodMb_nm)(iS) = bloodTable.(bloodMb_nm)(sub_blood_idx);
    end
end % subject loop

%% test correlations
for iMb = 1:n_Mb
    bloodMb_nm = bloodMb_names{iMb};
    for iQuest = 1:nQuest
        quest_nm = questToCheck{iQuest};
        corr_nm = [quest_nm,'_f_',bloodMb_nm];
        goodSubs.(corr_nm) = (~isnan(quest_data.(quest_nm))).*(~isnan(bloodMb.(bloodMb_nm))) == 1;
        [beta.(corr_nm),~,stats_tmp] =...
            glmfit(bloodMb.(bloodMb_nm)(goodSubs.(corr_nm)),...
            quest_data.(quest_nm)(goodSubs.(corr_nm)), 'normal');
        pval.(corr_nm) = stats_tmp.p;
        % store significant p.values for slope
        if stats_tmp.p(2) < 0.05
            pval.signif.(corr_nm) = stats_tmp.p(2);
        elseif stats_tmp.p(2) > 0.05 && stats_tmp.p(2) < 0.1
            pval.almostSignif.(corr_nm) = stats_tmp.p(2);
        end
        bloodMb_sort.(corr_nm) = sort(bloodMb.(bloodMb_nm)(goodSubs.(corr_nm)));
        quest_fit.(corr_nm) = glmval(beta.(corr_nm), bloodMb_sort.(corr_nm), 'identity');
    end % questionnaire
end % metabolite

%% correlation and figure
if figDisp == 1
    lWidth = 3;
    pSize = 25;
    black = 'k';
    orange = [254 75 3]./255;
    %% show results
    for iQuest = 1:nQuest
        quest_nm = questToCheck{iQuest};
        fig1 = fig; j_fig1 = 0;
        fig2 = fig; j_fig2 = 0;
        for iMb = 1:n_Mb
            bloodMb_nm = bloodMb_names{iMb};
            corr_nm = [quest_nm,'_f_',bloodMb_nm];
            % define y.label for questionnaires
            switch quest_nm
                case 'JPI_RScore'
                    quest_nm_bis = 'JPI-R';
                case 'MADRS_SCorrected'
                    quest_nm_bis = 'MADRS-S';
                case 'MPSTEFSPhysicalTraitScore'
                    quest_nm_bis = 'MPSTEFS physical';
                case 'MPSTEFSMentalTraitScore'
                    quest_nm_bis = 'MPSTEFS mental';
                case 'PunishmentScore'
                    quest_nm_bis = 'PANAS P';
                case 'RewardScore'
                    quest_nm_bis = 'PANAS R';
                case 'IPAQ'
                    quest_nm_bis = 'IPAQ activity';
                case 'IPAQInactivity'
                    quest_nm_bis = 'IPAQ inactivity';
            end
            
            switch bloodMb_nm
                case {'Nam','NMN','NR','NAD',...
                        'NADH','NADP','NADPH','MeNam',...
                        'MeXPY'}
                    figure(fig1);
                    j_fig1 = j_fig1 + 1;
                    subplot(3,4,j_fig1);
                case {'NAD_div_NADH',...
                        'NADP_div_NADPH',...
                        'total_NAD_precursors',...
                        'total_NAD',...
                        'total_NAD_with_precursors',...
                        'total_NAD_with_byproducts',...
                        'total_NAD_byproducts'}
                    figure(fig2);
                    j_fig2 = j_fig2 + 1;
                    subplot(3,3,j_fig2);
            end
            %% figure
            hold on;
            scat_hdl = scatter(bloodMb.(bloodMb_nm)(goodSubs.(corr_nm)),...
                quest_data.(quest_nm)(goodSubs.(corr_nm)));
            plot_hdl = plot(bloodMb_sort.(corr_nm), quest_fit.(corr_nm));
            scat_hdl.LineWidth = lWidth;
            scat_hdl.MarkerEdgeColor = black;
            plot_hdl.Color = orange;
            plot_hdl.LineStyle = '--';
            [blood_labelname] = blood_label(bloodMb_nm);
            xlabel(blood_labelname);
            ylabel(quest_nm_bis);
            legend_size(pSize);
        end % metabolite loop
    end % questionnaire loop
end % fig disp