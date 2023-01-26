%% split subjects according to level of brain metabolites

%% extract force data
[study_nm, cond, subject_id, sub_males, sub_females,...
    MVC_all, PCSA_all,...
    MVC_m, PCSA_m,...
    MVC_f, PCSA_f] = grip_MVC_vs_PCSA(0);

%% split the data according to metabolites
[low_met_subs, high_met_subs,...
    metabolite_nm, MRS_ROI_nm, metabolite_allSubs] = medSplit_metabolites(study_nm, subject_id);
[low_met_male_subs, high_met_male_subs,...
    metabolite_m_nm, MRS_ROI_nm, metabolite_male_subs] = medSplit_metabolites(study_nm, sub_males);
[low_met_female_subs, high_met_female_subs,...
    metabolite_f_nm, MRS_ROI_nm, metabolite_female_subs] = medSplit_metabolites(study_nm, sub_females);

%% figure
pSize = 40;
lWidth = 3;
black = [0 0 0];
grey = [143 143 143]./255;
blue = [0 0 1];
green = [0 1 0];
lightBlue = [0 255 255]./255;
lightGreen = [106 204 132]./255;

% all subjects
fig;
scat_low_met_hdl = scatter(PCSA_all(low_met_subs), MVC_all(low_met_subs));
scat_low_met_hdl.LineWidth = lWidth;
scat_low_met_hdl.MarkerEdgeColor = black;
scat_high_met_hdl = scatter(PCSA_all(high_met_subs), MVC_all(high_met_subs));
scat_high_met_hdl.LineWidth = lWidth;
scat_high_met_hdl.MarkerEdgeColor = lightGreen;
legend([scat_high_met_hdl, scat_low_met_hdl],...
    {['high ',metabolite_nm],['low ',metabolite_nm]});
legend('boxoff');
legend('Location','NorthWest');
xlabel('Theoretical maximal force (N)');
ylabel('Calibrated force (N)');
legend_size(pSize);

%% split data by gender
fig;
% females
scat_low_met_f_hdl = scatter(PCSA_f(low_met_female_subs), MVC_f(low_met_female_subs));
scat_low_met_f_hdl.LineWidth = lWidth;
scat_low_met_f_hdl.MarkerEdgeColor = green;
scat_high_met_f_hdl = scatter(PCSA_f(high_met_female_subs), MVC_f(high_met_female_subs));
scat_high_met_f_hdl.LineWidth = lWidth;
scat_high_met_f_hdl.MarkerEdgeColor = lightGreen;
% males
scat_low_met_m_hdl = scatter(PCSA_m(low_met_male_subs), MVC_m(low_met_male_subs));
scat_low_met_m_hdl.LineWidth = lWidth;
scat_low_met_m_hdl.MarkerEdgeColor = blue;
scat_high_met_m_hdl = scatter(PCSA_m(high_met_male_subs), MVC_m(high_met_male_subs));
scat_high_met_m_hdl.LineWidth = lWidth;
scat_high_met_m_hdl.MarkerEdgeColor = lightBlue;

legend([scat_high_met_f_hdl, scat_low_met_f_hdl,...
    scat_high_met_m_hdl, scat_low_met_m_hdl],...
    {['high ',metabolite_f_nm,' - females'],['low ',metabolite_f_nm,' - females'],...
    ['high ',metabolite_m_nm,' - males'],['low ',metabolite_m_nm,' - males']});
legend('boxoff');
legend('Location','NorthWest');
xlabel('Theoretical maximal force (N)');
ylabel('Calibrated force (N)');
legend_size(pSize);