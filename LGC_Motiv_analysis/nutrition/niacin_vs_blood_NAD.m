% script to correlate niacin nutritional intake derived from the FFQ to
% blood measurements of the NAD metabolome.
%
%

%% working directory
pcRoot = LGCM_root_paths;
switch pcRoot
    case 'E:\'
        dataStoragePath = 'M:\human_data_private\analyzed_data\study1';
        gitPath = fullfile('C:','Users','clairis','Desktop','GitHub',...
            'LGC_motiv','LGC_Motiv_analysis');
    case 'L:\human_data_private\raw_data_subject\'
        dataStoragePath = 'L:\human_data_private\analyzed_data\study1';
        gitPath = fullfile('C:','Users','Loco','Documents','GitHub',...
            'LGC_motiv','LGC_Motiv_analysis');
end

%% select study and participants
study_nm = 'study1';
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm,condition);

%% extract all blood data
bloodFilePath = [dataStoragePath,filesep,'blood',filesep,...
    '20220915_preliminary_blood_NAD_report.xlsx'];
NAD_blood_table = readtable(bloodFilePath,...
        'Sheet','Final');
%% extract all niacin data
niacineFilePath = [gitPath, filesep, 'nutrition',filesep,'niacine_scoring.xlsx'];
niacin_table = readtable(niacineFilePath,...
        'Sheet','Sheet1');

%% loop through subjects to perform the correlation
blood_metabolites = {'Nam','MeXPY','NAD','NADP','NMN','MeNam','NR',...
    'NADH','NADPH','NAD_div_NADH','NADP_div_NADPH'};
nNADHmetab = length(blood_metabolites);
niacin = NaN(1,NS);
for iBloodM = 1:nNADHmetab
    blood.(blood_metabolites{iBloodM}) = NaN(1,NS);
end
for iS = 1:NS
    sub_nm = subject_id{iS};
    sub_nm_bis = ['CID',sub_nm];
    sub_blood_idx = find(strcmp(NAD_blood_table.ID, sub_nm_bis));
    sub_niacin_idx = find(strcmp(niacin_table.CID, sub_nm));
    if ~isempty(sub_niacin_idx) && size(sub_niacin_idx,1) == 1
        niacin(iS) = niacin_table.NiacineParSemaine_ug_(sub_niacin_idx);
    end
    for iBloodM = 1:nNADHmetab
        blood_m_nm = blood_metabolites{iBloodM};
        if ~isempty(sub_blood_idx)
            blood.(blood_m_nm)(iS) = NAD_blood_table.(blood_m_nm)(sub_blood_idx);
        end
    end
end % subject loop

%% test the correlations
for iBloodM = 1:nNADHmetab
        blood_m_nm = blood_metabolites{iBloodM};
        % check good subjects for each measurement
        goodSubs = (~isnan(blood.(blood_m_nm)).*(~isnan(niacin))) == 1;
        niacin_f_blood_nm = ['niacin_f_',blood_m_nm];
        [b_tmp,~,stats_tmp] = glmfit(...
            niacin(goodSubs), blood.(blood_m_nm)(goodSubs),'normal');
        corr_vals.b.(niacin_f_blood_nm) = b_tmp(2);
        corr_vals.p.(niacin_f_blood_nm) = stats_tmp.p(2);
        fit_blood.(blood_m_nm) = NaN(1,NS);
        fit_blood.(blood_m_nm)(goodSubs) = glmval(b_tmp, niacin(goodSubs), 'identity');
end
%% display the correlation matrix
lWidth = 3;
pSize = 30;

% display each correlation
fig;
for iBloodM = 1:nNADHmetab
    blood_m_nm = blood_metabolites{iBloodM};
    % check good subjects for each measurement
    goodSubs = (~isnan(blood.(blood_m_nm)).*(~isnan(niacin))) == 1;
        
    subplot(3,4,iBloodM);
    hold on;
    hdl = scatter(niacin, blood.(blood_m_nm));
    [niacin_sorted, niacin_sort_idx] = sort(niacin(goodSubs));
    fit_blood_tmp = fit_blood.(blood_m_nm)(goodSubs);
    fit_blood_sorted = fit_blood_tmp(niacin_sort_idx);
    fit_hdl = plot(niacin_sorted, fit_blood_sorted);
    xlabel('niacin (μg/week)');
    switch blood_m_nm
        case 'NAD_div_NADH'
            ylabel('NAD/NADH (μM)');
        case 'NADP_div_NADPH'
            ylabel('NADP/NADPH (μM)');
        otherwise
            ylabel([blood_m_nm,' (μM)']);
    end
    % esthetic improvements
    hdl.LineWidth = lWidth;
    hdl.MarkerEdgeColor = [0 0 0];
    fit_hdl.LineWidth = lWidth;
    fit_hdl.Color = [143 143 143]./255;
    fit_hdl.LineStyle = '--';
    legend_size(pSize);
end % blood metabolite