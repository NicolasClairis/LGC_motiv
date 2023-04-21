%% check where the behavioral outliers are located in terms of metabolite distribution

%% load all subjects
study_nm = 'study1';
condition = 'behavior';
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% load saturators
figGrpDisp = 0;
[choiceND_perRun, saturationSubs] = choiceNDproportion_perRun_group(figGrpDisp);
% extract non-outliers distribution
okSubs = ~ismember(subject_id, saturationSubs.fullNonDef.subList);
goodSubs = subject_id(okSubs);
% extract outliers and type of outlier
badSubs_names = saturationSubs.fullNonDef.subList;
EpSatSubs = false(1,NS);
EmSatSubs = false(1,NS);
EpANDEmSatSubs = false(1,NS);
for iBadS = 1:length(badSubs_names)
    badSub_nm = badSubs_names{iBadS};
    jSub_idx = strcmp(badSub_nm, subject_id);
    [runs, n_runs] = runs_definition(study_nm, badSub_nm, condition);
    satRunNames_tmp = fieldnames(saturationSubs.satRunsPerSub.(['CID',badSub_nm]));
    for iSatRun = 1:length(satRunNames_tmp)
        jSatRun = str2double(satRunNames_tmp{iSatRun}(end));
        if strcmp(runs.tasks{jSatRun},'Ep') && ~strcmp(badSub_nm, EpSatSubs)
            EpSatSubs(jSub_idx) = true;
        elseif strcmp(runs.tasks{jSatRun},'Em') && ~strcmp(badSub_nm, EmSatSubs)
            EmSatSubs(jSub_idx) = true;
        end
    end % saturation runs loop
    
    if strcmp(badSub_nm, EpSatSubs) && strcmp(badSub_nm, EmSatSubs)
        EpANDEmSatSubs(jSub_idx) = true;
    end
end % bad subjects loop

%% load all metabolites
[metabolites] = metabolite_load(subject_id);
MRS_ROIs = fieldnames(metabolites);
n_MRS_ROIs = length(MRS_ROIs);

%% display distribution
normalCol = [0 0 255]./255;
EpSatCol = [0 255 0]./255;
EmSatCol = [143 55 0]./255;
EpANDEmSatCol = [255 143 0]./255;
pSize = 50;
bWidth = 0.1;
nBins = 6;
for iMRS_ROI = 1:n_MRS_ROIs
    MRS_ROI_nm = MRS_ROIs{iMRS_ROI};
    metabolite_names = fieldnames(metabolites.(MRS_ROI_nm));
    n_metabolites = length(metabolite_names);
    for iM = 1:n_metabolites
        metab_nm = metabolite_names{iM};
        
        % extract distribution on everybody
        [okSubsCounts, okSubsEdges] = histcounts(metabolites.(MRS_ROI_nm).(metab_nm)(okSubs), nBins);
        [EpSatSubsCounts] = histcounts(metabolites.(MRS_ROI_nm).(metab_nm)(EpSatSubs), okSubsEdges);
        [EmSatSubsCounts] = histcounts(metabolites.(MRS_ROI_nm).(metab_nm)(EmSatSubs), okSubsEdges);
        [EpANDEmSubsCounts] = histcounts(metabolites.(MRS_ROI_nm).(metab_nm)(EpANDEmSatSubs), okSubsEdges);
        % average
        okSubsX = NaN(1,nBins);
        for iBin = 1:nBins
            okSubsX(iBin) = mean([okSubsEdges(iBin), okSubsEdges(iBin+1)]);
        end
        
        %% display graph
        fig;
        hold on;
        bar_hdl = bar(okSubsX,...
            [okSubsCounts; EpSatSubsCounts; EmSatSubsCounts; EpANDEmSubsCounts],...
            'stacked');
        
        % fix colours
        bar_hdl(1).FaceColor = normalCol;
        bar_hdl(2).FaceColor = EpSatCol;
        bar_hdl(3).FaceColor = EmSatCol;
        bar_hdl(4).FaceColor = EpANDEmSatCol;
        xlabel([MRS_ROI_nm,' - ',metab_nm]);
        ylabel('Number');
        legend([bar_hdl(1), bar_hdl(2), bar_hdl(3), bar_hdl(4)],...
            {'normal','Ep','Em','Ep+Em'});
        legend('Location','Northeast');
        legend('boxoff');
        legend_size(pSize);
    end % metabolite loop
end % MRS ROI loop