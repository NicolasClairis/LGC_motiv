function[] = metabolites_f_sleep()
% check correlation between level of metabolites in a given brain area and
% sleep previous day of the experiment.

%% define subjects to check
condition = subject_condition;
study_nm = 'study1';
[subject_id, NS] = LGCM_subject_selection(study_nm,condition,'all');

%% define metabolite of interest
switch study_nm
    case 'study1'
        %% which ROI?
        ROIs = {'dmPFC','aIns'};
        nROIs = length(ROIs);
        ROI_idx = spm_input('Metabolites in which brain area?',1,'m',...
            ROIs,1:nROIs,0);
        ROI_nm = ROIs{ROI_idx};
        %% select metabolite of interest
        metabolites = {'Mac','Ala','Asp','PCho','Cr','PCr','GABA',...
            'Gln','Glu','GSH','Gly','Ins','Lac','NAA','Scyllo','Tau',...
            'Asc','Glc','NAAG','GPC','PE','Ser',...
            'NAA_NAAG','Glu_Gln','GPC_PCho','Cr_PCr','Gly_Ins','Gln_div_Glu'};
    otherwise
        error(['not ready yet for ',study_nm]);
end
n_met = length(metabolites);
metabolite_idx = spm_input('Which metabolite to focus on?',1,'m',...
    metabolites,1:n_met,0);
metabolite_nm = metabolites{metabolite_idx};

%% load sleep
[excelReadGeneralFile] = load_gal_data();
prevDaySleepTable = excelReadGeneralFile.HeureDeSommeilLaVeilleDeL_exp_rience;
avgSleepTable = excelReadGeneralFile.HeuresDeSommeil_enMoyenne_;
sleepCID = excelReadGeneralFile.IDParticipant;
% replace 'h' of hours by ':'
prevDaySleepTable = strrep(prevDaySleepTable, 'h',':');
avgSleepTable = strrep(avgSleepTable, 'h',':');
% extract the data
[avgSleep, prevDaySleep, metabolite] = deal(NaN(1,NS));
for iS = 1:NS
    sub_nm = subject_id{iS};

    %% load sleep
    % loop through subjects to find the correct one
    for iLine = 1:size(sleepCID,1)
        isItThisSubject = strcmp(sub_nm, sleepCID{iLine}(4:6));
        if isItThisSubject == true
            sub_idx = iLine;
        end
    end
    if ~strcmp(avgSleepTable{sub_idx},'NaN')
        avgSleep(iS) = datenum(avgSleepTable{sub_idx},'hh:mm');
    end
    if ~strcmp(prevDaySleepTable{sub_idx},'NaN')
        prevDaySleep(iS) = datenum(prevDaySleepTable{sub_idx},'hh:mm');
    end

    %% load metabolite
    [metabolites] = metabolite_load({sub_nm});
    % focus on metabolite and brain area selected
    metabolite(iS) = metabolites.(ROI_nm).(metabolite_nm);
end % subject loop

% look also at the difference between previous day and average
delta_PrevDay_AvgSleep = prevDaySleep - avgSleep;

%% perform correlations

end % function