function[avgSleep, prevDaySleep, delta_PrevDay_AvgSleep] = extract_sleep(study_nm, subject_id, NS)
% [avgSleep, prevDaySleep, delta_PrevDay_AvgSleep] = extract_sleep(study_nm, subject_id, NS)
% extract_sleep will extract the average amount of sleep hours, the
% previous day amount of sleep hours and the delta between previous day and
% average sleep hours (all expressed in minutes as unity)
%
% INPUTS
% study_nm: study name
%
% subject_id: list of subjects
%
% NS: number of subjects
%
% OUTPUTS
% avgSleep: average sleep hours (in minutes)
%
% prevDaySleep: previous day sleep (in minutes)
%
% delta_PrevDay_AvgSleep: difference between previous day and average
% amount of sleep (in minutes)

%% load sleep
[excelReadGeneralFile] = load_gal_data_bis(study_nm);
prevDaySleepTable = excelReadGeneralFile.HeuresDeSommeilLaVeilleDeL_exp_rience;
avgSleepTable = excelReadGeneralFile.HeuresDeSommeil_enMoyenne_;
sleepCID = excelReadGeneralFile.CID;

% replace 'h' of hours by ':'
prevDaySleepTable = strrep(prevDaySleepTable, 'h',':');
avgSleepTable = strrep(avgSleepTable, 'h',':');

% extract the data
[avgSleep, prevDaySleep] = deal(NaN(1,NS));

%% loop through subjects to extract nutrition + sleep
for iS = 1:NS
    sub_nm = subject_id{iS};
    % load sleep
    % loop through subjects to find the correct one
    for iLine = 1:size(sleepCID,1)
        isItThisSubject = strcmp(sub_nm, sleepCID{iLine}(4:6));
        if isItThisSubject == true
            sub_sleep_idx = iLine;
        end
    end
    avgSleep(iS) = minutes(duration(avgSleepTable{sub_sleep_idx},'InputFormat','hh:mm'));
    prevDaySleep(iS) = minutes(duration(prevDaySleepTable{sub_sleep_idx},'InputFormat','hh:mm'));
end % subject loop

%% look also at the difference between previous day and average
delta_PrevDay_AvgSleep = prevDaySleep - avgSleep;

end % function