function[training_dur] = training_duration()
% [training_dur] = training_duration()
% training_duration will extract the average duration of the training
% period.
%
% OUTPUTS
% training_dur: structure will include mean, SEM, SD

%% subject selection
[study_nm, condition, gender, subject_id, NS] = sub_id;

%% working directory
training_folder = [fullfile('M:','human_data_private',...
    'raw_data_subject',study_nm), filesep];

%% load the data
training_dur_table = readtable([training_folder,filesep,'General1.xlsx'],...
    'Sheet','training duration');

%% loop through subjects
training_dur.allSubs = NaN(1,NS);
for iS = 1:NS
    sub_nm = subject_id{iS};
    sub_idx = strcmp(training_dur_table.IDParticipant, ['CID',sub_nm]);
    if ~isempty(sub_idx) && size(sub_idx,2) == 1
        training_dur.allSubs(iS) = training_dur_table.timeConvertedInMinutes(sub_idx);
    elseif size(sub_idx,2) > 1
        error(['Subject ',sub_nm,' has more than one presence.']);
    end
end % subject loop

%% stats
[training_dur.mean,...
    training_dur.sem,...
    training_dur.sd] = mean_sem_sd(training_dur.allSubs, 2);

end % function