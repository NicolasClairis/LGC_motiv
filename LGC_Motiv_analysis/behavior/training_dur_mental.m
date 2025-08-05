function[mental_training_dur] = training_dur_mental()
% [mental_training_dur] = training_dur_mental()
% training_dur_mental will check the duration of the mental training in
% terms of number of trials
%
% OUTPUTS
% mental_training_dur: structure containing the number of trials
% (n_learning_trials) it took to finalize the mental learning phase.

%% subject selection
[study_nm, ~, ~, subject_id, NS] = sub_id;

%% working directory
% training_folder = [filesep,filesep,fullfile('sv-nas1.rcp.epfl.ch','human_data_private',...
%     'raw_data_subject',study_nm), filesep];

%% loop through subjects
mental_training_dur.n_learning_trials.allSubs = NaN(1,NS);
for iS = 1:NS
    sub_nm = subject_id{iS};
    % load the data
    if ~strcmp(sub_nm,'017') % training data missing for CID017
        mental_training_dur.n_learning_trials.allSubs(iS) = getfield(load(fullfile(training_folder,['CID',sub_nm],...
            'behavior',['training_data_CID',sub_nm,'.mat']),'jLearningTrial'),'jLearningTrial');
    end
end % subject loop

end % function