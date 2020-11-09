
% script for convolution of pupil data with all the onsets

% subject selection

% working directories
root.root = ['B:',filesep,'resultats',filesep];

% main parameters
nRuns_per_Task = 3; % 3 runs per task
SR = 1000; % sampling rate: 1kHz for EyeLink (=1000 datapoints/sec = 1 datapoint/ms)

% loop through studies (eventually)
for iStudy = 1:2
    if iStudy == 1
        path = 'MBB_june2016';
    elseif iStudy == 2
        path = 'multiseq_april_may2016';
    end
    
    root.subjects = [root.root, path, filesep];
    
    
    % loop through subjects
    for iSubject = 1:NS
        subName = subject_id{iSubject};
        root.subject = [rootSubjests, subName, filesep];
        root.subject.behavior = [root_subject, 'behavior',filesep];
        root.subject.preprocPupil = [root_subject, 'eye_analysis', filesep, 'preprocessed_files', filesep];
        root.subject.onsetsfMRIFolder = [root_subject, 'fMRI_analysis',filesep];
        
        
        % loop through tasks
        for iTask = 1:3
            
            % loop through runs
            for iRun = (1:3) + iTask*nRuns_per_Task
                runName = num2str(iRun);
                
                % load pupil data
                load([root_subject_preprocPupil, 'eye_preproc_',subName,'_',taskAbrev,'_run',runName,'.mat'],...
                    'pupil_diam','z_pupil_diam','nanEpisodes','eyeTime');
                
                % remove any data point before T0 (= first fMRI trigger)
                pupil_diam      = pupil_diam(eyeTime >= 0);
                z_pupil_diam    = z_pupil_diam(eyeTime >= 0);
                eyeTime = eyeTime(eyeTime >= 0);
                nanEpisodes.blink_episodes = nanEpisodes.blink_episodes(eyeTime >= 0);
                
                % load the onsets and the duration
                load([root_subject_behavior,'MBB_battery_',task,'_onsets_sub',subid,'_sess',runName,'.mat'],...
                    'onset','duration','rt');
                
                % normalize the onsets by T0 (should be aligned with the
                % pupil data) + convert from seconds to ms
                load([root_subject_behavior,'TTL_sub',subid,'_sess',runName,'.mat'],'T0');
                onset_fieldnames = fieldnames(onset);
                for iField = 1:length(onset_fieldnames)
                    onset.(onset_fieldnames{iField})    = (onset.(onset_fieldnames{iField}) - T0)*SR;
                    duration.(onset_fieldnames{iField}) = duration.(onset_fieldnames{iField})*SR;
                end
                
                %% delete too slow/fast samples
                if iTask == 1
                elseif iTask == 2
                    
                elseif iTask == 3
                    fieldsToFilter = {'display_optionsRE','ChoiceRE'};
                    load([root.subject.onsetsfMRIFolder,'onsets_RE_sub',subid,'_sess',runName,'.mat'],'RE');
                    trialN = RE.trialN_RE; % extract trial index filtered
                end
                for iField = 1:length(fieldsToFilter)
                    for iTrial = length(onset.(fieldsToFilter{iField})):-1:1
                        if ismember(iTrial,trialN) == 0
                            onset.(fieldsToFilter{iField})(iTrial) = [];
                            duration.(fieldsToFilter{iField})(iTrial) = [];
                        end
                    end
                end
                
                %% delete trials mainly based on interpolation but where
                % signal not acquired
                for iField = 1:length(onset_fieldnames)
                    eventName = onset_fieldnames{iField};
                    if strcmp(eventName,'firstTTL') == 0 &&... % ignore TTL
                            strcmp(eventName,'cross_ITT') == 0 % ignore intertask cross for RE task
                        current_onset       = onset.(eventName);
                        current_event_dur   = duration.(eventName); % duration of the corresponding event (in seconds)
                        current_endset      = current_onset + current_event_dur;
                        if ~isempty(current_onset)
                            
                            for iSample = length(current_onset):-1:1
                                % estimate percentage of signal for this
                                % sample
                                sampleNaNvalues = nanEpisodes.blink_episodes(floor(current_onset(iSample):current_endset(iSample)));
                                nanPercentage = (sum(sampleNaNvalues == 0)/length(sampleNaNvalues))*100;
                                % delete samples where more than [X]% of
                                % NaN values
                                if nanPercentage > 70
                                    onset.(eventName)(iSample)      = [];
                                    duration.(eventName)(iSample)   = [];
                                end
                            end
                        end
                    end
                end
                
                %% extract convolution for all events
                for iField = 1:length(onset_fieldnames)
                    eventName = onset_fieldnames{iField};
                    if strcmp(eventName,'firstTTL') == 0 &&... % ignore TTL
                            strcmp(eventName,'cross_ITT') == 0 &&...% ignore intertask cross for RE task
                            strcmp(eventName,'cross_ITI') == 0    % temporarily ignore ITI until finding a way to include it properly
                        current_onset       = onset.(eventName);
                        current_event_dur   = duration.(eventName); % duration of the corresponding event (in seconds)
                        current_endset      = current_onset + current_event_dur;
                        
                        %                         if strcmp(eventName,'cross_ITI') || strcmp(eventName,'cross_ITT')
                        %                             flipGamma1 = 0; % expect dilation
                        %                         else
                        flipGamma1 = 1; % expect constriction
                        %                         end
                        
                        if ~isempty(current_onset)
                            
                            for iSample = 1:length(current_onset)
                                timeScale = floor(current_onset(iSample):current_endset(iSample));
                                Y_pupil = z_pupil_diam(timeScale);
                                [PRF_kernel.(eventName).(['sample',num2str(iSample)]),...
                                    PRF_kernel_part1.(eventName).(['sample',num2str(iSample)]),...
                                    PRF_kernel_part2.(eventName).(['sample',num2str(iSample)])] = esti_PRF(Y_pupil, 0, flipGamma1);
                            end
                        end
                    end
                end
                
            end
            
        end
    end
end