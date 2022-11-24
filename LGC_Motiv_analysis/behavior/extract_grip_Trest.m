function[Trest1, Trest2] = extract_grip_Trest(subBehaviorFolder, sub_nm, run_nm)
% [Trest1, Trest2] = extract_grip_Trest(subBehaviorFolder, sub_nm, run_nm)
% extract_grip_Trest will extract several parameters related to the time
% spent "resting" (ie not squeezing the handgrip).
%
% INPUTS
% subBehaviorFolder: folder where data is stored
%
% sub_nm: string with subject CID name
%
% run_nm: string with run name
%
% OUTPUTS
% Trest1: time resting defined as feedback + pre-choice cross jitter 
% duration between two efforts
%
% Trest2: similar to Trest1, but when the easy effort was made in the
% previous trial, it will also include the whole previous trial duration
% as time resting
%
% N. Clairis - november 2022

%% general parameters
nTrialsPerRun = 54;
[Trest1,...
    Trest2] = deal(NaN(1,nTrialsPerRun));
%% load the data
behaviorStruct = load([subBehaviorFolder,...
    'CID',sub_nm,'_session',run_nm,'_physical_task.mat']);
durations = behaviorStruct.physicalPerf.durations;
% extract choice information
[choice_highE] = extract_choice_hE(subBehaviorFolder, sub_nm, run_nm, 'physical');

for iTrial = 1:nTrialsPerRun
    if iTrial > 1 % no point to rest before trial 1
        prevTrial = iTrial - 1;
        %% time resting = time between last effort and choice
        Trest1(iTrial) = durations.fbk(prevTrial) +...
            durations.preChoiceCross(iTrial);

        %% time resting = time between last non-default effort and choice
        % (ie ignore low effort choices)

        % was the previous choice a low or a high effort choice
        prevChoice_hE = choice_highE(prevTrial) == 1;
        %% adapt Trest2 according to the choice made before
        if prevChoice_hE == 1 % high effort chosen previously
            Trest2(iTrial) = durations.fbk(prevTrial) +...
                durations.preChoiceCross(iTrial);
        elseif prevChoice_hE == 0 % low effort chosen previously
            %% trial 2 = only 1 trial before
            if iTrial == 2
                Trest2(iTrial) = durations.preChoiceCross(prevTrial) +...
                    durations.dispChoiceOptions(prevTrial) +...
                    durations.dispChoice(prevTrial) +...
                    durations.preEffortCross(prevTrial) +...
                    durations.effortPeriod(prevTrial) +...
                    durations.fbk(prevTrial) +...
                    durations.preChoiceCross(iTrial);
            else % add feedback time of previous-previous trial
                prevPrevTrial = iTrial - 2;
                Trest2(iTrial) = durations.fbk(prevPrevTrial) +...
                    durations.preChoiceCross(prevTrial) +...
                    durations.dispChoiceOptions(prevTrial) +...
                    durations.dispChoice(prevTrial) +...
                    durations.preEffortCross(prevTrial) +...
                    durations.effortPeriod(prevTrial) +...
                    durations.fbk(prevTrial) +...
                    durations.preChoiceCross(iTrial);
            end
        end
    end % trial filter
end % trial loop


end % function