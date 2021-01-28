function[n_mental_max_perTrial, n_mental_max_allCalib] = LGCM_mental_calib(scr, stim)
%[n_mental_max_perTrial] = LGCM_mental_calib(scr, stim)
% LGCM_mental_calib will extract maximum number of subsequent correct
% answers participants can provide in the limited amount of time that is
% available for them to answer.
%
% INPUTS
%
% OUTPUTS
% n_mental_max_perTrial: maximum number of subsequent correct answers per
% trial
%
% n_mental_max_allCalib: maximum number of subsequent correct answers
% provided for each trial during the calibration
%

n_tests = 5;
n_nb_toStart = 10;

for iTest = 1:n_tests
    switch iTest
        case 1 % first trial just keep asking without any reference
            
        otherwise % next trials use initial participant max + 3
            
    end
end % number of tests to try to get max

%% get maximum for the participant
n_mental_max_perTrial = nanmax(n_mental_max_allCalib);

end % function