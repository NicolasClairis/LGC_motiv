function [numberVector_calib] = mental_calibNumberVector(n_calibTrials, n_calibMax)
%[numberVector_calib] = mental_calibNumberVector(n_calibTrials, n_calibMax)
%mental_calibNumberVector defines the vector of numbers to use for all
%participants. Vectors have been defined to be as much as possible hard and
%diverse.
%
% INPUTS
% n_calibTrials: number of calibration trials
%
% n_calibMax: number of correct answers to reach
%
% OUTPUTS
% numberVector_calib: matrix with 1 vector of numbers per line
% corresponding to each calibration trial.

numberVector_calib = NaN(n_calibTrials, 8*4);

if n_calibMax <= 6 % need to have a few more numbers in case of mistakes during calibration (but no need to have too many since there is a time limit
    
    switch n_calibTrials
        case 1
            numberVector_calib(1,:) = [1 7 6 4 9 2 3 8, 4 2 8 9 6 3 1 7, 1 7 2 9 8 3 4 6, 1 9 2 6 3 8 4 7];
        case 2
            numberVector_calib(1,:) = [1 7 6 4 9 2 3 8, 4 2 8 9 6 3 1 7, 1 7 2 9 8 3 4 6, 1 9 2 6 3 8 4 7];
            numberVector_calib(2,:) = [9 8 3 7 1 2 6 4, 3 2 7 1 8 9 6 4, 6 7 2 4 1 9 8 3, 2 7 1 8 6 9 3 4];
        case 3
            numberVector_calib(1,:) = [1 7 6 4 9 2 3 8, 4 2 8 9 6 3 1 7, 1 7 2 9 8 3 4 6, 1 9 2 6 3 8 4 7];
            numberVector_calib(2,:) = [9 8 3 7 1 2 6 4, 3 2 7 1 8 9 6 4, 6 7 2 4 1 9 8 3, 2 7 1 8 6 9 3 4];
            numberVector_calib(3,:) = [3 6 1 4 7 9 2 8, 1 3 4 6 9 7 8 2, 9 4 8 7 3 2 6 1, 6 4 3 7 2 1 9 8];
        case 5
            numberVector_calib(1,:) = [1 7 6 4 9 2 3 8, 4 2 8 9 6 3 1 7, 1 7 2 9 8 3 4 6, 1 9 2 6 3 8 4 7];
            numberVector_calib(2,:) = [9 8 3 7 1 2 6 4, 3 2 7 1 8 9 6 4, 6 7 2 4 1 9 8 3, 2 7 1 8 6 9 3 4];
            numberVector_calib(3,:) = [3 6 1 4 7 9 2 8, 1 3 4 6 9 7 8 2, 9 4 8 7 3 2 6 1, 6 4 3 7 2 1 9 8];
            numberVector_calib(4,:) = [2 9 1 8 6 4 7 3, 7 2 6 4 3 9 8 1, 7 9 8 4 3 2 1 6, 8 4 7 6 2 9 1 3];
            numberVector_calib(5,:) = [7 2 8 3 9 6 1 4, 6 2 3 7 4 9 8 1, 6 1 4 3 9 8 7 2, 9 8 1 3 6 2 7 4];
        otherwise
            error('case not ready yet');
    end
    
else
    error(['case where n_calibMax = ',num2str(n_calibMax),' not ready yet.']);
end

end % function