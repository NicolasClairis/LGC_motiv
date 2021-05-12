function [F_now_Voltage, timeSample, sampleOk] = F_read(dq, t_readWait)
%[F_now_Voltage, timeSample, sampleOk] = F_read(dq)
% F_read reads the current force level.
%
% INPUTS
% dq: NI device identifier for handgrip input
%
% t_readWait: time for pause
%
% OUTPUTS
% F_now_Voltage: force level in Voltage
%
% timeSample: time when sample was obtained
%
% sampleOk: was the sample read (1) or not?

%% define baseline correction for force
F_baseline = -0.11;

%% read force level
timeSample = GetSecs;
F_now_Voltage_tmp = read(dq,'all','OutputFormat','Matrix');

%% add a pause after read to ensure read works properly
pause(t_readWait);

%% extract force level
if ~isempty(F_now_Voltage_tmp)
    F_now_Voltage = F_now_Voltage_tmp(end) - F_baseline;
    % sample should be ok
    sampleOk = 1;
    
else % record when the output of read was empty to know when the force level was kept equal because of read failure
    sampleOk = 0;
    F_now_Voltage = NaN;% by default at zero at beginning if no output of read
end


end % function