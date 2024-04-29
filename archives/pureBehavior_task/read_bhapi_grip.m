function [Fvalue] = read_bhapi_grip(libname)
% [Fvalue] = read_bhapi_grip(libname)
%
% INPUTS
% libname: name of the dll ('MPDEV') for the handgrip
%
% OUTPUTS
% Fvalue: force value in realtime
%

numRead = 0;
% numValuesToRead = samplingFrequency*nChannels;
numValuesToRead = 1; % 1 sample at a time
% remainingSamplesToRead = nSamples*nChannels;
remainingSamplesToRead = 1; % how many samples to read before retrieving the value
tbuff(1:numValuesToRead) = double(0); % initialize the correct amount of data
offset = 1;

while (remainingSamplesToRead > 0)
    if numValuesToRead > remainingSamplesToRead
        numValuesToRead = remainingSamplesToRead;
    end
    
    [retval, tbuff, numRead]  = calllib(libname, 'receiveMPData',tbuff, numValuesToRead, numRead);
    
    if ~strcmp(retval,'MPSUCCESS')
        fprintf(1,'Failed to receive MP data.\n');
        calllib(libname, 'disconnectMPDev');
        return
    else
        Fvalue(offset:offset+double(numRead(1))-1) = tbuff(1:double(numRead(1)));
        % careful if more than 1 channel is used, you should adapt this
        % part of the script
    end
    offset = offset + double(numValuesToRead);
    remainingSamplesToRead = remainingSamplesToRead - double(numValuesToRead);
end % while there are not the samples required, keep waiting

end % function