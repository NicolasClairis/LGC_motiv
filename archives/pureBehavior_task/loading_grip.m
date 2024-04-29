function [libname] = loading_grip(samplingRate)
% [libname] = loading_grip(samplingRate)
%
% loading_grip will load MP36R handgrip to prepare to acquire data based on the sampling
% rate initialized in samplingRate parameter.
%
% INPUTS
% samplingRate = number of samples to extract within 1 second
%
% OUTPUTS
% libname: library name (mpdev) to use
% 

%% define paths
libname = 'mpdev';
mpdevDLLPath = 'C:\Program Files (x86)\BIOPAC Systems, Inc\BIOPAC Hardware API 2.2.3 Research\x64\mpdev.dll';
mpdevHeaderPath = 'C:\Program Files (x86)\BIOPAC Systems, Inc\BIOPAC Hardware API 2.2.3 Research\';

%% loading dll and library
% which device
mp150code = 101;
mp36R = 102;
mpDevice = mp36R;

% which method of communication with the device
mpmethod = 10;

% serial number
sn = 'auto';

loadlibrary(mpdevDLLPath, strcat(mpdevHeaderPath,'mpdev.h'));
fprintf(1,'\nMPDEV.DLL LOADED!!!\n');

%% connecting the grip
retval = calllib(libname,'connectMPDev',mpDevice,mpmethod,sn);

if ~strcmp(retval,'MPSUCCESS')
    fprintf(1,'Failed to Connect.\n');
    calllib(libname, 'disconnectMPDev');
    sca;
    return
end

fprintf(1,'Connected\n');

%% set sampling rate
retval = calllib(libname, 'setSampleRate', samplingRate);

if ~strcmp(retval,'MPSUCCESS')
    fprintf(1,'Failed to Set Sample Rate.\n');
    calllib(libname, 'disconnectMPDev');
    sca;
    return
else
    fprintf(1,['Sampling rate updated to ',num2str((1/samplingRate)*1000),' Hz']);
end

%% set acquisition channel 1
% define which channel is to be recorded
acqChannel1 = 1;
acqChannel2 = 0;
acqChannel3 = 0;
acqChannel4 = 0;

switch mpDevice
    case mp150code % MP150
        aCH = [int32(acqChannel1),int32(acqChannel2),int32(acqChannel3),int32(acqChannel4),int32(0),int32(0),int32(0),int32(0),int32(0),int32(0),int32(0),int32(0),int32(0),int32(0),int32(0),int32(0)];
    case mp36R % mp36R
        aCH = [int32(acqChannel1),int32(acqChannel2),int32(acqChannel3),int32(acqChannel4)];
    otherwise
        error(['mpDevice ',num2str(mpDevice),' not on the list']);
end

retval = calllib(libname, 'setAcqChannels',aCH);

if ~strcmp(retval,'MPSUCCESS')
    fprintf(1,'Failed to Set Acq Channels.\n');
    calllib(libname, 'disconnectMPDev');
    sca;
    return
end

fprintf(1,'Acquisition on Channel 1 Set\n');

end % function