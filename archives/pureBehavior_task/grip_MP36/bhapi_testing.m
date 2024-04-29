% function to test BHAPI

% go to correct folder
cd 'C:\Program Files (x86)\BIOPAC Systems, Inc\BIOPAC Hardware API 2.2.3 Research\SampleProjects\MATLAB\mpdevdemo';

% define paths
% mpdevDLLPath = 'C:\Windows\System32\mpdev.dll';
% mpdevHeaderPath = 'C:\Windows\System32\';
mpdevDLLPath = 'C:\Program Files (x86)\BIOPAC Systems, Inc\BIOPAC Hardware API 2.2.3 Research\x64\mpdev.dll';
mpdevHeaderPath = 'C:\Program Files (x86)\BIOPAC Systems, Inc\BIOPAC Hardware API 2.2.3 Research\';

% which device
mp150code = 101;
mp36R = 102;
mpDevice = mp36R;

% which method of communication with the device
mpmethod = 10;

% serial number
sn = 'auto';

% launch demo
mpdevdemo(mpdevDLLPath, mpdevHeaderPath, mpDevice, mpmethod, sn);