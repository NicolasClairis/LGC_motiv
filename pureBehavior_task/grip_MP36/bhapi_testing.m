% function to test BHAPI

% go to correct folder
cd 'C:\Program Files (x86)\BIOPAC Systems, Inc\BIOPAC Hardware API 2.2.3 Research\SampleProjects\MATLAB\mpdevdemo';

% define paths
% mpdevDLLPath = 'C:\Windows\System32\mpdev.dll';
% mpdevHeaderPath = 'C:\Windows\System32\';
mpdevDLLPath = 'C:\Program Files (x86)\BIOPAC Systems, Inc\BIOPAC Hardware API 2.2.3 Research\x64\mpdev.dll';
mpdevHeaderPath = 'C:\Program Files (x86)\BIOPAC Systems, Inc\BIOPAC Hardware API 2.2.3 Research\';

% launch demo
mpdevdemo(mpdevDLLPath, mpdevHeaderPath, 103,10,'auto');