% function to test BHAPI

% go to correct folder
cd 'C:\Program Files (x86)\BIOPAC Systems, Inc\BIOPAC Hardware API 2.2 Research\SampleProjects\MATLAB\mpdevdemo';

% define paths
% % mpdevDLLPath = 'C:\Users\clairis\Desktop\GitHub\LGC_motiv\pureBehavior_task\grip_MP36\BIOPAC Hardware API 2.2 Research\Win32\mpdev.dll';
% % mpdevHeaderPath = 'C:\Users\clairis\Desktop\GitHub\LGC_motiv\pureBehavior_task\grip_MP36\BIOPAC Hardware API 2.2 Research\';
% mpdevDLLPath = 'C:\Windows\System32\mpdev.dll';
% mpdevHeaderPath = 'C:\Windows\System32\';
mpdevDLLPath = 'C:\Program Files (x86)\BIOPAC Systems, Inc\BIOPAC Hardware API 2.2 Research\Win32\mpdev.dll';
mpdevHeaderPath = 'C:\Program Files (x86)\BIOPAC Systems, Inc\BIOPAC Hardware API 2.2 Research\';

% launch demo
mpdevdemo(mpdevDLLPath, mpdevHeaderPath, 103,10,'auto');