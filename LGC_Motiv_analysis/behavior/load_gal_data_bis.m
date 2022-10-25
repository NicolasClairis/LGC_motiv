function[excelReadGeneralFile] = load_gal_data_bis()
% [excelReadGeneralFile] = load_gal_data_bis()
%load_gal_data_bis will load summary_participants_infos.xlsx data in 
% matlab as a big table allowing you to then extract the data that you want
% from it.

%% define path
list_pcs = {'Lab','Home'};
which_pc_idx = listdlg('PromptString',{'Lab or home pc?'},...
    'SelectionMode','single','ListString',list_pcs);
switch list_pcs{which_pc_idx}
    case 'Lab'
        gitPath = fullfile('C:','Users','clairis','Desktop');
    case 'Home'
        gitPath = fullfile('C:','Users','Loco','Documents');
end
pcPath = fullfile(gitPath,'Github','LGC_motiv','LGC_Motiv_results');

%% extract the data
excelReadGeneralFile = readtable([pcPath,filesep,...
    'summary_participants_infos.xlsx']);

end % function