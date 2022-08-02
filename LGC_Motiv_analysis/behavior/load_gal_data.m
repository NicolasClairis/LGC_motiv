function[excelReadGeneralFile] = load_gal_data()
% [excelReadGeneralFile] = load_gal_data()
%load_gal_data will load General1.xlsx data in matlab as a big table.

%% define path
list_pcs = {'Lab','Home'};
which_pc_idx = listdlg('PromptString',{'Lab or home pc?'},...
    'SelectionMode','single','ListString',list_pcs);
switch list_pcs{which_pc_idx}
    case 'Lab'
        hardDisk = 'M:';
    case 'Home'
        hardDisk = 'L:';
end
pcPath = fullfile(hardDisk,'human_data_private','raw_data_subject','study1');
cd(pcPath);

%% extract the data
excelReadGeneralFile = readtable([pcPath,filesep,'General1.xlsx'],...
        'Sheet','MAIN');

end % function