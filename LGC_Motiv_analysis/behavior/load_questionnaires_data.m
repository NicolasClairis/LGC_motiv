function[excelReadQuestionnairesFile] = load_questionnaires_data()
% [excelReadQuestionnairesFile] = load_questionnaires_data()
%load_questionnaires_data will load Questionnaires-scores3.xlsx data in matlab as a big table.

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
excelReadQuestionnairesFile = readtable([pcPath,filesep,'Questionnaire-scores-final.xlsx'],...
        'Sheet','Sheet1');

end % function