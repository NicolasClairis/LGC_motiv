function[excelReadQuestionnairesFile, sub_CID_list] = load_questionnaires_data()
% [excelReadQuestionnairesFile, sub_CID_list] = load_questionnaires_data()
%load_questionnaires_data will load Questionnaires-scores3.xlsx data in matlab as a big table.
%
% OUTPUT
% excelReadQuestionnairesFile: table with all data in columns
%
% sub_CID_list: correct CID as a string

%% define path
list_pcs = {'Lab','Home'};
which_pc_idx = listdlg('PromptString',{'Lab or home pc?'},...
    'SelectionMode','single','ListString',list_pcs);
switch list_pcs{which_pc_idx}
    case 'Lab'
        pcPath = fullfile('M:','human_data_private',...
            'raw_data_subject','study1');
    case 'Home'
        pcPath = 'P:\boulot\postdoc_CarmenSandi\results\questionnaires';
end
% cd(pcPath);

%% extract the data
excelReadQuestionnairesFile = readtable([pcPath,filesep,'Questionnaire-scores-final.xlsx'],...
    'Sheet','Sheet1');

%% fix name of subjects
NS = length(excelReadQuestionnairesFile.CID);
sub_CID_list = cell(1,NS);
for iS = 1:NS
    sub_idx = excelReadQuestionnairesFile.CID(iS);
    if sub_idx < 10
        sub_CID_list{iS} = ['00',num2str(sub_idx)];
    elseif sub_idx >= 10 && sub_idx < 100
        sub_CID_list{iS} = ['0',num2str(sub_idx)];
    elseif sub_idx >= 100
        sub_CID_list{iS} = num2str(sub_idx);
    end
end % subject loop

end % function