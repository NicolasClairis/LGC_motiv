function saveDataExcel(all,nbTrial,file_nm)
% saveDataExcel(all,nbTrial,file_nm)
% converts matlab structure of the LGC motivation task with all the relevant
% informations into an excel file
% 
% INPUTS
%
% See also LGCM_main


% Ugly big function to put the data in the correct order to save. Prepare 4 excel sheet, one for each block
for iBlock = 1:4
    
    % for each trial they did
    for jBlock = 1:nbTrial(iBlock)
        
        % extract the data from the cells to put them back in the order excel wants
        thesignal = all.signals{iBlock,jBlock}; % Extract second cell.
        theWin(jBlock) = all.win{iBlock,jBlock};
        theVCmax(jBlock) = all.VCmax{iBlock,jBlock};
        theTrialLength(jBlock) = all.trialLength{iBlock,jBlock};
        theIncentive(jBlock) = all.incentive{iBlock,jBlock};
        theFirstT2(jBlock) = all.firstT2{iBlock,jBlock};
        
        % as this is the only vector to save, do a for loop on it
        for iRow = 1 : size(thesignal, 2)
            tableData{iRow,jBlock} = thesignal(iRow);
        end
        
    end
    
    % Unique value, so we save it appart
    theUsedMVC = all.usedMVC{iBlock,1};
    
    % Prepare the data in right order for cell conversion to table
    for iRow = 1: nbTrial(iBlock)
        tableData{iRow,nbTrial(iBlock)+1} = theWin(iRow);
        tableData{iRow,nbTrial(iBlock)+2} = theVCmax(iRow);
        tableData{1,nbTrial(iBlock)+3} = theUsedMVC;
        tableData{iRow,nbTrial(iBlock)+4} = theTrialLength(iRow);
        tableData{iRow,nbTrial(iBlock)+5} = theIncentive(iRow);
        tableData{iRow,nbTrial(iBlock)+6} = theFirstT2(iRow);
    end
    
    % transform data into table
    tableData = cell2table(tableData);
    
    % prepare table names for columns
    trial_name = {};
    for iTtrial = 1 :nbTrial(iBlock)
        trial_name{iTtrial} = strcat('trial_signal_',num2str(iTtrial));
    end
    trial_name{nbTrial(iBlock)+1} = 'doWin';
    trial_name{nbTrial(iBlock)+2} = 'VCmax';
    trial_name{nbTrial(iBlock)+3} = 'usedVCmax';
    trial_name{nbTrial(iBlock)+4} = 'duration_Nbframe';
    trial_name{nbTrial(iBlock)+5} = 'incentive_value';
    trial_name{nbTrial(iBlock)+6} = 'FirstT2_sec';
    tableData.Properties.VariableNames = trial_name;
    
    % put the table into the sheet and clear the table to prepare it for a new sheet (due to cells being annoying)
    writetable(tableData, [file_nm,'.xlsx'],'Sheet',iBlock)
    clear tableData;
end

% Move the data in the correct file
movefile(['D:' filesep 'Matlab codes' filesep 'Experiment Motivation Original' filesep 'Matlab Code' filesep subjectCodeName '_data.xlsx'],['D:' filesep 'Matlab codes' filesep 'Experiment Motivation Original' filesep 'SavedData'])

end

