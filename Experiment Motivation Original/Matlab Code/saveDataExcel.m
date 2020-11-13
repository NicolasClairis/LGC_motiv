function saveDataExcel(all,nbTrial,file_nm)
% saveDataExcel(all,nbTrial,file_nm)


% Ugly big function to put the data in the correct order to save. Prepare 4 excel sheet, one for each block
for iBlock = 1:4
    
    % for each trial they did
    for j = 1:nbTrial(iBlock)
        
        % extract the data from the cells to put them back in the order excel wants
        thesignal = all.signals{iBlock,j}; % Extract second cell.
        theWin(j) = all.win{iBlock,j};
        theVCmax(j) = all.VCmax{iBlock,j};
        theTrialLength(j) = all.trialLength{iBlock,j};
        theIncentive(j) = all.incentive{iBlock,j};
        theFirstT2(j) = all.firstT2{iBlock,j};
        
        % as this is the only vector to save, do a for loop on it
        for row = 1 : size(thesignal, 2)
            tableData{row,j} = thesignal(row);
        end
        
    end
    
    % Unique value, so we save it appart
    theUsedMVC = all.usedMVC{iBlock,1};
    
    % Prepare the data in right order for cell conversion to table
    for row = 1: nbTrial(iBlock)
        tableData{row,nbTrial(iBlock)+1} = theWin(row);
        tableData{row,nbTrial(iBlock)+2} = theVCmax(row);
        tableData{1,nbTrial(iBlock)+3} = theUsedMVC;
        tableData{row,nbTrial(iBlock)+4} = theTrialLength(row);
        tableData{row,nbTrial(iBlock)+5} = theIncentive(row);
        tableData{row,nbTrial(iBlock)+6} = theFirstT2(row);
    end
    
    % transform data into table
    tableData = cell2table(tableData);
    
    % prepare table names for columns
    trial_name = {};
    for trial_i = 1 :nbTrial(iBlock)
        trial_name{trial_i} = strcat('trial_signal_',num2str(trial_i));
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

