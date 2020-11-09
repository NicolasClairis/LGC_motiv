% developed by Emmanuelle Bioud - october 2017
clear all
rootDir = 'C:\Users\emmanuelle.bioud\Documents\MBB_DB_IRMf\Data fMRI\NIFTI\';
storageFolder = 'C:\Users\emmanuelle.bioud\Documents\MBB_DB_IRMf\Movement plots\';

[subject_id, NS] = MS2_subject_id_selection('fMRI');
nSub = size(listSubFolder, 1);

yScale = [-3 3];

for iS = 1:NS
    
    path_SubFolders{iS} = [rootDir, subject_id{iS}, filesep];
    cd(path_SubFolders{iS});
    listRunFolders = ls('*RUN*');
    
    for i = 1:size(listRunFolders, 1)
        idxRefblips(i) = ~isempty(strfind(listRunFolders(i,:), 'REFBLIP'));
    end
    
    if ~isempty(idxRefblips)
        listRunFolders(idxRefblips,:) = [];
    end
    
    
    for iRun = 1:size(listRunFolders, 1)
        
        cd(listRunFolders(iRun,:));
        rpFile = ls('rp_*');
        
        
        printfig = figure;
        set(printfig, 'Name', ['Motion parameters: Subject ', num2str(iS), ', Run ', num2str(iRun) ], 'Visible', 'on');
        loadmot = load(deblank(rpFile));
        subplot(2,1,1);
        plot(loadmot(:,1:3));
        legend('x', 'y', 'z    .', 'Location', 'NorthEastOutside')
        grid on;
        ylim(yScale);  % enable to always scale between fixed values as set above
        title('Motion parameters: shifts (top, in mm) and rotations (bottom, in dg)', 'interpreter', 'none');
        subplot(2,1,2);
        plot(loadmot(:,4:6)*180/pi);
        legend('pitch', 'roll', 'yaw', 'Location', 'NorthEastOutside')
        grid on;
        ylim(yScale);   % enable to always scale between fixed values as set above
        
        
        figName = [storageFolder, 'motionPlot_sub', num2str(iS), '_run', num2str(iRun), '.png'];
        print(printfig, '-dpng', '-noui', '-r100', figName);  % enable to print to file
        close(printfig);   % enable to close graphic window
        
        
        cd(path_SubFolders{iS});
        
    end
    
end