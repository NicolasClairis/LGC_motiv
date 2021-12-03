% developed by Emmanuelle Bioud - october 2017
% adapted by Nicolas Clairis - december 2021

computer_root = LGCM_root_paths();
if ~exist('study_nm','var') || isempty(study_nm)
    study_names = {'fMRI_pilots','study1','study2'};
    study_nm_idx = listdlg('ListString',study_names);
    study_nm = study_names{study_nm_idx};
end
switch study_nm
    case 'fMRI_pilots'
        root = fullfile(computer_root,'fMRI_pilots');
    case 'study1'
        root = fullfile(computer_root,'study1');
    case 'study2'
        root = fullfile(computer_root,'study2');
end

storageFolder = [root, filesep 'movement_summary' filesep];
if ~exist(storageFolder,'dir')
    mkdir(storageFolder);
end

[subject_id, NS] = LGCM_subject_selection(study_nm);

yScale = [-3 3];

for iS = 1:NS
    sub_nm = subject_id{iS};
    path_SubFolders = [root, filesep, 'CID', sub_nm, filesep 'fMRI_scans'];
    cd(path_SubFolders);
    listRunFolders = ls('*run*');
    
    n_files = size(listRunFolders, 1);
    n_runs = size(listRunFolders, 1);
    
    for iRun = 1:n_runs
        
        cd(listRunFolders(iRun,:));
        rpFile = ls('rp_*');
        
        
        printfig = fig();
        set(printfig, 'Name', ['Motion parameters: subject CID ', sub_nm, ', Run ', num2str(iRun) ], 'Visible', 'on');
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
        
        
        figName = [storageFolder, 'motionPlot_CID', sub_nm, '_run', num2str(iRun), '.png'];
        print(printfig, '-dpng', '-noui', '-r100', figName);  % enable to print to file
        close(printfig);   % enable to close graphic window
        
        
        cd(path_SubFolders);
        
    end
    
end