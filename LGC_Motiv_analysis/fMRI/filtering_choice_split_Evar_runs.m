function[outliers, allSubs] = filtering_choice_split_Evar_runs
% [outliers, allSubs] = filtering_choice_split_Evar_runs
% will aim at identifying the runs to exclude when we split high vs low
% effort choices and want to have a high effort regressor
%
% OUTPUTS
% outliers: structure indicating the subjects and the runs which lack
% diversity (at least 2 effort levels) to perform a regression on hE
% levels. The structure will also indicate you whether the lack of
% diversity is due to a lack of hE levels in the low effort choices
% (n_hE_lEch <= 1) or in the high effort choices (n_hE_hEch <= 1).
%
% allSubs: structure with number of hE levels for high effort choices
% (n_hE_hEch) and for low effort choices (n_hE_lEch) for all subjects.
%
% See also choiceNDproportion_perRun_group.m for the saturation runs and
% subjects

%% subjects to include
[study_nm, condition, gender, subject_id, NS] = sub_id;

%% working directory
rootPath = fullfile('E:',study_nm);

%% initialize variable of interest
outliers = struct;
n_hE_levels = 3;

%% loop over subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    CID_nm = ['CID',sub_nm];
    subFolder = [fullfile(rootPath,CID_nm,'behavior'),filesep];
    [runs, n_runs] = runs_definition(study_nm, sub_nm, condition);
    
    for iR = 1:n_runs
        jR = runs.runsToKeep(iR);
        run_nm = num2str(jR);
        switch runs.tasks{iR}
            case 'Em'
                task_fullName = 'mental';
            case 'Ep'
                task_fullName = 'physical';
        end
        
        %% extract choice made and hE level for the current run
        [choice_highE] = extract_choice_hE(subFolder, sub_nm, run_nm, task_fullName);
        [hE_level] = extract_hE_level(subFolder, sub_nm, run_nm, task_fullName);
        %% extract high effort (hE) level for the high effort choices (hEch)
        % and for the low effort choices (lEch)
        hE_hEch = hE_level(choice_highE == 1);
        hE_lEch = hE_level(choice_highE == 0);
        
        %% extract number of levels for each
        n_hE_hEch = size(unique(hE_hEch),2);
        n_hE_lEch = size(unique(hE_lEch),2);
        
        %% filter subjects where there is only one level
        if (n_hE_hEch <= 1) || (n_hE_lEch <= 1) % lack of diversity in hE levels
            outliers.(CID_nm).(['run',run_nm]).n_hE_hEch = n_hE_hEch;
            outliers.(CID_nm).(['run',run_nm]).n_hE_lEch = n_hE_lEch;
        end % filter
        
        %% store data for all subjects
        allSubs.(CID_nm).(['run',run_nm]).n_hE_hEch = n_hE_hEch;
        allSubs.(CID_nm).(['run',run_nm]).n_hE_lEch = n_hE_lEch;
        
        %% also store the number of trials for each occurrence
        [allSubs.(CID_nm).(['run',run_nm]).hE_hEch,...
            allSubs.(CID_nm).(['run',run_nm]).hE_lEch] = deal(NaN(1,n_hE_levels));
        for iE = 1:n_hE_levels
            allSubs.(CID_nm).(['run',run_nm]).hE_hEch(iE) = sum(hE_hEch == iE);
            allSubs.(CID_nm).(['run',run_nm]).hE_lEch(iE) = sum(hE_lEch == iE);
        end
    end % run loop
end % subject loop