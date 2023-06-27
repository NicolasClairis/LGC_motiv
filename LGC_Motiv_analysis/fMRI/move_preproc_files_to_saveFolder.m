function[] = move_preproc_files_to_saveFolder(root, study_nm,...
    subject_id, NS, smKernel, biasFieldCorr)
% move_preproc_files_to_saveFolder(root, study_nm,...
%     subject_id, NS, smKernel, biasFieldCorr)
% move_preproc_files_to_saveFolder will move the files generated after
% performing preprocessing to a corresponding folder in the initial subject
% scan folders allowing to use different smoothing kernels or with/without
% bias field correction without interfering with previous preprocessing
% pipelines.
%
% INPUTS
% root: path where data is stored
%
% study_nm: string with study name 'study1'/'study2'/pilots
%
% subject_id: list of subjects to move
%
% NS: number of subjects to move
%
% smKernel: smoothing kernel used during the preprocessing
%
% biasFieldCorr: (0) if not bias-field corrected or (1) if bias-field
% corrected

for iS = 1:NS
    sub_nm = subject_id{iS};
    switch study_nm
        case {'study1','study2'}
            sub_fullNm = ['CID',sub_nm];
        case {'fMRI_pilots','study2_pilots'}
            sub_fullNm = sub_nm;
    end
    subj_scans_folder = [root, sub_fullNm, filesep,'fMRI_scans'];
    subj_scan_folders_names = ls([subj_scans_folder,'*run*']); % takes all functional runs folders
    % remove AP/PA top-up corrective runs when they were performed (only 2
    % first pilots)
    if strcmp(study_nm,'fMRI_pilots') && ismember(sub_nm,{'pilot_s1','pilot_s2'})
        [subj_scan_folders_names] = clear_topup_fromFileList(subj_scan_folders_names);
    end
    %% define number of sessions to analyze
    if strcmp(study_nm,'fMRI_pilots') &&...
            ismember(sub_nm, {'pilot_s1','pilot_s2','pilot_s3'}) % only 2 sessions for these pilots
        n_runs = 2;
    elseif strcmp(study_nm,'study1') &&...
            ismember(sub_nm,{'040'}) % fMRI had to be crashed during run 3
        n_runs = 2;
    else
        n_runs = 4;
    end
    
    for iRun = 1:n_runs % loop through runs for 3 ratings, 3 choices 1D, 3 choices 2D runs
        runPath = [subj_scan_folders_names(iRun,:),filesep]; % go to run folder
        switch biasFieldCorr
            case 0 % no bias-field correction
                preproc_newFolder_nm = ['preproc_sm_',num2str(smKernel),'mm'];
                prefixToUse = 'swr';
                % s is for smoothed data
                % 'w' for normalized data
                % and 'r' for realigned
            case 1 % with bias-field correction
                preproc_newFolder_nm = ['preproc_sm_',num2str(smKernel),'mm_with_BiasFieldCorrection'];
                prefixToUse = 'swbr'; % 'b' is for bias-field corrected
        end
        if ~exist([runPath,preproc_newFolder_nm],'dir')
            mkdir([runPath,preproc_newFolder_nm]);
        else
            error(['preprocessing folder ',preproc_newFolder_nm,' already exists for subject ',sub_fullNm,' run ',num2str(iRun)]);
            % note: something should be done above to avoid re-doing the
            % preprocessing for the subjects where it was already done
        end
        
        if strcmp(study_nm,'fMRI_pilots')
            if ismember(sub_nm,{'pilot_s1'})
                filenames = ls([runPath,'*',prefixToUse,'LGCM*.nii']);
            elseif ismember(sub_nm,{'pilot_s2'})
                filenames = ls([runPath,'*',prefixToUse,'run*.nii']);
            elseif ismember(sub_nm,{'pilot_s3'})
                filenames = ls([runPath,'*',prefixToUse,'ABNC*.img']);
                filenames = [filenames; ls([runPath,'*',prefixToUse,'ABNC*.hdr'])];
            else
                filenames = ls([runPath,'*',prefixToUse,'CID*.nii']);
            end
        elseif strcmp(study_nm,'study2_pilots')
            if ismember(sub_nm,'fMRI_pilot1_AC')
                filenames = ls([runPath,'*',prefixToUse,'AC*.nii']);
            end
        else
            filenames = ls([runPath,'*',prefixToUse,'CID*.nii']);
            %             error('please check the format (nii/img) and the start of the name of each run because it has to be stabilized now...');
        end
        
        % move files
        for iFile = 1:length(filenames)
            movefile([runPath,filenames(iFile,:)],...
                [runPath,preproc_newFolder_nm]);
        end
    end % run loop
    disp(['Subject ',num2str(iS),'/',num2str(NS),' preprocessing files ',...
        'have been correctly moved']);
end % subject loop

disp('Preprocessing moving of files is done');

end % function