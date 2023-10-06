function[ROI_trial_b_trial] = extract_ROI_betas_onsets_only(computerRoot,...
    study_nm, subject_id, condition, GLM, ROI_infos)
% [ROI_trial_b_trial] = extract_ROI_betas_onsets_only(computerRoot,...
%   study_nm, subject_id, condition, GLM, ROI_infos)
% extract_ROI_betas_onsets_only will extract the activity of the ROI
% defined trial/trial and for each time point of the GLM that has been
% used.
%
% INPUTS
% computerRoot: path for computer root (will be asked if left empty)
%
% study_nm: study name ('study1' by default)
%
% subject_id: list of subjects (will be asked if left empty)
%
% condition: define which subjects to include (will be asked if left empty)
%
% GLM: GLM number to use for the trial/trial ROI activity extraction
%
% ROI_infos: structure containing ROI information (matrix coordinates, ROI
% name, number of ROIs to extract, etc.). Can be left empty and will be
% asked with ROI_selection.m.
%
% OUTPUTS
% ROI_trial_b_trial: structure with ROI BOLD deconvoluted trial/trial and
% for each phase of the GLM.

%% extract working directories
if ~exist('computerRoot','var') || isempty(computerRoot)
    [computerRoot] = LGCM_root_paths();
end

if ~exist('study_nm','var') || isempty(study_nm)
    study_nm = 'study1';
end
dataRoot = [computerRoot,filesep,study_nm,filesep];
% git folder only required if you need to define the ROI to use
if ~exist('ROI_infos','var') || isempty(ROI_infos)
    switch computerRoot
        case 'E:\' % lab computer
            gitFolder = fullfile('C:','Users','clairis','Desktop','GitHub','LGC_motiv','Matlab_DIY_functions','ROI');
        case 'L:\human_data_private\raw_data_subject\' % home computer
            gitFolder = fullfile('C:','Users','Loco','Documents','GitHub','LGC_motiv','Matlab_DIY_functions','ROI');
        otherwise
            error(['gitFolder path with computerRoot = ',computerRoot,' not ready yet.'])
    end
end
%% define subject list
if ~exist('condition','var') || isempty(condition)
        condition = subject_condition;
end
if ~exist('subject_id','var') || isempty(subject_id)
    [subject_id, NS] = LGCM_subject_selection(study_nm, condition);
else
    NS = length(subject_id);
end
% store list of subjects
ROI_trial_b_trial.subject_id = subject_id;

%% select the ROI to use
if ~exist('ROI_infos','var') || isempty(ROI_infos)
    % do you want to use the MRS ROI or a literature/GLM-fMRI-based ROI?
    ROI_type = spm_input('ROI category?',1,'m',...
        ['literature/fMRI-based ROI |'...
        'dmPFC MRS voxel |'...
        'AI MRS voxel'], ...
        1:3, 1);
    switch ROI_type
        case 1 % literature/fMRI-based ROI => need to select
            [ ROI_xyz, ~, ROI_names, n_ROIs ] = ROI_selection(gitFolder);
        case 2 % MRS dmPFC voxel
            ROI_names.ROI_1_shortName = 'MRS_dmPFC';
            ROI_names.ROI_1 = 'MRS_dmPFC';
            ROI_names.fullpath.ROI_1 = 'subject_specific';
            n_ROIs = 1;
        case 3 % MRS anterior insula voxel
            ROI_names.ROI_1_shortName = 'MRS_aINS';
            ROI_names.ROI_1 = 'MRS_aINS';
            ROI_names.fullpath.ROI_1 = 'subject_specific';
            n_ROIs = 1;
    end
else
    ROI_xyz = ROI_infos.ROI_xyz;
    ROI_names = ROI_infos.ROI_nm;
    n_ROIs = ROI_infos.n_ROIs;
end

%% which GLM
if ~exist('GLM','var') || isempty(GLM)
    listOfAllOnsetsOnlyGLM = [64, 65, 70,...
        90, 93, 94, 95, 128, 129,...
        158, 186, 201];
    nPossibleGLMs = size(listOfAllOnsetsOnlyGLM, 2);
    listGLM = ['GLM',num2str(listOfAllOnsetsOnlyGLM(1))];
    if nPossibleGLMs > 1
        for iGLM = 2:nPossibleGLMs
            listGLM = [listGLM,' | GLM',num2str(listOfAllOnsetsOnlyGLM(iGLM))];
        end
    end
    GLM_idx = spm_input('GLM number?',1,'m',...
        listGLM,1:nPossibleGLMs,0);
    GLM = listOfAllOnsetsOnlyGLM(GLM_idx);
end
GLMstr = num2str(GLM);
if GLM == 64
    error(['GLM ',GLMstr,' seems to be badly defined as regressors are not ',...
        'orthogonal (probably you can''t model all events in the same GLM']);
end
GLMprm = which_GLM(GLM);
% check GLM parameters
add_drv = GLMprm.gal.add_drv; % if derivative has been added, need to double the number of regressors
if add_drv > 0
    error('add_drv > 0 not ready yet for onsets_only extraction');
end
% is the GLM an onsets-only GLM?
if GLMprm.gal.onsets_only == 0
    error(['You selected GLM ',num2str(GLM),' where onsets_only = 0.'])
end
n_mvmt = 6;

%% preprocessing smoothing kernel
% preproc_sm_kernel  = spm_input('preprocessing smoothing kernel?',1,'r');
preproc_sm_kernel = 8; % by default

%% bias-field
biasFieldCorr = 0;

%% initialize variables of interest
nTrialsPerRun = 54;
potentialTimePeriods = {'preChoiceCross',...
    'choice','chosen',...
    'preEffortCross','Eperf','fbk'};
nPotentialTimePeriods = length(potentialTimePeriods);
tasks = {'Ep','Em'}; nTasks = length(tasks);
nRunsPerTask = 2;
for iROI = 1:n_ROIs
    ROI_nm1 = ROI_names.(['ROI_',num2str(iROI)]);
    for iTask = 1:nTasks
        task_nm = tasks{iTask};
        for iRun = 1:nRunsPerTask
            run_str = ['run',num2str(iRun)];
            for iTperiod = 1:nPotentialTimePeriods
                timePeriod_nm = potentialTimePeriods{iTperiod};
                curr_onset_nm = GLMprm.model_onset.(task_nm).(timePeriod_nm);
                if ~strcmp(curr_onset_nm,'none') && ismember(curr_onset_nm,{'stick','boxcar'})
                    ROI_trial_b_trial.(ROI_nm1).(task_nm).(run_str).(timePeriod_nm) = NaN(nTrialsPerRun, NS);
                elseif ~ismember(curr_onset_nm,{'none','stick','boxcar'})
                    error(['problem with time modulation = ',curr_onset_nm])
                end
            end % time period
        end % run loop
    end % task
end % ROI loop

nTimePeriods = 0;
timePeriods = {};
for iTperiod = 1:nPotentialTimePeriods
    timePeriod_nm = potentialTimePeriods{iTperiod};
    curr_onset_nm_Ep = GLMprm.model_onset.Ep.(timePeriod_nm);
    curr_onset_nm_Em = GLMprm.model_onset.Em.(timePeriod_nm);
    if strcmp(curr_onset_nm_Ep, curr_onset_nm_Em) &&...
            ismember(curr_onset_nm_Ep,{'stick','boxcar'}) &&...
            ismember(curr_onset_nm_Em,{'stick','boxcar'})
        nTimePeriods = nTimePeriods + 1;
        timePeriods = [timePeriods, timePeriod_nm];
    end
end % time period

%% extract the ROI for each event
for iROI = 1:n_ROIs
    disp(['starting ROI ',num2str(iROI),'/',num2str(n_ROIs)]);
    ROI_nb_nm = ['ROI_',num2str(iROI)];
    ROI_nm2 = ROI_names.(ROI_nb_nm);
    % load ROI coordinates
    if ROI_type == 1
        sxyz_ROI = ROI_xyz.(ROI_nb_nm);
    end
    
    for iS = 1:NS
        sub_nm = subject_id{iS};
        sub_fullNm = ['CID',sub_nm];
        %         if ~strcmp(sub_nm,'036')
        subj_folder = [fullfile(dataRoot, sub_fullNm),filesep];
        switch biasFieldCorr
            case 0
                sub_folderName = [fullfile(subj_folder, 'fMRI_analysis','functional',...
                    ['preproc_sm_',num2str(preproc_sm_kernel),'mm']),filesep];
            case 1
                sub_folderName = [fullfile(subj_folder, 'fMRI_analysis','functional',...
                    ['preproc_sm_',num2str(preproc_sm_kernel),'mm_with_BiasFieldCorrection']),filesep];
        end
        sub_fMRI_path = fMRI_subFolder(sub_folderName, GLM, condition);
        if ~exist(sub_fMRI_path,'dir')
            error(['ROI extraction not ready for ',condition,' yet.']);
        end
        subj_behavior_folder = [subj_folder, 'behavior' filesep];
        [subRuns, n_subRuns] = runs_definition(study_nm, sub_nm, condition);
        
        % load individual ROI coordinates when MRS voxel was selected
        if ismember(ROI_type,[2,3])
            sxyz_ROI = load_MRS_ROI_coords(study_nm, sub_nm, ROI_nm2);
            ROI_xyz.(ROI_nb_nm).(sub_fullNm) = sxyz_ROI;
        end
        
        if ismember(ROI_type,[1,2]) || (ROI_type == 3 && ~isempty(sxyz_ROI)) % filter case where no insula could be extracted in MRS
            
            %% extract index of runs and trials that were included in the GLM
            totalTrials_idx = []; % index of trials kept
            runPerTrial = [];
            
            for iRun = 1:n_subRuns
                jRun = subRuns.runsToKeep(iRun);
                run_str = num2str(jRun);
                run_nm = ['run',run_str];
                % extract trials where no choice was made
                currRunBehaviorFileName = ls([subj_behavior_folder,...
                    '*_session',run_str,'_*_task.mat']);
                if size(currRunBehaviorFileName,1) > 1
                    error(['problem file identification: too many files popping out with run number ',run_str]);
                end
                if strcmp(currRunBehaviorFileName(16:23),'physical') ||...
                        strcmp(currRunBehaviorFileName(17:24),'physical')
                    task_behavioral_id = 'Ep';
                    task_behavioral_id_bis = 'physicalPerf';
                elseif strcmp(currRunBehaviorFileName(16:21),'mental') ||...
                        strcmp(currRunBehaviorFileName(17:22),'mental')
                    task_behavioral_id = 'Em';
                    task_behavioral_id_bis = 'mentalE_perf';
                else
                    error('problem in identifying task type because file name doesn''t match');
                end
                behavioralDataStruct = load([subj_behavior_folder, currRunBehaviorFileName]);
                T0 = behavioralDataStruct.onsets.T0;
                choiceOnsets = behavioralDataStruct.(task_behavioral_id_bis).onsets.choice - T0;
                
                choiceMissedTrials = isnan(choiceOnsets);
%                 warning('careful cause 1st level has been changed and missed trials are now included => just be careful to match');
%                 choiceMissedTrials = false(size(choiceOnsets));
                
                % trialN: index of the trials
                trialN = 1:nTrialsPerRun;
                trialN(choiceMissedTrials) = [];
                totalTrials_idx = [totalTrials_idx; repmat(trialN,1,nTimePeriods)']; % multiply vector by as many conditions
                runPerTrial = [runPerTrial; repmat(jRun, 1, nTimePeriods*length(trialN))'];
                % ignore movement
                totalTrials_idx = [totalTrials_idx; NaN(n_mvmt,1)];
                runPerTrial = [runPerTrial; NaN(n_mvmt,1)];
            end % run loop
            
            totalTrials_idx = [totalTrials_idx; NaN(n_subRuns,1)]; % ignore run constant regressors
            runPerTrial = [runPerTrial; NaN(n_subRuns,1)]; % ignore run constant regressors
            n_totalTrialsForFirstLevel = length(totalTrials_idx); % number of trials kept for the 1st level
            
            %% check whether the number of trials corresponds to the
            % number of betas extracted by the GLM
            betas_nm = ls([sub_fMRI_path,'beta_*.nii']); % extract name of all betas
            % extract number of betas
            n_betas = size(betas_nm,1);
            % there should be as many betas as there are trials + movement
            % betas + run constant betas
            if n_betas ~= n_totalTrialsForFirstLevel
                error(['Number of betas extracted by the GLM (',num2str(n_betas),')',...
                    ' and number of trials extracted through onsets (',num2str(n_totalTrialsForFirstLevel),')',...
                    ' are not consistent for subject ',sub_nm,'. Please check what happened and come back.']);
            end
            
            %% extract beta values for each trial
            cd(sub_fMRI_path);
            beta_idx    = 0;
            
            %% loop through trials
            kRun = 0;
            jTimePhase = 0;
            for iTrial = 1:n_totalTrialsForFirstLevel
                jRunTrial = totalTrials_idx(iTrial);
                % convert run number
                switch runPerTrial(iTrial)
                    case {1,2}
                        run_nm = 'run1';
                    case {3,4}
                        run_nm = 'run2';
                end
                if ~isnan(jRunTrial) % ignore movement and run constant regressors
                    beta_idx = beta_idx + 1;
                    % identify which run and time period we are at now
                    if iTrial == 1
                        kRun = 1;
                        jTimePhase = 1;
                    elseif ~isnan(totalTrials_idx(iTrial-1)) && ~isnan(totalTrials_idx(iTrial)) &&...
                            (totalTrials_idx(iTrial) - totalTrials_idx(iTrial-1)) < 0 % start of new task phase (trial number gets re-initialized)
                        jTimePhase = jTimePhase + 1;
                    elseif isnan(totalTrials_idx(iTrial-1)) && ~isnan(totalTrials_idx(iTrial)) % start of new run
                        kRun = kRun + 1;
                        jTimePhase = 1;
                    end
                    task_nm = subRuns.tasks{kRun};
                    timePeriod_nm = timePeriods{jTimePhase};
                    
                    % extract con/scon files
                    if beta_idx < 10
                        betaZeros = '000';
                    elseif beta_idx >= 10 && beta_idx < 100
                        betaZeros = '00';
                    elseif beta_idx >= 100 && beta_idx < 1000
                        betaZeros = '0';
                    elseif beta_idx >= 1000
                        betaZeros = '';
                    end
                    
                    betaNum     = strcat(['beta_', betaZeros, num2str(beta_idx), '.nii']);
                    betaVol     = spm_vol(betaNum); % extracts header to read con file
                    betadata    = spm_read_vols(betaVol); % reads the con file
                    
                    % extract indices inside the con file based on the
                    % MNI coordinates entered previously for the
                    % sphere/mask
                    vxyz        = unique(floor((inv(betaVol.mat) * sxyz_ROI')'), 'rows'); % converts from MNI (mm) to voxel-space
                    vi          = sub2ind(betaVol.dim, vxyz(:, 1), vxyz(:, 2), vxyz(:, 3)); % extracts coordinates of the sphere in the con (voxel-space)
                    beta_value  = mean(betadata(vi),'omitnan'); % extracts mean beta for the selected ROI sphere/mask
                    ROI_trial_b_trial.(ROI_nm2).(task_nm).(run_nm).(timePeriod_nm)(jRunTrial, iS) = beta_value; % save mean beta for the selected ROI inside the big resulting matrix
                else % skip movement regressors from beta extraction
                    beta_idx = beta_idx + 1;
                    
                end % ignore mvmt and run cstt
            end % trial number loop
            %         else
            %             warning('check 036 with last upgrades.');
            %         end % filter subject 036 who for some reason doesn't work
        end % filter for anterior insula MRS voxel which was not extracted in all participants
        
        %% indicator subject done
        disp(['Subject ',num2str(iS),'/',num2str(NS),' extracted']);
    end % subject loop
    disp(['ROI ',num2str(iROI),'/',num2str(n_ROIs),' done']);
end % ROI loop

save([dataRoot,'results',filesep,'ROI',filesep,...
    'ROI_BOLDperTrial_GLM',GLMstr,'_',num2str(NS),'subs.mat'],...
    'ROI_trial_b_trial','ROI_names','ROI_xyz');
end % function