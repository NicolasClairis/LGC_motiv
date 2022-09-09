function[ROI_trial_b_trial] = extract_ROI_betas_onsets_only(computerRoot, study_nm, subject_id, condition, GLM)
% [ROI_trial_b_trial] = extract_ROI_betas_onsets_only(computerRoot, study_nm, subject_id, condition, GLM)
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
switch computerRoot
    case 'E:\' % lab computer
        gitFolder = fullfile('C:','Users','clairis','Desktop','GitHub','LGC_motiv','Matlab_DIY_functions','ROI');
    case 'L:\human_data_private\raw_data_subject\' % home computer
        gitFolder = fullfile('C:','Users','Loco','Documents','GitHub','LGC_motiv','Matlab_DIY_functions','ROI');
    otherwise
        error(['computer root ',computerRoot,' not ready yet.'])
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

%% select the ROI to use
[ ROI_xyz, ~, ROI_names, n_ROIs ] = ROI_selection(gitFolder);

%% which GLM
if ~exist('GLM','var') || isempty(GLM)
    listOfAllOnsetsOnlyGLM = [64];
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

%% initialize variables of interest
nTrialsPerRun = 54;
potentialTimePeriods = {'preChoiceCross',...
    'choice','chosen',...
    'preEffortCross','Eperf','fbk'};
nPotentialTimePeriods = length(potentialTimePeriods);
tasks = {'Ep','Em'}; nTasks = length(tasks);
nRunsPerTask = 2;
for iROI = 1:n_ROIs
    ROI_nm = ROI_names.(['ROI_',num2str(iROI)]);
    for iTask = 1:nTasks
        task_nm = tasks{iTask};
        for iRun = 1:nRunsPerTask
            run_str = ['run',num2str(iRun)];
            for iTperiod = 1:nPotentialTimePeriods
                timePeriod_nm = potentialTimePeriods{iTperiod};
                curr_onset_nm = GLMprm.model_onset.(task_nm).(timePeriod_nm);
                if ~strcmp(curr_onset_nm,'none') && ismember(curr_onset_nm,{'stick','boxcar'})
                    ROI_trial_b_trial.(ROI_nm).(task_nm).(run_str).(timePeriod_nm) = NaN(nTrialsPerRun, NS);
                else
                    error(['problem with ',curr_onset_nm])
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
    ROI_nm = ROI_names.(ROI_nb_nm);
    sxyz_ROI = ROI_xyz.(ROI_nb_nm);

    for iS = 1:NS
        sub_nm = subject_id{iS};
        subj_folder = [fullfile(dataRoot, ['CID',sub_nm]),filesep];
        sub_fMRI_path = [fullfile(subj_folder, 'fMRI_analysis','functional',...
            ['preproc_sm_',num2str(preproc_sm_kernel),'mm'],['GLM',GLMstr]),filesep];
        subj_behavior_folder = [subj_folder, 'behavior' filesep];
        [subRuns, n_subRuns] = runs_definition(study_nm, sub_nm, condition);

        %% extract index of runs and trials that were included in the GLM
        totalTrials_idx = []; % index of trials kept

        for iRun = subRuns.runsToKeep
            run_str = num2str(iRun);
            run_nm = ['run',run_str];
            % extract trials where no choice was made
            currRunBehaviorFileName = ls([subj_behavior_folder,'*_session',run_str,'_*_task.mat']);
            if size(currRunBehaviorFileName,1) > 1
                error(['problem file identification: too many files popping out with run number ',run_str]);
            end
            if strcmp(currRunBehaviorFileName(16:23),'physical') ||...
                    strcmp(currRunBehaviorFileName(17:24),'physical')
                task_behavioral_id = 'physical';
                task_behavioral_id_bis = 'physicalPerf';
            elseif strcmp(currRunBehaviorFileName(16:21),'mental') ||...
                    strcmp(currRunBehaviorFileName(17:22),'mental')
                task_behavioral_id = 'mental';
                task_behavioral_id_bis = 'mentalE_perf';
            else
                error('problem in identifying task type because file name doesn''t match');
            end
            behavioralDataStruct = load([subj_behavior_folder, currRunBehaviorFileName]);
            T0 = behavioralDataStruct.onsets.T0;
            choiceOnsets = behavioralDataStruct.(task_behavioral_id_bis).onsets.choice - T0;
            choiceMissedTrials = isnan(choiceOnsets);

            % trialN: index of the trials
            trialN = 1:nTrialsPerRun;
            trialN(choiceMissedTrials) = [];
            totalTrials_idx = [totalTrials_idx; repmat(trialN,1,nTimePeriods)']; % multiply vector by as many conditions
            % ignore movement
            totalTrials_idx = [totalTrials_idx; NaN(n_mvmt,1)];
        end % run loop

        totalTrials_idx = [totalTrials_idx; NaN(n_subRuns,1)]; % ignore run constant regressors
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
        jRun = 0;
        jTimePhase = 0;
        for iTrial = 1:n_totalTrialsForFirstLevel
            jRunTrial = totalTrials_idx(iTrial);
            if ~isnan(jRunTrial) % ignore movement and run constant regressors
                beta_idx = beta_idx + 1;
                % identify which run and time period we are at now
                if iTrial == 1
                    jRun = 1;
                    jTimePhase = 1;
                elseif ~isnan(totalTrials_idx(iTrial-1)) && ~isnan(totalTrials_idx(iTrial)) &&...
                        (totalTrials_idx(iTrial) - totalTrials_idx(iTrial-1)) < 0 % start of new task phase (trial number gets re-initialized)
                    jTimePhase = jTimePhase + 1;
                elseif isnan(totalTrials_idx(iTrial-1)) && ~isnan(totalTrials_idx(iTrial)) % start of new run
                    jRun = jRun + 1;
                    jTimePhase = 1;
                end
                task_nm = subRuns.tasks{jRun};
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
                ROI_trial_b_trial.(ROI_nm).(task_nm).(run_nm).(timePeriod_nm)(jRunTrial, iS) = beta_value; % save mean beta for the selected ROI inside the big resulting matrix
            else % skip movement regressors from beta extraction
                beta_idx = beta_idx + 1;

            end % ignore mvmt and run cstt
        end % trial number loop
        %% indicator subject done
        disp(['Subject ',num2str(iS),'/',num2str(NS),' extracted']);
    end % subject loop
    disp(['ROI ',num2str(iROI),'/',num2str(n_ROIs),' done']);
end % ROI loop

save(['ROI_BOLDperTrial_GLM',GLM_str,'.mat'],...
    'ROI_trial_b_trial','ROI_names','ROI_xyz')
% end % function