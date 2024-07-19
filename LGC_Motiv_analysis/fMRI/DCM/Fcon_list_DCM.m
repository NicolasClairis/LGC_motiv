function[fcon_names, fcon_matrix] = Fcon_list_DCM(study_nm, sub_nm, GLM, DCM_mode, computer_root, preproc_sm_kernel, condition, biasFieldCorr)
% [fcon_names, fcon_matrix] = Fcon_list_DCM(study_nm, sub_nm, GLM, DCM_mode, computer_root, preproc_sm_kernel, condition, biasFieldCorr)
% [fcon_names, fcon_matrix] = LGCM_contrasts(study_nm, sub_nm, GLM,...
%   computer_root, preproc_sm_kernel, condition, biasFieldCorr, DCM_GLM)
% LGCM_contrasts will define the contrast names and contrast vector for the
% subject and study entered in input.
%
% INPUTS
% computer_root: path where data is stored
%
% study_nm: definition of the study on which you want to analyze the data
% 'fMRI_pilots': pilots
% 'study1': first study (dmPFC + AI)
% 'study2': second study (clinical trial)
%
% sub_nm: subject name
%
% GLM: GLM number
%
% DCM_mode:
% (1) all sessions modeled independently like in a classic univariate GLM
% => hard to manipulate for DCM but could be useful for testing
% session-specific effects or comparing sessions
% (2) sessions pooled within each task (ex: session 1 and 3 of physical
% effort will be concatenated into one single regressor) but each task will
% be modeled separately
% (3) all sessions pooled together
% (4) all trial periods are pooled together across sessions except for
% choice and effort which are modeled independently for each task (but
% pooled across sessions of the same task)
% (5) all trial periods are pooled together across sessions except for
% the effort period which is modeled independently for each task (but
% pooled across sessions of the same task)
%
% preproc_sm_kernel: kernel used in preprocessing for smoothing the data
%
% condition: define subjects and runs to include
% 'fMRI': all subjects where fMRI ok
% 'fMRI_no_move': remove runs with too much movement
%
% biasFieldCorr: binary variable indicating whether results are based on
% bias-field corrected files (1) or not (0)
%
% DCM_GLM: use preprocessing adapted for DCM or not (the main difference is
% that the preprocessing will be including or not a slice-timing correction)
% (0) use basic preprocessing
% (1) use preprocessing including slice-timing correction
%
% OUTPUTS
% fcon_names: list of contrast names
%
% fcon_matrix: nLines*nRegs*nContrasts matrix with all contrasts

%% working directories
switch study_nm
    case 'fMRI_pilots'
        root = fullfile(computer_root,'fMRI_pilots');
    case 'study1'
        root = fullfile(computer_root,'study1');
    case 'study2'
        root = fullfile(computer_root,'study2');
end
subj_folder             = [root, filesep, 'CID',sub_nm];
switch biasFieldCorr
    case 0
        biasField_sufix = '';
    case 1
        biasField_sufix = '_with_BiasFieldCorrection';
    otherwise
        error('biasFieldCorr not defined in the input. Please enter a value');
end
subj_analysis_folder    = [subj_folder, filesep, 'fMRI_analysis' filesep, 'functional',filesep,...
    'preproc_sm_',num2str(preproc_sm_kernel),'mm',biasField_sufix,'_DCM',filesep];
[resultsFolderName] = fMRI_subFolder_DCM(subj_analysis_folder, GLM, condition, DCM_mode);

%% extract run information
[runs, n_runs] = runs_definition(study_nm, sub_nm, condition);

%% extract GLM informations
dispGLM = 0;
[reg_names, n_totalRegs] = GLM_details_DCM(GLM, DCM_mode, dispGLM, runs, n_runs);

%% control number of regressors
betaFiles = ls([resultsFolderName,'beta_*.nii']);
n_1stLevelBetas = size(betaFiles,1);
if n_1stLevelBetas ~= n_totalRegs
    error(['Problem: number of regressors predicted by contrasts = ',num2str(n_totalRegs),...
        ' while number of betas produced by 1st level = ',num2str(n_1stLevelBetas),...
        ' for subject ',sub_nm]);
% else
%     disp(['CID',sub_nm,' number of regressors and number of betas produced in 1st level - ok'])
end

%% check DCM_mode
switch DCM_mode
    case 1
        error(['DCM_mode = ',num2str(DCM_mode),' but contrasts not adapted yet to regroup runs together']);
    case {2,4,5}
        warning(['DCM_mode = ',num2str(DCM_mode),' meaning Ep/Em modeled independently ',...
            ' but they are not regrouped in the contrasts for now.',...
            '=> consider pooling them']);
end

%% initialize variables of interest
fcon_names = {};
fcon_matrix = struct;

%% extract "effects-of-interest" F.contrast

% identify regressors of 'no-interest'
mvmt_reg_nm = {'movement'};
run_cstt_reg_nm = {'r1 constant','r2 constant','r3 constant','r4 constant'};
allCross_reg_nm = {'ONSET fixation cross r1','ONSET fixation cross r2',...
    'ONSET fixation cross r3','ONSET fixation cross r4',...
    'ONSET fixation cross Ep','ONSET fixation cross Em',...
    'ONSET fixation cross'};
preChoiceCross_reg_nm = {'ONSET preChoice white fixation cross r1','ONSET preChoice white fixation cross r2',...
    'ONSET preChoice white fixation cross r3','ONSET preChoice white fixation cross r4',...
    'ONSET preChoice white fixation cross Ep','ONSET preChoice white fixation cross Em',...
    'ONSET preChoice white fixation cross'};
RP_names = {'RP','R','P'};
n_RP = length(RP_names);
E_names = {'E',...
    'E1','E2','E3',...
    'Ech0','Ech1','Ech2','Ech3',...
    'lEch','hEch'};
n_E = length(E_names);
n_reg_types = 7*n_RP*n_E;
[preECross_reg_nm,...
    dispC_reg_nm,...
    fbk_reg_nm] = deal(cell(1,n_reg_types));
for iRP = 1:n_RP
    RP_nm = RP_names{iRP};
    for iE = 1:n_E
        E_nm = E_names{iE};
        reg_idx = (1:7) + 7*(iE-1) + 7*n_E*(iRP-1);
        preECross_reg_nm(reg_idx) = {['ONSET preEffort black cross ',RP_nm,' ',E_nm,' r1'],...
            ['ONSET preEffort black cross ',RP_nm,' ',E_nm,' r2'],...
            ['ONSET preEffort black cross ',RP_nm,' ',E_nm,' r3'],...
            ['ONSET preEffort black cross ',RP_nm,' ',E_nm,' r4'],...
            ['ONSET preEffort black cross ',RP_nm,' ',E_nm,' Ep'],...
            ['ONSET preEffort black cross ',RP_nm,' ',E_nm,' Em'],...
            ['ONSET preEffort black cross ',RP_nm,' ',E_nm,'']};
        dispC_reg_nm(reg_idx) = {['ONSET chosen ',RP_nm,' ',E_nm,' r1'],...
            ['ONSET chosen ',RP_nm,' ',E_nm,' r2'],...
            ['ONSET chosen ',RP_nm,' ',E_nm,' r3'],...
            ['ONSET chosen ',RP_nm,' ',E_nm,' r4'],...
            ['ONSET chosen ',RP_nm,' ',E_nm,' Ep'],...
            ['ONSET chosen ',RP_nm,' ',E_nm,' Em'],...
            ['ONSET chosen ',RP_nm,' ',E_nm,'']};
        fbk_reg_nm(reg_idx) = {['ONSET feedback ',RP_nm,' ',E_nm,' r1'],...
            ['ONSET feedback ',RP_nm,' ',E_nm,' r2'],...
            ['ONSET feedback ',RP_nm,' ',E_nm,' r3'],...
            ['ONSET feedback ',RP_nm,' ',E_nm,' r4'],...
            ['ONSET feedback ',RP_nm,' ',E_nm,' Ep'],...
            ['ONSET feedback ',RP_nm,' ',E_nm,' Em'],...
            ['ONSET feedback ',RP_nm,' ',E_nm,'']};
    end % E loop
end % RP loop

% set to 1 everything but movement and run constants (to clean DCM)
iCon = 1;
con_nm = ['fcon_',conNumber2conName(iCon)];
fcon_names{iCon} = 'effects_of_interest_Fcon1';
fcon_matrix.(con_nm) = eye(n_totalRegs);
% remove movement + run constants
regs_to_ignore = ismember(reg_names,[mvmt_reg_nm, run_cstt_reg_nm]);
fcon_matrix.(con_nm)(regs_to_ignore,:) = [];

% set to 1 everything but movement, run constants and fixation cross events (to clean DCM)
iCon = iCon + 1;
con_nm = ['fcon_',conNumber2conName(iCon)];
fcon_names{iCon} = 'effects_of_interest_Fcon2';
fcon_matrix.(con_nm) = eye(n_totalRegs);
% remove movement + run constants
regs_to_ignore = ismember(reg_names,[mvmt_reg_nm, run_cstt_reg_nm,...
    allCross_reg_nm, preChoiceCross_reg_nm, preECross_reg_nm]);
fcon_matrix.(con_nm)(regs_to_ignore,:) = [];

% set to 1 everything but movement, run constants, fixation cross, dispChosen+ fbk events (to clean DCM)
iCon = iCon + 1;
con_nm = ['fcon_',conNumber2conName(iCon)];
fcon_names{iCon} = 'effects_of_interest_Fcon3';
fcon_matrix.(con_nm) = eye(n_totalRegs);
% remove movement + run constants
regs_to_ignore = ismember(reg_names,[mvmt_reg_nm, run_cstt_reg_nm,...
    allCross_reg_nm, preChoiceCross_reg_nm, preECross_reg_nm,...
    dispC_reg_nm, fbk_reg_nm]);
fcon_matrix.(con_nm)(regs_to_ignore,:) = [];

%% extract each contrast of interest
for iReg = 1:n_totalRegs
    reg_nm = reg_names{iReg};
    if ~ismember(reg_nm,{'movement',...
            'r1 constant','r2 constant','r3 constant','r4 constant'})
        iCon = iCon + 1;
        con_nm = ['fcon_',conNumber2conName(iCon)];
        fcon_names{iCon} = reg_nm;
        fcon_matrix.(con_nm)(1,:) = zeros(1, n_totalRegs);
        fcon_matrix.(con_nm)(1,iReg) = 1;
    end % regressors to ignore
end % loop over regressors

%% combine regressors depending on DCM_mode (to be done)

end % function