function[pb_subs] = debug_2nd_level_batch_metabolites_bis(GLM, condition, gender, iCon)
% [pb_subs] = debug_2nd_level_batch_metabolites_bis(GLM, condition, gender, iCon)
% script to launch contrast iCon but removing all subjects one by one
% in order to detect which subjects are problematic and to perform those 
% which are not problematic.
%
% INPUTS
% GLM: number of the GLM
%
% condition: define subjects and runs to include (you can also use
% subject_condition.mat to select the condition)
%
% gender:
% 'all': all subjects by default
% 'males': remove females from the list
% 'females': remove males from the list
%
% checking:
% (0) launch script directly
% (1) show GUI and allow user to check before launching
%
% iCon: contrast number
%
% pb_subs: list of subjects which, when removed, remove the bug

%% clear workspace
close all; clc;

%% define study and list of subjects to include
% define study
if ~exist('study_nm','var')
    %     study_names = {'fMRI_pilots','study1','study2'};
    %     study_nm_idx = listdlg('ListString',study_names);
    %     study_nm = study_names{study_nm_idx};
    study_nm = 'study1'; % by default
end

% define subjects
if ~exist('condition','var') || ~strcmp(condition(1:4),'fMRI')
    condition = subject_condition;
end
if ~exist('gender','var') ||...
        isempty(gender) ||...
        ~ismember(gender,{'all','males','females'})
    gender = 'all';
end
[subject_id, NS] = LGCM_subject_selection(study_nm, condition, gender);

%% ask which metabolite to focus on?
[low_met_subs, high_met_subs, metabolite_nm,...
    MRS_ROI_nm, metabolite_allSubs] = medSplit_metabolites(study_nm, subject_id);
NS_low = sum(low_met_subs);
NS_high = sum(high_met_subs);
%% remove subjects who may have NaN values for the metabolite selected
badSubs = isnan(metabolite_allSubs);
n_badSubs = sum(badSubs);
if n_badSubs > 0
    badSub_idx = find(badSubs ~= 0);
    for iBS = 1:n_badSubs
        badSub_nm = subject_id{badSub_idx(iBS)};
        disp(['Subject ',badSub_nm,' had to be removed as ',...
            MRS_ROI_nm,' ',metabolite_nm,' was NaN.']);
    end
    subject_id(badSub_idx) = [];
    NS = NS - n_badSubs;
    low_met_subs(badSub_idx) = [];
    high_met_subs(badSub_idx) = [];
end
%% string with number of subjects
NS_str = num2str(NS);
%% initiate SPM
spm('defaults','fmri');
spm_jobman('initcfg');

%% check the batch before launching the script?
if ~exist('checking','var') ||...
        isempty(checking) ||...
        ~ismember(checking,[0,1])
    checking = 0;
end

%% value of the smoothing during preprocessing?
preproc_sm_kernel = 8;

%% working directories
% computer_root = LGCM_root_paths();
computer_root = ['E:',filesep];
% scripts_folder = fullfile(computer_root,'GitHub','LGC_motiv','LGC_Motiv_analysis','fMRI');
% addpath(scripts_folder);
switch study_nm
    case 'fMRI_pilots'
        studyRoot = fullfile(computer_root,'fMRI_pilots');
    case 'study1'
        studyRoot = fullfile(computer_root,'study1');
    case 'study2'
        studyRoot = fullfile(computer_root,'study2');
end

%% define GLM
if ~exist('GLM','var') || isempty(GLM) || GLM <= 0
    GLM = spm_input('GLM number?',1,'e',[]);
end
GLM_str = num2str(GLM);
GLMprm = which_GLM(GLM);

%% create results folder
mainPath = [studyRoot,filesep,'Second_level',filesep];
MRS_split_nm = [MRS_ROI_nm,'_',metabolite_nm,'MedSplit'];
switch condition
    case {'fMRI'}
        results_folder = [mainPath,...
            'GLM',GLM_str,'_',NS_str,'subs_',...
            'preprocSm',num2str(preproc_sm_kernel),'mm_',MRS_split_nm,filesep];
    case 'fMRI_noSatRunSub'
        results_folder = [mainPath,...
            'GLM',GLM_str,'_',NS_str,'subs_',...
            'preprocSm',num2str(preproc_sm_kernel),'mm_noSatRunSubs_',MRS_split_nm,filesep];
    case 'fMRI_noSatTaskSub'
        results_folder = [mainPath,...
            'GLM',GLM_str,'_',NS_str,'subs_',...
            'preprocSm',num2str(preproc_sm_kernel),'mm_noSatTaskSubs_',MRS_split_nm,filesep];
    case 'fMRI_noMoveSub'
        results_folder = [mainPath,...
            'GLM',GLM_str,'_',NS_str,'subs_',...
            'preprocSm',num2str(preproc_sm_kernel),'mm_noFullMvmtSub_',MRS_split_nm,filesep];
    case 'fMRI_noMoveSub_bis'
        results_folder = [mainPath,...
            'GLM',GLM_str,'_',NS_str,'subs_',...
            'preprocSm',num2str(preproc_sm_kernel),'mm_noMvmtSubLenient_',MRS_split_nm,filesep];
    case 'fMRI_noMoveSub_ter'
        results_folder = [mainPath,...
            'GLM',GLM_str,'_',NS_str,'subs_',...
            'preprocSm',num2str(preproc_sm_kernel),'mm_noMvmtSubStringent_',MRS_split_nm,filesep];
    case 'fMRI_noSatTaskSub_noMove_bis_Sub'
        results_folder = [mainPath,...
            'GLM',GLM_str,'_',NS_str,'subs_',...
            'preprocSm',num2str(preproc_sm_kernel),'mm_noSatTaskNoMvmtSub_',MRS_split_nm,filesep];
    case 'fMRI_noSatTask'
        results_folder = [mainPath,...
            'GLM',GLM_str,'_',NS_str,'subs_',...
            'preprocSm',num2str(preproc_sm_kernel),'mm_noSatTask_',MRS_split_nm,filesep];
    case 'fMRI_noSatRun'
        results_folder = [mainPath,...
            'GLM',GLM_str,'_',NS_str,'subs_',...
            'preprocSm',num2str(preproc_sm_kernel),'mm_noSatRun_',MRS_split_nm,filesep];
    case 'fMRI_noSatRun_choiceSplit_Elvl'
        results_folder = [mainPath,...
            'GLM',GLM_str,'_',NS_str,'subs_',...
            'preprocSm',num2str(preproc_sm_kernel),'mm_noSatRun_bis_',MRS_split_nm,filesep];
    case 'fMRI_noMove_bis'
        results_folder = [mainPath,...
            'GLM',GLM_str,'_',NS_str,'subs_',...
            'preprocSm',num2str(preproc_sm_kernel),'mm_noMvmtRunLenient_',MRS_split_nm,filesep];
    case 'fMRI_noMove_ter'
        results_folder = [mainPath,...
            'GLM',GLM_str,'_',NS_str,'subs_',...
            'preprocSm',num2str(preproc_sm_kernel),'mm_noMvmtRunStringent_',MRS_split_nm,filesep];
    case 'fMRI_noSatTask_noMove_bis'
        results_folder = [mainPath,...
            'GLM',GLM_str,'_',NS_str,'subs_',...
            'preprocSm',num2str(preproc_sm_kernel),'mm_noSatTaskNoMvmtRun_',MRS_split_nm,filesep];
    otherwise
        error(['folder not ready yet for the condition ',condition]);
end

% create folder to store the results
if ~exist(results_folder,'dir')
    mkdir(results_folder);
else
%     error(['Folder ',results_folder,' already exists.']);
end

%% 1) take mean anatomy across participants
batch_idx = 0;

%% perform the second level

% load all contrasts of interest
con_names = LGCM_contrasts(study_nm, subject_id{1}, GLM,...
    computer_root, preproc_sm_kernel, condition);

% loop over contrasts
okSubs = zeros(1,NS);
for iS = 1:NS
    
    % remove selectionned subject
    subject_id_bis = subject_id;
    low_met_subs_bis = low_met_subs;
    high_met_subs_bis = high_met_subs;
    subject_id_bis(iS) = [];
    NS_bis = NS - 1;
    low_met_subs_bis(iS) = [];
    high_met_subs_bis(iS) = [];
    
    current_con_nm = con_names{iCon};
    
    % directory for concatenated contrast
    conFolder_nm = strrep(current_con_nm,' ','_');
    conFolder_nm = strrep(conFolder_nm,':','');
    conFolder_nm = [results_folder, conFolder_nm,'_without_CID',subject_id{iS}];
    if ~exist(conFolder_nm,'dir')
        mkdir(conFolder_nm);
        
        % start second level:
        matlabbatch{1}.spm.stats.factorial_design.dir = {conFolder_nm};
        % list of inputs
        conlist_low = cell(NS_low,1); % 1 con per EPI-subject
        conlist_high = cell(NS_high,1); % 1 con per EPI-subject
        % re-initialize counters
        jS_low = 0;
        jS_high = 0;
        % extract contrasts per subject
        for iS_bis = 1:NS_bis
            sub_nm = subject_id_bis{iS_bis};
            
            % check incompatibility between some GLM and some subjects to see
            % if you need to redefine the list of subjects included in the
            % analysis
            checkGLM_and_subjectIncompatibility(study_nm, sub_nm, condition, GLMprm);
            subject_main_folder = [studyRoot,filesep,'CID',sub_nm, filesep, 'fMRI_analysis' filesep,...
                'functional' filesep, 'preproc_sm_',num2str(preproc_sm_kernel),'mm',filesep];
            subject_main_folder = fMRI_subFolder(subject_main_folder, GLM, condition);
            if isempty(subject_main_folder)
                error(['condition ',condition,' not planned yet. Please add it.']);
            end
            
            %% adapt contrast index since some conditions and contrasts are missing for some subjects
            con_names_perSub = LGCM_contrasts(study_nm, sub_nm, GLM,...
                computer_root, preproc_sm_kernel, condition);
            if sum(strcmp(current_con_nm, con_names_perSub)) > 0
                % extract index (for this subject) of the current contrast
                jCon = find(strcmp(current_con_nm, con_names_perSub));
                
                % extract name for this particular subject of the contrast of
                % interest
                [con_str] = conNumber2conName(jCon);
                if low_met_subs_bis(iS_bis) == 1 && high_met_subs_bis(iS_bis) == 0
                    jS_low = jS_low + 1;
                    conlist_low(jS_low) = {[subject_main_folder,'con_',con_str,'.nii,1']};
                elseif low_met_subs_bis(iS_bis) == 0 && high_met_subs_bis(iS_bis) == 1
                    jS_high = jS_high + 1;
                    conlist_high(jS_high) = {[subject_main_folder,'con_',con_str,'.nii,1']};
                else
                    error('Who''s that dude with weird metabolites?');
                end
            end % in case contrast exists for the current subject
        end % subject loop
        
        % remove empty subjects from the list for the current contrast
        conlist_low(cellfun(@isempty,conlist_low)) = [];
        conlist_high(cellfun(@isempty,conlist_high)) = [];
        
        matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = conlist_low;
        matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = conlist_high;
        
        % default parameters
        matlabbatch{1}.spm.stats.factorial_design.des.t2.dept = 0;
        matlabbatch{1}.spm.stats.factorial_design.des.t2.variance = 1;
        matlabbatch{1}.spm.stats.factorial_design.des.t2.gmsca = 0;
        matlabbatch{1}.spm.stats.factorial_design.des.t2.ancova = 0;
        matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
        matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
        matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
        matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
        matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
        
        % model estimation
        matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File',...
            substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
        matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
        matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
        % contrasts
        matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File',...
            substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
        % contrast t.test high>low
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = [current_con_nm,'_H_min_l',metabolite_nm];
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = [-1 1];
        matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
        % contrast t.test low>high
        matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = [current_con_nm,'_l_min_H',metabolite_nm];
        matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = [1 -1];
        matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
        % contrast: F.test
        % TO BE ADDED
        matlabbatch{3}.spm.stats.con.delete = 0;
        try
            spm_jobman('run',matlabbatch);
        catch
            disp(['subject ',subject_id{iS},' doesn''t solve the pb']);
            okSubs(iS) = 1; % keep track of the subjects who are not creating the bug
        end
    end
end % loop over subjects

% extract problematic subjects
pb_subs = subject_id(okSubs == 0);


end % function