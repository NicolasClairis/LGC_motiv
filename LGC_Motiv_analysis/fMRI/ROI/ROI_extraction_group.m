function[con_vec_all,...
    con_avg, con_sem, con_sd,...
    con_names,...
    ROI_coords, ttest_ROI] = ROI_extraction_group(study_nm, GLM,...
    subject_id, condition, fig_disp)
% [con_vec_all,...
%     con_avg, con_sem, con_sd,...
%     con_names,...
%     ROI_coords, ttest_ROI] = ROI_extraction_group(study_nm, GLM,...
%     subject_id, condition, fig_disp)
% ROI_extraction_group will serve to extract the ROI data across all
% participants for a given GLM number.
%
% INPUTS
% study_nm:
% 'study1': dmPFC/aIns study
% 'study2': clinical trial
%
% GLM: GLM number for which you want to get the ROI
%
% subject_id: list of subject names, if left empty it will select all
% subjects
%
% condition: define subjects and runs to include
%
% fig_disp:
%(0) no display
%(1) display figure
%
% OUTPUTS
% con_vec_all: contrast*participant*ROI matrix with all the data
%
% con_avg: mean ROI contrast value for each participant
%
% con_sem: sem ROI contrast value for each participant
%
% con_sd: standard deviation ROI contrast value for each participant
%
% con_names: list of contrast names
%
% ROI_coords: structure with informations regarding ROI coordinates and
% name
%
% ttest_ROI: structure with summary statistics over ROI


%% extract working directories
[computerRoot] = LGCM_root_paths();

if ~strcmp(study_nm,'study1')
    error('case not ready yet');
end

dataRoot = [computerRoot,filesep,study_nm,filesep];

switch computerRoot
    case 'E:\' % lab computer
        gitFolder = fullfile('C:','Users','clairis','Desktop','GitHub','LGC_motiv','Matlab_DIY_functions','ROI');
    case 'L:\human_data_private\raw_data_subject\' % home computer
        gitFolder = fullfile('C:','Users','Loco','Documents','GitHub','LGC_motiv','Matlab_DIY_functions','ROI');
end
ROI_path = [dataRoot,filesep,'results',filesep,'ROI',filesep];

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
[ ROI_xyz, ~, ROI_nm, n_ROIs ] = ROI_selection(gitFolder);

%% which GLM
if ~exist('GLM','var') || isempty(GLM)
    GLM = spm_input('fMRI GLM number?',1,'e');
end
GLMstr = num2str(GLM);
GLMprm = which_GLM(GLM);

%% beta or t value?
% beta_or_t_values = {'beta_value','t_value'};
% b_or_t = spm_input('beta value or t value?',1,'b','beta value| t.value',[1 2]);
% beta_or_t_value  = beta_or_t_values{b_or_t};
beta_or_t_value = 'beta_value';

%% preprocessing smoothing kernel
% preproc_sm_kernel  = spm_input('preprocessing smoothing kernel?',1,'r');
preproc_sm_kernel = 8; % by default

%% prepare contrasts
[con_names] = LGCM_contrasts(study_nm, subject_id{1}, GLM,...
    computerRoot, preproc_sm_kernel, condition);
n_max_con = length(con_names);

%% how many figures do you want to plot
if fig_disp == 1
    n_figs = spm_input('How many figures?',1,'e');
    
    % prepare contrast names for question
    list_con_question = con_names{1};
    for iCon = 2:n_max_con
        list_con_question = [list_con_question,' | ',con_names{iCon}];
    end
    list_con_question = [list_con_question,' | If all wished contrasts for this figure selected click here'];
    % select contrasts to display at the end
    which_con = zeros(n_figs, n_max_con); % column with 0/1 to see which contrasts to display at the end
end

%% name + contrasts to display for each figure
if fig_disp == 1
    conName = cell(n_figs,1);
    for iFig = 1:n_figs
        %% what name for each figure?
        conName{iFig} = spm_input(['Name for figure ',num2str(iFig),' please?'],1,'s');
        
        %% select which contrasts you want to display for each figure
        stop_con_loop = 0;
        while stop_con_loop == 0
            selectedContrast = spm_input(['What contrast for fig.',num2str(iFig),':',conName{iFig},' ?'],1,'m',...
                list_con_question, ...
                1:(n_max_con+1), 0);
            if selectedContrast < n_max_con + 1
                which_con(iFig,selectedContrast) = 1; % 1 for selected contrasts
            elseif selectedContrast == n_max_con + 1 % stop contrast selection loop
                stop_con_loop = 1;
            end
        end % contrast loop
        
        %% select which contrasts you want to display for each figure
        n_con_fig_tmp = sum(which_con(iFig,:));
        figCon_idx = find(which_con(iFig,:));
        figConNames.(['fig',num2str(iFig)]) = cell(1,n_con_fig_tmp);
        for iCon = 1:n_con_fig_tmp
            figConNames.(['fig',num2str(iFig)]){iCon} = spm_input([con_names{figCon_idx(iCon)},...
                ' short name:'],1,'s');
        end % contrast loop
    end % figure loop
end % figure display

%% create big matrix to store ROI values
con_vec_all = NaN(n_max_con, NS, n_ROIs);
    
%% extract the ROI for each contrast of the selected subjects
for iROI = 1:n_ROIs
    fprintf('Region %d\n',iROI);
    
    sxyz_ROI = ROI_xyz.(['ROI_',num2str(iROI)]);
    %% loop through subjects for each study
    for iS = 1:NS
        sub_nm = subject_id{iS};
        subPath = fullfile(dataRoot,['CID',sub_nm],...
            'fMRI_analysis','functional',...
            ['preproc_sm_',num2str(preproc_sm_kernel),'mm']);
        [sub_fMRI_path] = fMRI_subFolder(subPath, GLM, condition);
        % extract contrasts for the current subject
        [con_names_currSub] = LGCM_contrasts(study_nm, sub_nm, GLM,...
            computerRoot, preproc_sm_kernel, condition);
        
        % extract contrast for the current subject in the current ROI
        % for each contrast
        for iCon = 1:n_max_con
            jCon = find(strcmp(con_names(iCon), con_names_currSub));
            if ~isempty(jCon) % if participant doesn't have the corresponding contrast, just ignore this subject
                [ con_value ] = ROI_extraction( '', jCon, sub_fMRI_path, sxyz_ROI, beta_or_t_value );
                con_vec_all(iCon, iS, iROI) = con_value; % save mean beta for the selected ROI inside the big resulting matrix
            end
        end % contrast
        
        %% indicator subject done
        disp(['Subject ',num2str(iS),'/',num2str(NS),' extracted']);
    end % subject loop
    
    %% extract coordinates informations
    ROI_coords.(['ROI_',num2str(iROI)]).xyz = sxyz_ROI;
    ROI_coords.ROI_nm = ROI_nm;
    ROI_coords.n_ROIs = n_ROIs;
end % roi loop

%% average across subjects
[con_avg, con_sem, con_sd] = deal(NaN(n_max_con, n_ROIs));
for iCon = 1:n_max_con
    for iROI = 1:n_ROIs
        [con_avg(iCon, iROI),...
            con_sem(iCon, iROI),...
            con_sd(iCon, iROI)] = mean_sem_sd(con_vec_all(iCon,:,iROI), 2);
    end
end

%% t.test significativity of each contrast
ttest_pval = NaN(n_max_con, n_ROIs);
for iROI = 1:n_ROIs
    for iCon = 1:n_max_con
        curr_con = NaN(1,NS);
        curr_con(1,:) = con_vec_all(iCon, :, iROI);
        [~, ttest_pval(iCon, iROI),~,stats_tmp] = ttest( curr_con );
        ttest_ROI.p_value(iCon, iROI) = ttest_pval(iCon, iROI);
        ttest_ROI.t_value(iCon, iROI) = stats_tmp.tstat;
    end % contrast loop
end % ROI loop

%% save all the data
filename = [ROI_path,beta_or_t_value,...
    '_GLM',GLMstr,'_',num2str(n_ROIs),'ROIs_',num2str(NS),'subs_',condition];
if ~exist([filename,'.mat'],'file')
        save([filename,'.mat'])
    else
        while exist([filename,'.mat'],'file')
            filename = [filename,'_bis'];
        end
        save([filename,'.mat'])
end

%% loop through figures as defined when script launched at the
% beginning
if fig_disp == 1
    for iFig = 1:n_figs
        conFig = which_con(iFig,:);
        selectedCon = find(conFig == 1); % extract index of contrasts selected
        % display figure
        [roi_fig] = roi_graph(selectedCon, n_ROIs,...
            con_avg, con_sem, figConNames.(['fig',num2str(iFig)]), ttest_pval);
        % save image
        cd(ROI_path);
        img_name = ['GLM',num2str(GLM),'_',beta_or_t_value,...
            '_',conName{iFig},'_',num2str(n_ROIs),'ROIs_',num2str(NS),'_subs.png'];
        if exist(img_name,'file') ~= 2
            set(gcf,'PaperPosition',[0 0 1 1]);
            set(gcf,'PaperPositionMode','auto');
            saveas(roi_fig,img_name);
        else
            warning(['There is already a file with the name ', img_name,' so the figure has not been saved yet.',...
                ' Please do it manually or change the script to be able to do it automatically in the future.']);
        end
    end % figure loop
end % figure display

end % function