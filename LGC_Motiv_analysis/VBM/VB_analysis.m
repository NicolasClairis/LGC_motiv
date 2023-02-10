clc
clear all
close all

%% iniate spm
spm('defaults','fmri');
spm_jobman('initcfg');

%% checking = 0 if you want to run directly or = 1 if you want to look at the batch first
checking = 0;

%% Setup the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 8);

% Specify sheet and range
opts.Sheet = "Sheet1";
opts.DataRange = "A2:H73";

% Specify column names and types
opts.VariableNames = ["CID", "kR", "kP", "Kpe", "kme", "bias", "kpf", "kce"];
opts.VariableTypes = ["char", "double", "double", "double", "double", "double", "double", "double"];

% Specify variable properties
opts = setvaropts(opts, "CID", "WhitespaceRule", "preserve");
opts = setvaropts(opts, "CID", "EmptyFieldRule", "auto");

% Import the data
CIDanatforbehavior = readtable("D:\anat_analysis_study1\CID_anat_for_behavior.xlsx", opts, "UseExcel", false);
clear opts
kR = double(table2array(CIDanatforbehavior(:,2)));


%% define working directory
main_folder = 'D:\anat_analysis_study1';
nb_files = length(kR);
matlabbatch{1}.spm.stats.factorial_design.dir = {'D:\anat_analysis_study1\VBM_kR'};

%% prepare the which files to select
j_files = 1;
for i_files = 1:length(kR)
    
    cd(main_folder)
    all_files = dir('CID*');
    if 1-isnan(kR(i_files))
        %load anat filenames
        CID_filename = all_files(i_files).name;
        cd(CID_filename)
        % select the UNI-DEN
        smwc_files = dir('*smwc1*');
        smwc_filename = smwc_files(1).name;
        
        anaPath{j_files,1} = strcat(main_folder,'\',CID_filename,'\',smwc_filename,',1');
        kR_wo_nan(j_files) = kR(i_files);
        j_files=j_files+1;
    end
end
cd(main_folder)
%%
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = anaPath;
matlabbatch{1}.spm.stats.factorial_design.cov.c = kR_wo_nan';
%%
matlabbatch{1}.spm.stats.factorial_design.cov.cname = 'kR';
matlabbatch{1}.spm.stats.factorial_design.cov.iCFI = 1;
matlabbatch{1}.spm.stats.factorial_design.cov.iCC = 1;
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {''};
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;


%% display spm batch before running it or run it directly
if checking == 1
    spm_jobman('interactive',matlabbatch);
    %     spm_jobman('run',matlabbatch);
elseif checking == 0
    spm_jobman('run',matlabbatch);
end