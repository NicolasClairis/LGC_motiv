function[] = ROI_mask_generator()
% function to create a ROI mask based on multiple .nii images (can be more
% than one)
%
% 1) a) asks you the number of images you want to pool for the conjunction
% 1) b) asks you if you want to sum images (useful if you want to pool left
% and right image in particular) or if you prefer to make a conjunction
% (extract only a common area between two clusters, like if one ios a
% functional and the other an anatomical cluster)
%
% 2) asks you to select the images you want to pool
%
% 3) a) asks where you want to store the resulting file and under which
% name
% 3) b) does the job for you and saves them in the path you want with the name
% you want
%
% Remarks:
% 1) This script only works with .nii images currently. If you have .hdr images
% and want to use this script, please see img2nii_spm.m function developed
% by Kiyotaka Nemoto or any other function to make the conversion for you.
% 2) This scripts require SPM to be installed in Matlab paths
%
% Nicolas Clairis - august 2017


% iniate spm
spm('defaults','fmri');
spm_jobman('initcfg');


%% 1) how many images to pool?
% nb_raw_img = spm_input('How many images do you want to pool?',1,'e');
nb_raw_img = 2; % by default

% sum or conjunction?
sum_or_cunj = spm_input('sum, conjunction or AAL extraction?',1,'m','sum | conjunction | AAL extraction',[1 2 3], 0);
% identification number for the AAL area you want to extract
if sum_or_cunj == 3
    AAL_area_id_nber = spm_input('AAL area number? X1 X2',1,'e');
    % if you want to pool two areas (left right for ex) just enter the two values of the AAL atlas
    % (if of course they are proximal numbers, otherwise you have to change
    % the way the script is made)
    % otherwise, you can just enter one value
    %
    % ex: I want to pool left (33) and right (34) insula, I enter [33 34]
    % but if i just want to extract the left insula for ex. i can just
    % enter [33].
end

%% 2) select the images you want to pool
% directory containing all the raw images
root = ['B:',filesep,'resultats',filesep,'ROI_definition'];
if sum_or_cunj ~= 3
    [raw_imgs, select_ok] = spm_select(nb_raw_img,'.nii','Select images you want to pool',{},root);
    if select_ok == 0
        warning('problem in file selection');
        return;
    end
elseif sum_or_cunj == 3 % automatically selects the AAL atlas to avoid this useless manual step
    AAL_folder = [root, filesep, 'AAL'];
    raw_imgs = [AAL_folder, filesep, 'aal', filesep, 'aal2.nii']; % path to AAL atlas (change according to your own organization)
end
raw_imgs = cellstr(raw_imgs);

%% 3) creates the mask + asks you a name and a path for it

% ask path where to store resulting file
if sum_or_cunj ~= 3
    saveFolder = spm_select(1,'dir','Select the folder where you want to store resulting file',{},root);
elseif sum_or_cunj == 3
    saveFolder = [AAL_folder, filesep, 'ROI_masks', filesep];
end

% ask name for resulting file
fileName = spm_input('file name?',1,'s');


% use ImCalc for conjunction
conjBatch_idx = 1;
matlabbatch{conjBatch_idx}.spm.util.imcalc.input = raw_imgs;
matlabbatch{conjBatch_idx}.spm.util.imcalc.output = fileName;
matlabbatch{conjBatch_idx}.spm.util.imcalc.outdir = {saveFolder};
% spm calculation of the conjunction
if sum_or_cunj == 1
    % matlabbatch{conjBatch_idx}.spm.util.imcalc.expression = 'i1+i2';
    % matlabbatch{conjBatch_idx}.spm.util.imcalc.options.dmtx = 0;
    matlabbatch{conjBatch_idx}.spm.util.imcalc.expression = 'sum(X)';
    matlabbatch{conjBatch_idx}.spm.util.imcalc.options.dmtx = 1;
elseif sum_or_cunj == 2
    matlabbatch{conjBatch_idx}.spm.util.imcalc.expression = '(i2>0).*(i1>0)'; %%% BE CAREFUL: depending on the images you use and the order you enter them, the value may have an impact.
    % Simpler version leave at such (any value higher than zero included.)
    % problem = some of the images are probabilistic => many values
    matlabbatch{conjBatch_idx}.spm.util.imcalc.options.dmtx = 0;
elseif sum_or_cunj == 3
    matlabbatch{conjBatch_idx}.spm.util.imcalc.expression = ['(i1>',num2str(min(AAL_area_id_nber)-1),').*(i1<',num2str(max(AAL_area_id_nber)+1),')']; % keep only the voxels of this particular area
    matlabbatch{conjBatch_idx}.spm.util.imcalc.options.dmtx = 0;
end
matlabbatch{conjBatch_idx}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{conjBatch_idx}.spm.util.imcalc.options.mask = 0; % treats all zeros as NaN values if set to 1
matlabbatch{conjBatch_idx}.spm.util.imcalc.options.interp = 1;
matlabbatch{conjBatch_idx}.spm.util.imcalc.options.dtype = 4;

spm_jobman('run',matlabbatch);
% spm_jobman('interactive',matlabbatch);

close all;

end