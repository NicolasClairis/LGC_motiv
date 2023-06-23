function spm_biascorrect(epi_path)
% this function does automatic bias correction for an EPI series by
% segmenting the first image and applying the bias field to all images
% FORMAT spm_biascorrect(epi) where epi can be a path containing an EPI
% series, or a cellstr containing the EPI images in question
% Dominik R Bach & Guillaume Flandin 3.3.2010

% Disclaimer - This routine has not been fully validated or published.
% _________________________________________________________________________
% Customization by A. Lutti, 2023.
% Copyright (C) 2023 Laboratory for Neuroimaging Research
% Laboratory for Neuroimaging Research, Lausanne University Hospital, Switzerland

% check input
% -------------------------------------------------------------------------
if nargin < 1
    warning('No path specified');
    return;
end

% get EPI images
% -------------------------------------------------------------------------
if ischar(epi_path) && exist(epi_path, 'dir')
%     % Realigned & Unwarped images:
    epi_all = cellstr(spm_select('FPList',epi_path,'^AC.*\.(img|nii)$'));
%     % Realigned images:
%     epi_all = cellstr(spm_select('FPList',epi,'^rf.*\.(img|nii)$'));
%     % Unprocessed images:
%     epi_all = cellstr(spm_select('FPList',epi,'^f.*\.(img|nii)$'));
    epipath = epi_path;
elseif iscell(epi_path)
    epi_all = epi_path;
    [epipath, foo, bar] = fileparts(epi_path{1});
else
    warning('Unknown input arguments'); return;
end
    

clear matlabbatch
matlabbatch{1}.spm.tools.preproc8.channel.vols = {epi_all{1}};
matlabbatch{1}.spm.tools.preproc8.channel.write = [1 0];
spm_jobman('interactive', matlabbatch);


%% read bias field into memory
% -------------------------------------------------------------------------
bf = fullfile(epipath, spm_select('List', epipath, '^BiasField.*\.nii$'));
bfv = spm_vol(bf);
BF = double(spm_read_vols(bfv));


%% apply bias field
% -------------------------------------------------------------------------
for f = 1:numel(epi_all)
    % read file
    fn = [epi_all{f}];
    V = spm_vol(fn);
    Y = spm_read_vols(V);
    % apply bias field
    Y = BF.*Y;
    % save file
    [pth, fnn, ext] = fileparts(fn);
    nfn = fullfile(pth, ['b', fnn, ext]);
    V.fname = nfn;
    spm_write_vol(V,Y);
end
    
return;

