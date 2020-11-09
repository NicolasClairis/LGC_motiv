function [ pval_X_signif_clusters, RFTbetas ] = RFT_extraction( data_to_test, RFT_prm )
%[ pval_X_signif_clusters, RFTbetas ] = RFT_extraction( data_to_test, RFT_prm )
% RFT_extraction performs the Random Field Test analysis (based on the VBA
% toolbox RFT_GLM_contrast function) on a continuous (and preferentially
% smoothed) 1D-data data_to_test. It outputs the index for the significant 
% p.values based on the corrected threshold you enter in pValCorr_threshold
% in the pval_X_signif_clusters output.
%
% RFTbetas contains the whole resume of the RFT_GLM_contrast.m output.
%
% INPUTS
% data_to_test: structure containing the 1D-data on which you want to perform
% the test inside data_to_test.(subfield_nm). The 'subfield_nm' is
% indicated inside RFT_prm.subfield_nm.
%
% RFT_prm contains the following subfields:
% {'subfield_nm','Xvector','con','conType','time_limits','pValCorr_threshold'}
% - RFT_prm.subfield_nm indicates which field of data_to_test you will test
% - RFTprm.Xvector indicates how the data is structured
% - RFT_prm.con and RFT_prm.conType which contrast to perform
% - RFT_prm.time_limits which subset of data_to_test.(current_beta_nm) to
% test
% - RFT_prm.pValCorr_threshold: p.value threshold to use for corrected data
% extraction.
%
% See also RFT_GLM_contrast.m


%% check if VBA toolbox is in the path already or not
% if the VBA toolbox is not installed yet, installs it so that the RFT can
% be made
wherever = pwd;
if exist('RFT_GLM_contrast.m','file') == 0
    VBA_folder = ['C:', filesep, 'Users', filesep, 'nicolas.clairis', filesep, 'Documents', filesep, 'GitHub', filesep, 'VBA-toolbox', filesep];
    cd(VBA_folder);
    VBA_setup;
    cd(wherever);
end

%% main RFT parameters
% X vector
if isfield(RFT_prm,'Xvector') && ~isempty(RFT_prm.Xvector)
    Xvector = RFT_prm.Xvector;
else
    error('Define RFT_prm.Xvector for the current test');
end

% which subfield of data_to_test (and which name to give to the output)
if isfield(RFT_prm,'subfield_nm') && ~isempty(RFT_prm.subfield_nm)
    subfield_nm = RFT_prm.subfield_nm;
else
    error(' Please enter data_to_test as a structure and define RFT_prm.subfield_nm');
end

% which contrast?
if isfield(RFT_prm,'con') && ~isempty(RFT_prm.con) &&...
        isfield(RFT_prm,'conType') && ~isempty(RFT_prm.conType)
    con     = RFT_prm.con;
    conType = RFT_prm.conType;
else % by default F test
    con         = 1; % could be 1 or -1
    conType     = 'F'; % ttest 't' or ftest 'F' (F.test allows to check both positive and negative correlations)
end

% display output or not?
if isfield(RFT_prm,'verbose') && ~isempty(RFT_prm.verbose)
    verbose = RFT_prm.verbose;
else % by default just not show the resulting graphs
    verbose = 0; % display (1) or not (0) the result in a separate figure
end

% test threshold
if isfield(RFT_prm,'pValCorr_threshold') && ~isempty(RFT_prm.pValCorr_threshold)
    pValCorr_threshold = RFT_prm.pValCorr_threshold;
else
    pValCorr_threshold = 0.05; % by default p.value at 0.05 for corrected data
end

%% perform the RFT
% check relative position of columns and lines ok
if size( Xvector, 1) == size( data_to_test.(subfield_nm), 1)
    Yvector = data_to_test.(subfield_nm);
elseif size( Xvector, 1) ~= size( data_to_test.(subfield_nm), 1)
    Yvector = data_to_test.(subfield_nm)';
    warning(['data_to_test.' subfield_nm,' had to be flipped to make the script work. Please check if it''s ok like that.']);
else
    error(['Problem Xvector size does not match with data_to_test.' subfield_nm,' size']);
end

% check only on a specific time limit or on all the data?
if isfield(RFT_prm,'time_limits') && ~isempty(RFT_prm.time_limits)
    time_limits = RFT_prm.time_limits;
else % by default check on the whole data
    time_limits = 1:size(Yvector, 2);
end

[~, RFTbetas.(subfield_nm)] = RFT_GLM_contrast(Xvector, Yvector(:, time_limits), con, conType, [], verbose);

%% extract time points where RFT says there is a significant cluster

nClusters = size(RFTbetas.(subfield_nm).clusters.prft, 2); % extracts number of temporal clusters

if nClusters == 0 % no significant cluster at all
    pval_X_signif_clusters.(subfield_nm) = [];
    
else % if there is at least one significant cluster, extract all the indexes for the statistically significant values
    tmp_index = []; % temporary storage of the indexes
    for iCluster = 1:nClusters % loop through temporal clusters
        if RFTbetas.(subfield_nm).clusters.prft(iCluster) < pValCorr_threshold
            tmp_index = [tmp_index; RFTbetas.(subfield_nm).clusters.ind{iCluster,1}]; % pool all indexes in one variable
        end
    end
    pval_X_signif_clusters.(subfield_nm) = tmp_index;
end


end

