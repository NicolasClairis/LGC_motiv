function[mean_deltaIP_perSub_Ep, mean_deltaIP_perSub_Em] = extract_IP(subject_id, NS)
% [mean_deltaIP_perSub_Ep, mean_deltaIP_perSub_Em] = extract_IP(subject_id, NS)
% extract_IP will extract the indifference point (IP) for physical and
% mental effort
%
% INPUTS
% subject_id: list of subjects
%
% NS: number of subjects included
%
% OUTPUTS
% mean_deltaIP_perSub_Ep: mean delta for indifference point for physical
% effort
%
% mean_deltaIP_perSub_Em: mean delta for indifference point for mental
% effort

%% working directory
study_nm = 'study1';
dataRoot = [fullfile('E:',study_nm),filesep];

%% initialize variables
[mean_deltaIP_perSub_Ep, mean_deltaIP_perSub_Em] = deal(NaN(1,NS));

% loop over subjects
for iS = 1:NS
    sub_nm = ['CID',subject_id{iS}];
    subFolder = [dataRoot, sub_nm, filesep, 'behavior', filesep];
    IPdata = getfield(load([subFolder,'delta_IP_',sub_nm,'.mat'],'IP_variables'),'IP_variables');
    mean_deltaIP_perSub_Ep(iS) = IPdata.physicalDeltaIP;
    mean_deltaIP_perSub_Em(iS) = IPdata.mentalDeltaIP;
end % subject loop

end % function