function[] = extract_IP(subject_id, NS)
[mean_deltaIP_perSub_Ep, mean_deltaIP_perSub_Em] = deal(NaN(1,NS));

for iS = 1:NS
    sub_nm = ['CID',subject_id{iS}];
    subFolder = [dataRoot, sub_nm, filesep, 'behavior', filesep];
    IPdata = getfield(load([subFolder,'delta_IP_',sub_nm,'.mat'],'IP_variables'),'IP_variables');
    mean_deltaIP_perSub_Ep(iS) = IPdata.physicalDeltaIP;
    mean_deltaIP_perSub_Em(iS) = IPdata.mentalDeltaIP;
end % subject loop