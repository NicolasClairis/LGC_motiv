function[ISCE_study_rank, years_studies] = subject_studies_level()
% [ISCE_study_rank, years_studies] = subject_studies_level()
% subject_studies_level will extract the rank and number of years of the
% subjects selected based on the SES questionnaire.
%
% OUTPUTS
% ISCE_study_rank: average rank of study based on the International
% Standard Classification of Education (ISCE)
%
% years_studies: average years of studies based on the category of studies
% (Bachelor, MD, PhD, etc.)

%% subject selection
[study_nm, ~, ~, subject_id, NS] = sub_id;

%% working directory
gitFolder = fullfile('C:','Users','clairis','Desktop');
dataFolder = [fullfile(gitFolder, 'GitHub','LGC_motiv','LGC_Motiv_results',study_nm),filesep];

%% initiate variables of interest
[ISCE_study_rank.allSubs, years_studies.allSubs] = deal(NaN(1,NS));

%% load data
SES_table = readtable([dataFolder,'SES_from_Coline.xlsx']);
for iS = 1:NS
    sub_nm = subject_id{iS};
    sub_idx = find(SES_table.CID == str2double(sub_nm));
    if size(sub_idx, 2) == 1
        ISCE_study_rank.allSubs(iS) = SES_table.niveauD__ducation_bas_SurISCE_(sub_idx);
        years_studies.allSubs(iS) = SES_table.nombreD_ann_esD__tudesApproximatif(sub_idx);
    elseif size(sub_idx,2) == 0
        error(['Missing subject ',sub_nm,' for some reason']);
    elseif size(sub_idx,2) > 1
        error(['Too many occurrences of subject ',sub_nm,' for some reason']);
    end
end % subject loop

%% extract mean/SEM/SD
[ISCE_study_rank.mean,...
    ISCE_study_rank.sem,...
    ISCE_study_rank.sd] = mean_sem_sd(ISCE_study_rank.allSubs, 2);
[years_studies.mean,...
    years_studies.sem,...
    years_studies.sd] = mean_sem_sd(years_studies.allSubs, 2);

end % function