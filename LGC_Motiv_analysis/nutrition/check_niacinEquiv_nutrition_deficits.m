% check whether participants had any deficit in niacin through nutrition.
% recommended dietary allowance (RDA) of 14mg/day of NE based on the 
% following document:
% "Dietary Reference Intakes for Thiamin, Riboflavin,
% Niacin, Vitamin B6, Folate, Vitamin B12, Pantothenic
% Acid, Biotin, and Choline"
% https://www.ncbi.nlm.nih.gov/books/NBK114310/pdf/Bookshelf_NBK114310.pdf

%% working directories
root = LGCM_root_paths;
% study
study_nm = 'study1';
studyPath = [root, filesep, study_nm, filesep];
% nutrition data
switch root
    case ['E:',filesep]
        gitPath = fullfile('C:','Users','clairis','Desktop');
    case {[fullfile('C:','Users','Loco','Downloads'),filesep],...
            [fullfile('L:','human_data_private','raw_data_subject'),filesep]}
        gitPath = fullfile('C:','Users','Loco','Downloads');
    otherwise
        error('case not ready yet');
end
nutritionPath = [fullfile(gitPath,'GitHub','LGC_motiv',...
    'LGC_Motiv_results',study_nm,'nutrition'),filesep];

%% subject definition
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% extract niacin equivalents data
niacineEquivFilePath = [nutritionPath,'niacine_equivalents_scoring.xlsx'];
niacinEquiv_table = readtable(niacineEquivFilePath,...
    'Sheet','Sheet1');

% niacin equivalents (in ug/week)
nutri.niacinEquiv = NaN(1,NS);

%% loop through subjects to extract nutrition
for iS = 1:NS
    sub_nm = subject_id{iS};
    
    %% load nutrition score
    sub_niacinEquiv_idx = find(strcmp(niacinEquiv_table.CID, sub_nm));
    % extract nutrition intake values
    % extract niacin equivalents values
    if ~isempty(sub_niacinEquiv_idx) && size(sub_niacinEquiv_idx,1) == 1
        nutri.niacinEquiv(iS) = niacinEquiv_table.x_quivalentsDeNiacineParSemaine__g_(sub_niacinEquiv_idx);
    end
end % subject loop

%% convert data from (μg/week) to (mg/day)
NE = (nutri.niacinEquiv./7)./1000;

%% recommended daily intake of niacin equivalents threshold
NE_min_RDA = 14;

% extract info about subjects with low NE
low_NE_subs = NE <= NE_min_RDA;
n_low_NE_subs = sum(low_NE_subs);

%% display graph
figure;
mHdl = scatter(1:NS, sort(NE));
mHdl.MarkerEdgeColor = 'k';
mHdl.LineWidth = 1;
hold on;
line(xlim(),[NE_min_RDA, NE_min_RDA],...
    'Color','k','LineWidth',1);
% mark subjects with deficit

mHdl = scatter(1:n_low_NE_subs, sort(NE(low_NE_subs)));
mHdl.MarkerEdgeColor = 'r';
mHdl.LineWidth = 1;
ylabel('Niacin equivalents (mg/day)');
legend_size(30);