%% check TESTO, CORT and TESTO/CORT distribution per gender

%% general parameters
study_nm = 'study1';
condition = 'fullList';

%% load data
[TESTO_data] = load_TESTO(study_nm);
[CORT_data] = load_CORT(study_nm);
TESTO_id = strrep(TESTO_data.CID,'CID','');
CORT_id = strrep(CORT_data.CID,'CID','');
NS_TESTO = size(TESTO_id,2);
NS_CORT = size(CORT_id,2);

%% load subjects of interest
males_id = LGCM_subject_selection(study_nm,condition,'males');
females_id = LGCM_subject_selection(study_nm,condition,'females');

%% extract relevant data
% extract gender index
[males_TESTO_id, females_TESTO_id] = deal(false(1,NS_TESTO));
for iS = 1:NS_TESTO
    if ismember(TESTO_id{iS},males_id)
        males_TESTO_id(iS) = true;
    elseif ismember(TESTO_id{iS},females_id)
        females_TESTO_id(iS) = true;
    end
end % subject loop
[males_CORT_id, females_CORT_id] = deal(false(1,NS_CORT));
for iS = 1:NS_CORT
    if ismember(CORT_id{iS},males_id)
        males_CORT_id(iS) = true;
    elseif ismember(CORT_id{iS},females_id)
        females_CORT_id(iS) = true;
    end
end % subject loop
% extract data for first sample
males_TESTO_A1 = TESTO_data.TESTO(1,males_TESTO_id);
females_TESTO_A1 = TESTO_data.TESTO(1,females_TESTO_id);
males_CORT_A1 = CORT_data.CORT(1,males_CORT_id);
females_CORT_A1 = CORT_data.CORT(1,females_CORT_id);
% compute ratio (after normalizing both in same units (μg/dL))
males_TESTO_CORT_A1 = (males_TESTO_A1.*0.0001)./males_CORT_A1;
females_TESTO_CORT_A1 = (females_TESTO_A1.*0.0001)./females_CORT_A1;

% extract data for last sample
males_TESTO_E = TESTO_data.TESTO(4,males_TESTO_id);
females_TESTO_E = TESTO_data.TESTO(4,females_TESTO_id);
males_CORT_E = CORT_data.CORT(4,males_CORT_id);
females_CORT_E = CORT_data.CORT(4,females_CORT_id);
% compute ratio (after normalizing both in same units (μg/dL))
males_TESTO_CORT_E = (males_TESTO_E.*0.0001)./males_CORT_E;
females_TESTO_CORT_E = (females_TESTO_E.*0.0001)./females_CORT_E;

%% display data
pSize = 30;
%% for first sample
fig_A1 = fig;

% testosterone
subplot(1,3,1);
hold on;
Violin({males_TESTO_A1}, 1);
Violin({females_TESTO_A1},2);
xticks(1:2);
xticklabels({'males','females'});
ylabel('Testosterone (pg/mL)');
legend_size(pSize);

% cortisol
subplot(1,3,2);
Violin({males_CORT_A1},1);
Violin({females_CORT_A1},2);
xticks(1:2);
xticklabels({'males','females'});
ylabel('Cortisol (μg/dL)');
legend_size(pSize);

% testosterone/cortisol
subplot(1,3,3);
Violin({males_TESTO_CORT_A1},1);
Violin({females_TESTO_CORT_A1},2);
xticks(1:2);
xticklabels({'males','females'});
ylabel('Testosterone/Cortisol');
legend_size(pSize);

%% for last sample
fig_E = fig;

% testosterone
subplot(1,3,1);
hold on;
Violin({males_TESTO_E}, 1);
Violin({females_TESTO_E},2);
xticks(1:2);
xticklabels({'males','females'});
ylabel('Testosterone (pg/mL)');
legend_size(pSize);

% cortisol
subplot(1,3,2);
Violin({males_CORT_E},1);
Violin({females_CORT_E},2);
xticks(1:2);
xticklabels({'males','females'});
ylabel('Cortisol (μg/dL)');
legend_size(pSize);

% testosterone/cortisol
subplot(1,3,3);
Violin({males_TESTO_CORT_E},1);
Violin({females_TESTO_CORT_E},2);
xticks(1:2);
xticklabels({'males','females'});
ylabel('Testosterone/Cortisol');
legend_size(pSize);