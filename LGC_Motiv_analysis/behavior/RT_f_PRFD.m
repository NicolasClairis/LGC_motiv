function[] = RT_f_PRFD(figDisp, n_bins)
% RT_f_PRFD(figDisp, n_bins)
% RT_f_PRFD tests whether there is a link between the score in the social
% dominance questionnaire (PRF-D) and reaction times across subjects

%% by default display figure
if ~exist('figDisp','var') || isempty(figDisp) || ~islogical(figDisp)
    figDisp = true;
end

%% define default number of bins
if ~exist(n_bins,'var') || isempty(n_bins) || ~isnumeric(n_bins)
    n_bins = 6;
end

%% subject selection
condition = subject_condition();
gender = 'all';
[subject_id, NS] = LGCM_subject_selection('study1',condition,'all');

%% define main variables
mRT = NaN(1,NS);
nInc = 6;
RT_f_inc = NaN(nInc, NS);

%% load questionnaire results

%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};

end % subject loop

%% average across subjects
mean_sem_sd

%% figure display
if figDisp == true

end % figure display

end % function