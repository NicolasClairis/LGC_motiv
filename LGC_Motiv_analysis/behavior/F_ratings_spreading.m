function[NS] = F_ratings_spreading()
% [NS] = F_ratings_spreading()
% F_ratings_spreading aims at displaying how the fatigue ratings
% spread within the selected population.
%
% OUTPUTS
% NS: number of subjects

%% subject selection
[study_nm, condition, subject_id, NS, genderFilter] = subject_selection;

%% load fatigue ratings
[~,...
    F1, F2,...
    F3, F4,...
    F4_min_F1, F4_min_F3] = extract_subjective_fatigue_ratings(study_nm, subject_id, NS);
fatigue.F1          = F1;
fatigue.F2          = F2;
fatigue.F3          = F3;
fatigue.F4          = F4;
fatigue.F4_min_F1   = F4_min_F1;
fatigue.F4_min_F3   = F4_min_F3;
fatigue_names = fieldnames(fatigue);
n_fatigue_measures = length(fatigue_names);

%% display spreading of values for each relevant questionnaire
dHist = 0.52; % space to add to visualize the borders in the histogram
for iF = 1:n_fatigue_measures
    F_nm = fatigue_names{iF};
    
    fig;
    histogram(fatigue.(F_nm));
    ylabel('N. subjects');
    % adapt x.scale to highlight spreading of values within the range of
    % each scale
    switch F_nm
        case {'F1','F2','F3','F4'}
            x_vals = [0 10];
            xlabel(F_nm);
        case {'F4_min_F1','F4_min_F3'}
            x_vals = [-10 10];
            xlabel(strrep(F_nm,'_min_','-'));
    end
    xlim([x_vals(1)-dHist x_vals(2)+dHist]); % add space for the borders in the histogram
end % questionnaire loop

end % function