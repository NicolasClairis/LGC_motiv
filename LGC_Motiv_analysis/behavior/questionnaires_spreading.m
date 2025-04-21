function[NS] = questionnaires_spreading()
% [NS] = questionnaires_spreading()
% questionnaires_spreading aims at displaying how the selected
% questionnaires spread within the selected population.
%
% OUTPUTS
% NS: number of subjects

%% subject selection
[study_nm, condition, subject_id, NS, genderFilter] = subject_selection;

%% load questionnaires
[questionnaires, categ_quests, n_categ, subject_id] = extract_questionnaires(study_nm, subject_id, NS);

%% select questionnaires of interest
% 1) select category of questionnaires of interest
selectedCateg_idx = listdlg('PromptString','Select category of the questionnaires to look',...
    'ListString',categ_quests);
categ_nm = categ_quests{selectedCateg_idx};

% 2) select questionnaires of interest within the selected category
questionnaire_list = fieldnames(questionnaires.(categ_nm));
quest_idx = listdlg('PromptString','Select category of the questionnaires to look',...
    'ListString',questionnaire_list);
n_selected_Q = length(quest_idx);

%% display spreading of values for each relevant questionnaire
dHist = 0.52; % space to add to visualize the borders in the histogram
for iQ = 1:n_selected_Q
    quest_nm = questionnaire_list{quest_idx(iQ)};
    quest_data = questionnaires.(categ_nm).(quest_nm);
    fig;
    histogram(quest_data);
    xlabel(strrep(quest_nm,'_',' '));
    ylabel('N. subjects');
    % adapt x.scale to highlight spreading of values within the range of
    % each scale
    switch quest_nm
        case {'MPSTEFS_physical_energy',...
                'MPSTEFS_mental_energy',...
                'MPSTEFS_physical_fatigue',...
                'MPSTEFS_mental_fatigue'} % energy and fatigue
            x_vals = [3 15];
        case 'JPIR' % energy
            x_vals = [0 20];
        case 'MADRS_S' % depression
            x_vals = [0 27];
        otherwise
            warning(['Consider updating scale for ',quest_nm,...
                ' to display full range of possible values']);
            x_vals = [min(quest_data), max(quest_data)];
    end
    xlim([x_vals(1)-dHist x_vals(2)+dHist]); % add space for the borders in the histogram
end % questionnaire loop

end % function