% test canonical correlation analysis (CCA) between plasma metabolites and brain
% metabolites

%% subject selection
[study_nm, condition, subject_id, NS, genderFilter] = subject_selection;

%% load blood metabolites and pool whole-blood + plasma metabolites
% load whole-blood and plasma metabolites
[wholeBlood_mb] = load_blood_NAD(study_nm, subject_id);
[plasmaM, plasma_mb_names, n_plasma_mb] = load_plasma_metabolites(subject_id);
% pool all data together
wholeB_mb_names = fieldnames(wholeBlood_mb);
n_wb_mb = length(wholeB_mb_names);
bloodMb_names = [wholeB_mb_names; plasma_mb_names];
blood_mb_all_possible_vars = [];
for iWB = 1:n_wb_mb
    wb_mb_nm = wholeB_mb_names{iWB};
    blood_mb_all_possible_vars = [blood_mb_all_possible_vars,...
        wholeBlood_mb.(wb_mb_nm)'];
end % whole-blood metabolites loop
for iPB = 1:n_plasma_mb
    plasma_mb_nm = plasma_mb_names{iPB};
    blood_mb_all_possible_vars = [blood_mb_all_possible_vars,...
        plasmaM.(plasma_mb_nm)'];
end % plasma metabolites loop
% select metabolites of interest
blood_mb_idx = listdlg('PromptString','Select blood metabolites',...
    'ListString',bloodMb_names);
blood_mb = blood_mb_all_possible_vars(:,blood_mb_idx); % nObservations*nVars

%% load brain metabolites
[~, ~, brainMetabolites] = metabolite_load(subject_id);
brain_MRS_areas = fieldnames(brainMetabolites);
n_brain_areas = length(brain_MRS_areas);

brain_mb_all_possible_vars = [];
brainMb_names = [];
for iBrain_area = 1:n_brain_areas
    brain_area_nm = brain_MRS_areas{iBrain_area};
    brain_mb_names = fieldnames(brainMetabolites.(brain_area_nm));
    n_brain_mb = length(brain_mb_names);
    for iBrain_mb = 1:n_brain_mb
        brain_mb_nm_tmp = brain_mb_names{iBrain_mb};
        brainMb_names = [brainMb_names;...
            {[brain_area_nm,' ',brain_mb_nm_tmp]}];
        brain_mb_all_possible_vars = [brain_mb_all_possible_vars,...
            brainMetabolites.(brain_area_nm).(brain_mb_nm_tmp)'];
    end % loop through brain metabolites
end % brain areas loop
% select metabolites of interest
brain_mb_idx = listdlg('PromptString','Select brain metabolites',...
    'ListString',brainMb_names);
brain_mb = brain_mb_all_possible_vars(:,brain_mb_idx); % nObservations*nVars

%% filter datapoints with NaN values
okSubs = false(NS, 1);
for iS = 1:NS
    if sum(isnan([blood_mb(iS,:), brain_mb(iS,:)])) > 0
        okSubs(iS) = false;
    elseif sum(isnan([blood_mb(iS,:), brain_mb(iS,:)])) == 0
        okSubs(iS) = true;
    end
end % subject loop

%% perform CCA
[A, B, r, U, V, stats] = canoncorr(blood_mb(okSubs,:), brain_mb(okSubs,:));
% A: canonical coefficients of blood_mb
% B: canonical coefficients for brain_mb
% r: canonical correlations (one for each canonical variable)
% U & V: canonical variables for blood_mb and brain_mb
% stats: statistical information

%% display the results
