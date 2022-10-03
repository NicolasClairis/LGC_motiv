function[betas, pval] = brain_GSH_f_nutrition_precursors(figDisp)
% [betas, pval] = brain_GSH_f_nutrition_precursors(figDisp)
%% brain_GSH_f_nutrition_precursors will test whether there is any 
% correlation between nutrition precursors of GSH (namely Gly, Cys and Glu)
% and brain levels of GSH in the dmPFC or aINS (for study 1)
%
% INPUTS
% figDisp: display figure (1) or not (0) ? Will be equal to 1 by default
% if not enterrd
%
% OUTPUTS
% betas: structure with betas for each test
%
% pval: structure with p.value for each test

%% main parameters
if ~exist('figDisp','var') || isempty(figDisp)
    figDisp = 1;
end
%% working directories
root = LGCM_root_paths;
% study
study_nm = 'study1';
studyPath = [root, filesep, study_nm, filesep];

% condition
condition = subject_condition;
[subject_id, NS] = LGCM_subject_selection(study_nm, condition);

%% initialize variables of interest
[dmPFC_GSH, aINS_GSH,...
    nutri.Gly, nutri.Cys, nutri.Glu,...
    nutri.Gly_div_totalCal,  nutri.Cys_div_totalCal, nutri.Glu_div_totalCal,...
    nutri.GlyCysGlu_sum,...
    nutri.GlyCysGlu_sum_div_totalCal] = deal(NaN());

%% loop through subjects
for iS = 1:NS
    sub_nm = subject_id{iS};
    subPath = [studyPath, sub_nm];
    
    %% load brain metabolites
    [metabolite_tmp] = metabolite_load({sub_nm});
    dmPFC_GSH(iS) = metabolite_tmp.dmPFC.GSH;
    aINS_GSH(iS) = metabolite_tmp.aIns.GSH;
    
    %% load nutrition score
    
end % subject loop

%% average

%% figure
if figDisp == 1
    % show result Glu/Cys/Gly alone
    fig;
    % Glutamate
    subplot(2,3,1);
    % Glutamate/total calories
    subplot(2,3,4);
    
    % Cysteine
    subplot(2,3,2);
    % Cysteine/total calories
    subplot(2,3,5);
    
    % Glycine
    subplot(2,3,3);
    % Glycine/total calories
    subplot(2,3,6);
    
end % figure display
end % function