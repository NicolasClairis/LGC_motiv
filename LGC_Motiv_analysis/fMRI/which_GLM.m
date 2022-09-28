function [GLMprm] = which_GLM(GLM)
% [GLMprm] = which_GLM(GLM)
% which_GLM defines GLM parameters (GLMprm) depending on the GLM number
% defined in GLM.
%
% INPUTs
% GLM: GLM number
%
% OUTPUTS
% GLMprm: GLM parameters = structure with all relevant informations for 1st
% and 2nd level analysis
%   .gal: general parameters
%       .add_drv: model spatial and/or a temporal derivative of the BOLD
%       (0) no derivative added
%       (1) temporal derivative
%       (2) temporal AND spatial derivative
%
%       .grey_mask:
%       (0) include all voxels
%       (1) grey-matter filter based on individual grey matter
%       (2) use SPM template grey matter mask
%       (3) grey-matter filter based on group average grey-matter
%
%       .orth_vars:
%       (0) don't orthogonalize the regressors
%       (1) orthogonalize regressors of the GLM (1)
%
%       .zPerRun:
%       (0) raw values (or whatever is defined for each regressor)
%       (1) values zscored per run
%
%       .onsets_only:
%       (0) regular GLM (by default)
%       (1) extract 1 beta/trial/condition
%
%   .model_onset: indicate for each task (Ep/Em: physical/mental) for each
%   event (preChoiceCross/choice/chosen/preEffortCross/Eperf/fbk) if it should be modelled as a
%   stick ('stick') as a boxcar ('boxcar') or not included in the GLM
%   ('none'), for some cases, 'boxcar_bis' corresponds to situation where
%   the boxcar also entails following periods of the task
%       .preChoiceCross: (white) fixation cross before choice period
%       .choice: choice period (when options are displayed on screen)
%       .chosen: moment when the chosen option is displayed on screen
%       .preEffortCross: (black) fixation cross before effort period
%       .Eperf: physical/mental effort performance period
%       .fbk: feedback period
%
%   .preChoiceCross regressors:
%       .preChoiceCross.(Ep/Em).RT
%       (1) raw reaction time
%       (2) reaction time zscored per run
%       (3) reaction time zscored per subject across all runs
%       (4-6) same as (1-3) but with RT as first regressor instead of last
%
%       .preChoiceCross.(Ep/Em).choiceHighE
%       (1) 0 when low effort/default option is selected and 1 when high
%       effort/non-default option is selected
%
%   .choice/chosen/preEffortCross/Eperf/fbk: for each event, for each task (Ep/Em) and 
%   for each condition (R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3) indicate if a given regressor should be 
%   included or not.
%   .Ep/Em: physical (Ep) or mental (Em) effort task
%   .R/P/RP/splitPerE/splitPerEch: reward only (R), punishment only (P) or reward and punishment
%   trials mixed (RP) or split per effort proposed (splitPerE) or per
%   effort chosen (splitPerEch)
%
%       .(choice/chosen).(Ep/Em).RPpool:
%       (0) split rewards and punishments as separate events
%       (1) pool reward and punishment trials
%
%       .(choice/chosen).(Ep/Em).splitPerE:
%       (0) pool all effort together
%       (1) split per level of effort proposed (E1/E2/E3)
%       (2) split per level of effort chosen (Ech0/Ech1/Ech2/Ech3)
%       (3) split per option chosen (low/default vs high/non-default effort chosen)
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).R_vs_P:
%       (1) 0 when punishment trial, 1 when reward trial
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).choiceHighE
%       (1) 0 when low effort/default option is selected and 1 when high
%       effort/non-default option is selected
%
%       .choice.(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).money_left
%       (1) money amount associated to left option
%
%       .choice.(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).money_right
%       (1) money amount associated to right option
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).money_chosen:
%       (1) money chosen amount
%       (2) |money chosen amount|
%       (3) money levels (-4/-3/-2/-1/1/2/3/4) chosen from highest loss
%       until highest gain
%       (4) |money levels (1/2/3/4) chosen|
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).money_unchosen:
%       (1) money unchosen amount
%       (2) |money unchosen amount|
%       (3) money levels (-4/-3/-2/-1/1/2/3/4) unchosen
%       (4) |money levels (1/2/3/4) unchosen|
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).money_varOption:
%       (1) money amount non-default option
%       (2) |money amount non-default option|
%       (3) money levels  (-3/-2/-1/1/2/3) non-default option
%       (4) |money levels (1/2/3) non-default option|
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).money_ch_min_unch:
%       (1) money chosen - money unchosen amount
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).money_ch_min_fixOption:
%       (1) money chosen - money default option - amount
%       (2) money chosen - money default option - levels (1/2/3/4 - 1)
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).money_sum:
%       (1) money default + money non-default amount
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).E_left
%       (1) effort level (0/1/2/3) associated to left option
%       (2) effort difficulty associated to left option (Ep: duration to hold; Em: nb answers to give)
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).E_right
%       (1) effort level (0/1/2/3) associated to right option
%       (2) effort difficulty associated to right option (Ep: duration to hold; Em: nb answers to give)
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).E_chosen
%       (1) effort level (0/1/2/3) associated to chosen option
%       (2) effort difficulty associated to chosen option (Ep: duration to hold; Em: nb answers to give)
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).E_unchosen
%       (1) effort level (0/1/2/3) associated to unchosen option
%       (2) effort difficulty associated to unchosen option (Ep: duration to hold; Em: nb answers to give)
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).E_varOption
%       (1) effort level (0/1/2/3) associated to the non-default option
%       (2) effort difficulty associated to the non-default option (Ep: duration to hold; Em: nb answers to give)
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).E_ch_min_unch
%       (1) effort level (0/1/2/3) associated to the chosen option minus
%       the effort level (0/1/2/3) associated to the unchosen option
%       (2) effort difficulty associated to the chosen option (Ep: duration to hold; Em: nb answers to give)
%       minus the effort difficulty associated to the unchosen option (Ep: duration to hold; Em: nb answers to give)
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).E_ch_min_fixOption:
%       (1) E chosen - E default option - levels (0/1/2/3)
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).E_sum
%       (1) sum of the effort levels (0/1/2/3) associated to both options
%       (2) sum of the effort difficulties (Ep: duration to hold; Em: nb answers to give) associated to both options
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).NV_chosen
%       (1) net value of the chosen option based on the model defined in
%       .(choice/chosen).(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).NV_mdl (='mdl_X' or 'bayesianModel_X')
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).NV_varOption
%       (1) net value of the chosen option based on the model defined in
%       .(choice/chosen).(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).NV_mdl (='mdl_X' or 'bayesianModel_X')
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).trialN
%       (1) trial number
%       (2) (trial number)*(E chosen - E non-chosen option)
%       (3) (trial number)*(E non-default - E default option)
%       (4) (trial number)*(E non-default)
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).confidence
%       (1) confidence level (0/1) given by the subject for each choice
%       (2) confidence inferred by the model (p(choice)-0.5)² for the model
%       defined in .(choice/chosen).(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).conf_mdl (='mdl_X' or 'bayesianModel_X')
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).RT: reaction time for choice
%       (1) raw reaction time
%       (2) reaction time zscored per run
%       (3) reaction time zscored per subject across all runs
%       (4-6) same as (1-3) but with RT as first regressor instead of last
%
%
%       .Eperf.(Ep/Em).RPpool
%       (0) split rewards and punishments as separate events
%       (1) pool reward and punishment trials
%
%       .Eperf.(Ep/Em).splitPerE:
%       (0) pool all effort together
%       (1) split per level of effort proposed (E1/E2/E3)
%       (2) split per level of effort chosen (Ech0/Ech1/Ech2/Ech3)
%       (3) split per option chosen (low/default vs high/non-default effort chosen)
%
%       .Eperf.(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).choiceHighE
%       (1) 0 when low effort/default option is selected and 1 when high
%       effort/non-default option is selected
%
%       .Eperf.(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).money_chosen
%       (1) money chosen amount
%       (2) |money chosen amount|
%       (3) money levels (-4/-3/-2/-1/1/2/3/4) chosen
%       (4) |money levels (1/2/3/4) chosen|
%
%       .Eperf.(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).E_chosen
%       (1) effort level (0/1/2/3) associated to chosen option
%       (2) effort difficulty associated to chosen option (Ep: duration to hold; Em: nb answers to give)
%
%       .Eperf.Ep.(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).F_peak: physical effort only: force peak
%       (1) force peak in newtons
%
%       .Eperf.Ep.(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).F_integral: physical effort only: force integral
%       (1) integral of effort performed during the effort period (sum of 
%       efforts in newtons)
%
%       .Eperf.Em.(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).RT_avg: mental effort only: average reaction
%       time for answering N-back task
%       (1) average reaction time for answering to questions (in seconds)
%
%       .Eperf.Em.(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).n_errors: mental effort only: number of errors
%       made per trial
%       (1) number of errors made
%
%       .Eperf.(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).NV_chosen
%       (1) net value of the chosen option based on the model defined in
%       .Eperf.(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).NV_mdl (='mdl_X' or 'bayesianModel_X')
%
%       .Eperf.(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).NV_varOption
%       (1) net value of the chosen option based on the model defined in
%       .Eperf.(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).NV_mdl (='mdl_X' or 'bayesianModel_X')
%
%       .Eperf.(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).RT_1stAnswer
%       (1) raw reaction time for first answer (force above threshold for Ep 
%       and first answer to first digit for Em)
%
%       .Eperf.(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).trialN
%       (1) trial number
%       (2) (trial number)*(E chosen - E non-chosen option)
%       (3) (trial number)*(E non-default - E default option)
%       (4) (trial number)*(E non-default)
%
%       .Eperf.(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).confidence
%       (1) confidence level (0/1) given by the subject for each choice
%       (2) confidence inferred by the model (p(choice)-0.5)² for the model
%       defined in .Eperf.(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).conf_mdl (='mdl_X' or 'bayesianModel_X')
%
%
%       .fbk.(Ep/Em).RPpool
%       (0) split rewards and punishments as separate events
%       (1) pool reward and punishment trials
%
%       .fbk.(Ep/Em).splitPerE:
%       (0) pool all effort together
%       (1) split per level of effort proposed (E1/E2/E3)
%       (2) split per level of effort chosen (Ech0/Ech1/Ech2/Ech3)
%       (3) split per option chosen (low/default vs high/non-default effort chosen)
%
%       .fbk.(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).win_vs_loss
%       (1) win (1) - loss (0) trials
%
%       .fbk.(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).choiceHighE
%       (1) 0 when low effort/default option is selected and 1 when high
%       effort/non-default option is selected
%
%       .fbk.(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).money_obtained
%       (1) money amount obtained at the end of the trial
%       (1) |money amount| obtained at the end of the trial (~saliency)
%
%       .fbk.(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).E_made
%       (1) level of effort performed during effort performance
%
%       .fbk.(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).trialN
%       (1) trial number
%       (2) (trial number)*(E chosen - E non-chosen option)
%       (3) (trial number)*(E non-default - E default option)
%       (4) (trial number)*(E non-default)
%
%       .fbk.(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).confidence
%       (1) confidence level (0/1) given by the subject for each choice
%       (2) confidence inferred by the model (p(choice)-0.5)² for the model
%       defined in .fbk.(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3).conf_mdl (='mdl_X' or 'bayesianModel_X')
%
% See also GLM_details.m
%
% Nicolas Clairis - 2021/2022


%% initialize all variables

% general parameters: derivative, grey matter filtering
[GLMprm.gal.add_drv,...
    GLMprm.gal.grey_mask,...
    GLMprm.gal.zPerRun,...
    GLMprm.gal.orth_vars,...
    GLMprm.gal.onsets_only] = deal(0);

% onsets: not modelled (none), modelled as stick function (stick) or as
% boxcar function (boxcar)
[GLMprm.model_onset.Ep.preChoiceCross,...
    GLMprm.model_onset.Ep.preEffortCross,...
    GLMprm.model_onset.Ep.choice,...
    GLMprm.model_onset.Ep.chosen,...
    GLMprm.model_onset.Ep.Eperf,...
    GLMprm.model_onset.Ep.fbk,...
    GLMprm.model_onset.Em.preChoiceCross,...
    GLMprm.model_onset.Em.preEffortCross,...
    GLMprm.model_onset.Em.choice,...
    GLMprm.model_onset.Em.chosen,...
    GLMprm.model_onset.Em.Eperf,...
    GLMprm.model_onset.Em.fbk] = deal('none');


% initialize all modulators of interest at 0
%
% use same conditions for physical and mental effort and pool reward and
% punishment trials or split reward and punishment trials
Ep_Em = {'Ep','Em'}; % apply same or different conditions to physical and mental effort
for iEpm = 1:length(Ep_Em)
    EpEm_nm = Ep_Em{iEpm};
    
    % no regressor during pre-choice cross (by default)
    GLMprm.preChoiceCross.(EpEm_nm).RT = 0;
    GLMprm.preChoiceCross.(EpEm_nm).choiceHighE = 0;
    
    % pool reward and punishment together (by default)
    GLMprm.choice.(EpEm_nm).RPpool = 1;
    GLMprm.chosen.(EpEm_nm).RPpool = 1;
    GLMprm.preEffortCross.(EpEm_nm).RPpool = 1;
    GLMprm.Eperf.(EpEm_nm).RPpool = 1;
    GLMprm.fbk.(EpEm_nm).RPpool = 1;
    RPconditions = {'R','P','RP'};
    
    % by default pool all effort trials
    GLMprm.choice.(EpEm_nm).splitPerE = 0;
    GLMprm.chosen.(EpEm_nm).splitPerE = 0;
    GLMprm.preEffortCross.(EpEm_nm).splitPerE = 0;
    GLMprm.Eperf.(EpEm_nm).splitPerE = 0;
    GLMprm.fbk.(EpEm_nm).splitPerE = 0;
    
    Econditions = {'E','E1','E2','E3',...
        'E1','E2','E3',...
        'Ech0','Ech1','Ech2','Ech3',...
        'lEch','hEch'};
    % loop through conditions
    for iRP = 1:length(RPconditions)
        RP_nm = RPconditions{iRP};
        
        for iEsplit = 1:length(Econditions)
            Econd_nm = Econditions{iEsplit};
            
            % by default all regressors are not included
            [GLMprm.choice.(EpEm_nm).(RP_nm).(Econd_nm).R_vs_P,...
                GLMprm.choice.(EpEm_nm).(RP_nm).(Econd_nm).choiceHighE,...
                GLMprm.choice.(EpEm_nm).(RP_nm).(Econd_nm).money_left,...
                GLMprm.choice.(EpEm_nm).(RP_nm).(Econd_nm).money_right,...
                GLMprm.choice.(EpEm_nm).(RP_nm).(Econd_nm).money_chosen,...
                GLMprm.choice.(EpEm_nm).(RP_nm).(Econd_nm).money_unchosen,...
                GLMprm.choice.(EpEm_nm).(RP_nm).(Econd_nm).money_varOption,...
                GLMprm.choice.(EpEm_nm).(RP_nm).(Econd_nm).money_ch_min_unch,...
                GLMprm.choice.(EpEm_nm).(RP_nm).(Econd_nm).money_ch_min_fixOption,...
                GLMprm.choice.(EpEm_nm).(RP_nm).(Econd_nm).money_sum,...
                GLMprm.choice.(EpEm_nm).(RP_nm).(Econd_nm).E_left,...
                GLMprm.choice.(EpEm_nm).(RP_nm).(Econd_nm).E_right,...
                GLMprm.choice.(EpEm_nm).(RP_nm).(Econd_nm).E_chosen,...
                GLMprm.choice.(EpEm_nm).(RP_nm).(Econd_nm).E_unchosen,...
                GLMprm.choice.(EpEm_nm).(RP_nm).(Econd_nm).E_varOption,...
                GLMprm.choice.(EpEm_nm).(RP_nm).(Econd_nm).E_ch_min_unch,...
                GLMprm.choice.(EpEm_nm).(RP_nm).(Econd_nm).E_ch_min_fixOption,...
                GLMprm.choice.(EpEm_nm).(RP_nm).(Econd_nm).E_sum,...
                GLMprm.choice.(EpEm_nm).(RP_nm).(Econd_nm).NV_chosen,...
                GLMprm.choice.(EpEm_nm).(RP_nm).(Econd_nm).NV_varOption,...
                GLMprm.choice.(EpEm_nm).(RP_nm).(Econd_nm).confidence,...
                GLMprm.choice.(EpEm_nm).(RP_nm).(Econd_nm).RT,...
                GLMprm.choice.(EpEm_nm).(RP_nm).(Econd_nm).trialN] = deal(0);
            [GLMprm.choice.(EpEm_nm).(RP_nm).(Econd_nm).NV_mdl,...
                GLMprm.choice.(EpEm_nm).(RP_nm).(Econd_nm).conf_mdl] = deal('');
            
            % chosen option display
            % by default all regressors are not included
            [GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).R_vs_P,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).choiceHighE,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).money_chosen,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).money_unchosen,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).money_varOption,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).money_ch_min_unch,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).money_ch_min_fixOption,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).money_sum,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).E_chosen,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).E_unchosen,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).E_varOption,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).E_ch_min_unch,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).E_ch_min_fixOption,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).E_sum,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).NV_chosen,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).NV_varOption,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).confidence,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).RT,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).trialN] = deal(0);
            [GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).NV_mdl,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).conf_mdl] = deal('');
            
            % pre-effort black cross
            % by default all regressors are not included
            [GLMprm.preEffortCross.(EpEm_nm).(RP_nm).(Econd_nm).choiceHighE,...
                GLMprm.preEffortCross.(EpEm_nm).(RP_nm).(Econd_nm).money_chosen,...
                GLMprm.preEffortCross.(EpEm_nm).(RP_nm).(Econd_nm).E_chosen,...
                GLMprm.preEffortCross.(EpEm_nm).(RP_nm).(Econd_nm).NV_chosen,...
                GLMprm.preEffortCross.(EpEm_nm).(RP_nm).(Econd_nm).NV_varOption,...
                GLMprm.preEffortCross.(EpEm_nm).(RP_nm).(Econd_nm).RT_1stAnswer,...
                GLMprm.preEffortCross.(EpEm_nm).(RP_nm).(Econd_nm).trialN,...
                GLMprm.preEffortCross.(EpEm_nm).(RP_nm).(Econd_nm).confidence] = deal(0);
            % specific variables for each effort type
            switch EpEm_nm
                case 'Ep'
                    % effort performance
                    [GLMprm.preEffortCross.Ep.(RP_nm).(Econd_nm).F_peak,...
                        GLMprm.preEffortCross.Ep.(RP_nm).(Econd_nm).F_integral] = deal(0);
                case 'Em'
                    % effort performance
                    [GLMprm.preEffortCross.Em.(RP_nm).(Econd_nm).RT_avg,...
                        GLMprm.preEffortCross.Em.(RP_nm).(Econd_nm).n_errors] = deal(0);
            end
            [GLMprm.preEffortCross.(EpEm_nm).(RP_nm).(Econd_nm).NV_mdl,...
                GLMprm.preEffortCross.(EpEm_nm).(RP_nm).(Econd_nm).conf_mdl] = deal('');
            
            % effort performance
            % by default all regressors are not included
            [GLMprm.Eperf.(EpEm_nm).(RP_nm).(Econd_nm).choiceHighE,...
                GLMprm.Eperf.(EpEm_nm).(RP_nm).(Econd_nm).money_chosen,...
                GLMprm.Eperf.(EpEm_nm).(RP_nm).(Econd_nm).E_chosen,...
                GLMprm.Eperf.(EpEm_nm).(RP_nm).(Econd_nm).NV_chosen,...
                GLMprm.Eperf.(EpEm_nm).(RP_nm).(Econd_nm).NV_varOption,...
                GLMprm.Eperf.(EpEm_nm).(RP_nm).(Econd_nm).RT_1stAnswer,...
                GLMprm.Eperf.(EpEm_nm).(RP_nm).(Econd_nm).trialN,...
                GLMprm.Eperf.(EpEm_nm).(RP_nm).(Econd_nm).confidence] = deal(0);
            % specific variables for each effort type
            switch EpEm_nm
                case 'Ep'
                    % effort performance
                    [GLMprm.Eperf.Ep.(RP_nm).(Econd_nm).F_peak,...
                        GLMprm.Eperf.Ep.(RP_nm).(Econd_nm).F_integral] = deal(0);
                case 'Em'
                    % effort performance
                    [GLMprm.Eperf.Em.(RP_nm).(Econd_nm).RT_avg,...
                        GLMprm.Eperf.Em.(RP_nm).(Econd_nm).n_errors] = deal(0);
            end
            [GLMprm.Eperf.(EpEm_nm).(RP_nm).(Econd_nm).NV_mdl,...
                GLMprm.Eperf.(EpEm_nm).(RP_nm).(Econd_nm).conf_mdl] = deal('');
            
            % feedback
            % by default all regressors are not included
            [GLMprm.fbk.(EpEm_nm).(RP_nm).(Econd_nm).choiceHighE,...
                GLMprm.fbk.(EpEm_nm).(RP_nm).(Econd_nm).money_obtained,...
                GLMprm.fbk.(EpEm_nm).(RP_nm).(Econd_nm).win_vs_loss,...
                GLMprm.fbk.(EpEm_nm).(RP_nm).(Econd_nm).E_made,...
                GLMprm.fbk.(EpEm_nm).(RP_nm).(Econd_nm).confidence,...
                GLMprm.fbk.(EpEm_nm).(RP_nm).(Econd_nm).trialN] = deal(0);
            GLMprm.fbk.(EpEm_nm).(RP_nm).(Econd_nm).conf_mdl = '';
        end % effort loop
    end % RP loop
end % task loop
    
%% define variables according to GLM number
Epm = {'Ep','Em'};
RP_conds = {'R','P'};
switch GLM
    case 1
        % general parameters
        GLMprm.gal.orth_vars = 1;
        % cross
        GLMprm.model_onset.Ep.preChoiceCross = 'stick';
        GLMprm.model_onset.Em.preChoiceCross = 'stick';
        % choice
        GLMprm.model_onset.Ep.choice = 'stick';
        GLMprm.model_onset.Em.choice = 'stick';
        GLMprm.choice.Ep.RP.E.money_sum = 1;
        GLMprm.choice.Em.RP.E.money_sum = 1;
        GLMprm.choice.Ep.RP.E.E_sum = 1;
        GLMprm.choice.Em.RP.E.E_sum = 1;
        GLMprm.choice.Ep.RP.E.RT = 1;
        GLMprm.choice.Em.RP.E.RT = 1;
        % disp chosen
        GLMprm.model_onset.Ep.chosen = 'stick';
        GLMprm.model_onset.Em.chosen = 'stick';
        GLMprm.chosen.Ep.RP.E.money_chosen = 1;
        GLMprm.chosen.Em.RP.E.money_chosen = 1;
        GLMprm.chosen.Ep.RP.E.E_chosen = 1;
        GLMprm.chosen.Em.RP.E.E_chosen = 1;
        % effort perf
        GLMprm.model_onset.Ep.Eperf = 'stick';
        GLMprm.model_onset.Em.Eperf = 'stick';
        % feedback
        GLMprm.model_onset.Ep.fbk = 'stick';
        GLMprm.model_onset.Em.fbk = 'stick';
        GLMprm.fbk.Ep.RP.E.money_obtained = 1;
        GLMprm.fbk.Em.RP.E.money_obtained = 1;
    case 2
        % general parameters
        GLMprm.gal.orth_vars = 1;
        % cross
        GLMprm.model_onset.Ep.preChoiceCross = 'stick';
        GLMprm.model_onset.Em.preChoiceCross = 'stick';
        % choice
        GLMprm.model_onset.Ep.choice = 'stick';
        GLMprm.model_onset.Em.choice = 'stick';
        GLMprm.choice.Ep.RP.E.RT = 1;
        GLMprm.choice.Em.RP.E.RT = 1;
        % disp chosen
        GLMprm.model_onset.Ep.chosen = 'stick';
        GLMprm.model_onset.Em.chosen = 'stick';
        GLMprm.chosen.Ep.RP.E.money_chosen = 1;
        GLMprm.chosen.Em.RP.E.money_chosen = 1;
        GLMprm.chosen.Ep.RP.E.E_chosen = 1;
        GLMprm.chosen.Em.RP.E.E_chosen = 1;
        % effort perf
        GLMprm.model_onset.Ep.Eperf = 'stick';
        GLMprm.model_onset.Em.Eperf = 'stick';
        % feedback
        GLMprm.model_onset.Ep.fbk = 'stick';
        GLMprm.model_onset.Em.fbk = 'stick';
        GLMprm.fbk.Ep.RP.E.win_vs_loss = 1;
        GLMprm.fbk.Em.RP.E.win_vs_loss = 1;
    case 3
        % general parameters
        GLMprm.gal.orth_vars = 1;
        % cross
        GLMprm.model_onset.Ep.preChoiceCross = 'stick';
        GLMprm.model_onset.Em.preChoiceCross = 'stick';
        % choice
        GLMprm.model_onset.Ep.choice = 'stick';
        GLMprm.model_onset.Em.choice = 'stick';
        GLMprm.choice.Ep.RP.E.R_vs_P = 1;
        GLMprm.choice.Em.RP.E.R_vs_P = 1;
        GLMprm.choice.Ep.RP.E.RT = 1;
        GLMprm.choice.Em.RP.E.RT = 1;
        % disp chosen
        GLMprm.model_onset.Ep.chosen = 'stick';
        GLMprm.model_onset.Em.chosen = 'stick';
        GLMprm.chosen.Ep.RP.E.money_chosen = 1;
        GLMprm.chosen.Em.RP.E.money_chosen = 1;
        GLMprm.chosen.Ep.RP.E.E_chosen = 1;
        GLMprm.chosen.Em.RP.E.E_chosen = 1;
        % effort perf
        GLMprm.model_onset.Ep.Eperf = 'stick';
        GLMprm.model_onset.Em.Eperf = 'stick';
        % feedback
        GLMprm.model_onset.Ep.fbk = 'stick';
        GLMprm.model_onset.Em.fbk = 'stick';
    case 4
        % general parameters
        GLMprm.gal.orth_vars = 1;
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_vs_P = 1;
            GLMprm.choice.(Epm_nm).RP.E.money_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            % effort perf
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
        end
    case 5
        % general parameters
        GLMprm.gal.orth_vars = 1;
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_vs_P = 1;
            % effort perf
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
        end
    case 6 % testing Vchosen - Vunchosen during choice + adding temporal derivative
        % general parameters
        GLMprm.gal.orth_vars = 1;
        GLMprm.gal.add_drv = 1;
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.money_ch_min_unch = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_ch_min_unch = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end
    case 7 % R/P split, Vch/Vunch R/P and E levels during choice
        GLMprm.gal.orth_vars = 1;
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % fixation cross
            GLMprm.model_onset.(Epm_nm).preChoiceCross = 'stick';
            GLMprm.model_onset.(Epm_nm).preEffortCross = 'stick';
            % choice - split R/P and use R/P and E levels
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RPpool = 0;
            for iRP = 1:length(RP_conds)
                RP_nm = RP_conds{iRP};
                GLMprm.choice.(Epm_nm).(RP_nm).E.money_chosen = 3;
                GLMprm.choice.(Epm_nm).(RP_nm).E.money_unchosen = 3;
                GLMprm.choice.(Epm_nm).(RP_nm).E.E_chosen = 1;
                GLMprm.choice.(Epm_nm).(RP_nm).E.E_unchosen = 1;
                GLMprm.choice.(Epm_nm).(RP_nm).E.RT = 1;
            end
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort performance
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback - split R/P
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
            GLMprm.fbk.(Epm_nm).RP.E.win_vs_loss = 1;
        end
    case 8 % R/P split, Vch/Vunch R/P and E levels during disp chosen
        GLMprm.gal.orth_vars = 1;
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % fixation cross
            GLMprm.model_onset.(Epm_nm).preChoiceCross = 'stick';
            GLMprm.model_onset.(Epm_nm).preEffortCross = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen - split R/P and use R/P and E levels
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            GLMprm.chosen.(Epm_nm).RPpool = 0;
            for iRP = 1:length(RP_conds)
                RP_nm = RP_conds{iRP};
                GLMprm.chosen.(Epm_nm).(RP_nm).E.money_chosen = 3;
                GLMprm.chosen.(Epm_nm).(RP_nm).E.money_unchosen = 3;
                GLMprm.chosen.(Epm_nm).(RP_nm).E.E_chosen = 1;
                GLMprm.chosen.(Epm_nm).(RP_nm).E.E_unchosen = 1;
            end
            % effort performance
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback - split R/P
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
            GLMprm.fbk.(Epm_nm).RP.E.win_vs_loss = 1;
        end
    case 9 % R/P split, levels of variable option in R/P and E during choice
        GLMprm.gal.orth_vars = 1;
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % fixation cross
            GLMprm.model_onset.(Epm_nm).preChoiceCross = 'stick';
            GLMprm.model_onset.(Epm_nm).preEffortCross = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RPpool = 0;
            for iRP = 1:length(RP_conds)
                RP_nm = RP_conds{iRP};
                GLMprm.choice.(Epm_nm).(RP_nm).E.money_varOption = 3;
                GLMprm.choice.(Epm_nm).(RP_nm).E.E_varOption = 1;
                GLMprm.choice.(Epm_nm).(RP_nm).E.RT = 1;
            end
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort performance
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback - split R/P
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
            GLMprm.fbk.(Epm_nm).RP.E.win_vs_loss = 1;
        end
    case 10 % model money amounts during performance instead of choice periode
        % (VS should be triggered by higher rewards)
        GLMprm.gal.orth_vars = 1;
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % fixation cross
            GLMprm.model_onset.(Epm_nm).preChoiceCross = 'stick';
            GLMprm.model_onset.(Epm_nm).preEffortCross = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort performance
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            GLMprm.Eperf.(Epm_nm).RPpool = 0;
            for iRP = 1:length(RP_conds)
                RP_nm = RP_conds{iRP};
                GLMprm.Eperf.(Epm_nm).(RP_nm).E.money_chosen = 3;
                GLMprm.Eperf.(Epm_nm).(RP_nm).E.E_chosen = 1;
            end
            % feedback - split R/P
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
            GLMprm.fbk.(Epm_nm).RP.E.win_vs_loss = 1;
        end
    case 11 % R/P split, Vch R/P and E levels during disp chosen, like GLM8 but without the Vunchosen
        GLMprm.gal.orth_vars = 1;
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % fixation cross
            GLMprm.model_onset.(Epm_nm).preChoiceCross = 'stick';
            GLMprm.model_onset.(Epm_nm).preEffortCross = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen - split R/P and use R/P and E levels
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            GLMprm.chosen.(Epm_nm).RPpool = 0;
            for iRP = 1:length(RP_conds)
                RP_nm = RP_conds{iRP};
                GLMprm.chosen.(Epm_nm).(RP_nm).E.money_chosen = 3;
                GLMprm.chosen.(Epm_nm).(RP_nm).E.E_chosen = 1;
            end
            % effort performance
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback - split R/P
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
            GLMprm.fbk.(Epm_nm).RP.E.win_vs_loss = 1;
        end
    case 12 % Vch/R-P/VE during dispChosen option
        GLMprm.gal.orth_vars = 1;
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen - split R/P and use R/P and E levels
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            GLMprm.chosen.(Epm_nm).RPpool = 1;
            GLMprm.chosen.(Epm_nm).RP.E.money_chosen = 1;
            GLMprm.chosen.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.chosen.(Epm_nm).RP.E.R_vs_P = 1;
            % effort performance
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback - split R/P
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end
    case 13 % Rch/Ech during dispChosen option splitting R and P trials and using levels instead of money
        GLMprm.gal.orth_vars = 1;
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen - split R/P and use R/P and E levels
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            GLMprm.chosen.(Epm_nm).RPpool = 0;
            for iRP = 1:length(RP_conds)
                RP_nm = RP_conds{iRP};
                GLMprm.chosen.(Epm_nm).(RP_nm).E.money_chosen = 3;
                GLMprm.chosen.(Epm_nm).(RP_nm).E.E_chosen = 1;
            end
            % effort performance
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback - split R/P
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end
    case 14 % same as GLM 13 but with temporal derivative
        GLMprm.gal.orth_vars = 1;
        GLMprm.gal.add_drv = 1;
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen - split R/P and use R/P and E levels
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            GLMprm.chosen.(Epm_nm).RPpool = 0;
            for iRP = 1:length(RP_conds)
                RP_nm = RP_conds{iRP};
                GLMprm.chosen.(Epm_nm).(RP_nm).E.money_chosen = 3;
                GLMprm.chosen.(Epm_nm).(RP_nm).E.E_chosen = 1;
            end
            % effort performance
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end
    case 15 % like GLM 4 but feedback period modelled as well
        % general parameters
        GLMprm.gal.orth_vars = 1;
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_vs_P = 1;
            GLMprm.choice.(Epm_nm).RP.E.money_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            % effort perf
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end
    case 16 % net value + RT during choice (pooling R and P)
        % general parameters
        GLMprm.gal.orth_vars = 1;
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.NV_varOption = 1;
            GLMprm.choice.(Epm_nm).RP.E.NV_mdl = 'mdl_2';
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end
    case 17 % like GLM 16 but R and P split
        % general parameters
        GLMprm.gal.orth_vars = 1;
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
                GLMprm.choice.(Epm_nm).RPpool = 0;
            for iRP = 1:length(RP_conds)
                RP_nm = RP_conds{iRP};
                GLMprm.choice.(Epm_nm).(RP_nm).E.NV_varOption = 1;
                GLMprm.choice.(Epm_nm).(RP_nm).E.NV_mdl = 'mdl_2';
                GLMprm.choice.(Epm_nm).(RP_nm).E.RT = 1;
            end
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end
    case 18 % like GLM 16 but during chosen period
        % general parameters
        GLMprm.gal.orth_vars = 1;
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            GLMprm.chosen.(Epm_nm).RPpool = 0;
            for iRP = 1:length(RP_conds)
                RP_nm = RP_conds{iRP};
                GLMprm.chosen.(Epm_nm).(RP_nm).E.NV_chosen = 1;
                GLMprm.chosen.(Epm_nm).(RP_nm).E.NV_mdl = 'mdl_2';
                GLMprm.chosen.(Epm_nm).(RP_nm).E.RT = 1;
            end
            % effort perf
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end
    case 19 % Rvar/Evar/trialN*E var/RT during choice
        % general parameters
        GLMprm.gal.orth_vars = 1;
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RPpool = 0;
            for iRP = 1:length(RP_conds)
                RP_nm = RP_conds{iRP};
                GLMprm.choice.(Epm_nm).(RP_nm).E.money_varOption = 3;
                GLMprm.choice.(Epm_nm).(RP_nm).E.E_varOption = 1;
                GLMprm.choice.(Epm_nm).(RP_nm).E.trialN = 3;
                GLMprm.choice.(Epm_nm).(RP_nm).E.RT = 1;
            end
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end
    case 20 % Rvar/Evar/trialN*E var/RT during chosen
        % general parameters
        GLMprm.gal.orth_vars = 1;
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            GLMprm.chosen.(Epm_nm).RPpool = 0;
            for iRP = 1:length(RP_conds)
                RP_nm = RP_conds{iRP};
                GLMprm.chosen.(Epm_nm).(RP_nm).E.money_varOption = 3;
                GLMprm.chosen.(Epm_nm).(RP_nm).E.E_varOption = 1;
                GLMprm.chosen.(Epm_nm).(RP_nm).E.trialN = 3;
                GLMprm.chosen.(Epm_nm).(RP_nm).E.RT = 1;
            end
            % effort perf
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end
    case 21 % NV during performance R/P split
        % general parameters
        GLMprm.gal.orth_vars = 0; % no orthogonalization (nothing to orthogonalize anyway)
        GLMprm.gal.zPerRun = 1; % zscore net value per run
        GLMprm.gal.add_drv = 1; % temporal derivative
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            GLMprm.Eperf.(Epm_nm).RPpool = 0;
            for iRP = 1:length(RP_conds)
                RP_nm = RP_conds{iRP};
                GLMprm.Eperf.(Epm_nm).(RP_nm).E.NV_chosen = 1;
                GLMprm.Eperf.(Epm_nm).(RP_nm).E.NV_mdl = 'mdl_4';
            end % R/P
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 22 % chosen option R/E/trial number/RT during choice
        % general parameters
        GLMprm.gal.orth_vars = 0; % no orthogonalization (nothing to orthogonalize anyway)
        GLMprm.gal.zPerRun = 1; % zscore net value per run
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RPpool = 0;
            for iRP = 1:length(RP_conds)
                RP_nm = RP_conds{iRP};
                GLMprm.choice.(Epm_nm).(RP_nm).E.money_chosen = 1;
                GLMprm.choice.(Epm_nm).(RP_nm).E.E_chosen = 1;
                GLMprm.choice.(Epm_nm).(RP_nm).E.trialN = 1;
                GLMprm.choice.(Epm_nm).(RP_nm).E.RT = 1;
            end % R/P
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 23 % chosen option R/E/trial number/RT during chosen
        % general parameters
        GLMprm.gal.orth_vars = 0; % no orthogonalization (nothing to orthogonalize anyway)
        GLMprm.gal.zPerRun = 1; % zscore net value per run
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            GLMprm.chosen.(Epm_nm).RPpool = 0;
            for iRP = 1:length(RP_conds)
                RP_nm = RP_conds{iRP};
                GLMprm.chosen.(Epm_nm).(RP_nm).E.money_chosen = 1;
                GLMprm.chosen.(Epm_nm).(RP_nm).E.E_chosen = 1;
                GLMprm.chosen.(Epm_nm).(RP_nm).E.trialN = 1;
                GLMprm.chosen.(Epm_nm).(RP_nm).E.RT = 1;
            end % R/P
            % effort perf
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 24 % chosen option R/E/trial number during effort performance
        % general parameters
        GLMprm.gal.orth_vars = 0; % no orthogonalization (nothing to orthogonalize anyway)
        GLMprm.gal.zPerRun = 1; % zscore net value per run
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            GLMprm.Eperf.(Epm_nm).RPpool = 0;
            for iRP = 1:length(RP_conds)
                RP_nm = RP_conds{iRP};
                GLMprm.Eperf.(Epm_nm).(RP_nm).E.money_chosen = 1;
                GLMprm.Eperf.(Epm_nm).(RP_nm).E.E_chosen = 1;
                GLMprm.Eperf.(Epm_nm).(RP_nm).E.trialN = 1;
            end % R/P
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 25 % NV chosen option/RT during choice
        % general parameters
        GLMprm.gal.orth_vars = 0; % no orthogonalization (nothing to orthogonalize anyway)
        GLMprm.gal.zPerRun = 1; % zscore net value per run
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RPpool = 0;
            for iRP = 1:length(RP_conds)
                RP_nm = RP_conds{iRP};
                GLMprm.choice.(Epm_nm).(RP_nm).E.NV_chosen = 1;
                GLMprm.choice.(Epm_nm).(RP_nm).E.NV_mdl = 'mdl_4';
                GLMprm.choice.(Epm_nm).(RP_nm).E.RT = 1;
            end
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 26 % NV chosen option/RT during chosen
        % general parameters
        GLMprm.gal.orth_vars = 0; % no orthogonalization (nothing to orthogonalize anyway)
        GLMprm.gal.zPerRun = 1; % zscore net value per run
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            GLMprm.chosen.(Epm_nm).RPpool = 0;
            for iRP = 1:length(RP_conds)
                RP_nm = RP_conds{iRP};
                GLMprm.chosen.(Epm_nm).(RP_nm).E.NV_chosen = 1;
                GLMprm.chosen.(Epm_nm).(RP_nm).E.NV_mdl = 'mdl_4';
                GLMprm.chosen.(Epm_nm).(RP_nm).E.RT = 1;
            end
            % effort perf
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 27 % NV chosen option during effort performance
        % general parameters
        GLMprm.gal.orth_vars = 0; % no orthogonalization (nothing to orthogonalize anyway)
        GLMprm.gal.zPerRun = 1; % zscore net value per run
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            GLMprm.Eperf.(Epm_nm).RPpool = 0;
            for iRP = 1:length(RP_conds)
                RP_nm = RP_conds{iRP};
                GLMprm.Eperf.(Epm_nm).(RP_nm).E.NV_chosen = 1;
                GLMprm.Eperf.(Epm_nm).(RP_nm).E.NV_mdl = 'mdl_4';
            end
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 28 % non-default option R/E/trial number/RT during choice
        % general parameters
        GLMprm.gal.orth_vars = 0; % no orthogonalization (nothing to orthogonalize anyway)
        GLMprm.gal.zPerRun = 1; % zscore net value per run
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RPpool = 0;
            for iRP = 1:length(RP_conds)
                RP_nm = RP_conds{iRP};
                GLMprm.choice.(Epm_nm).(RP_nm).E.money_varOption = 1;
                GLMprm.choice.(Epm_nm).(RP_nm).E.E_varOption = 1;
                GLMprm.choice.(Epm_nm).(RP_nm).E.trialN = 1;
                GLMprm.choice.(Epm_nm).(RP_nm).E.RT = 1;
            end
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 29 % non-default option R/E/trial number/RT during chosen
        % general parameters
        GLMprm.gal.orth_vars = 0; % no orthogonalization (nothing to orthogonalize anyway)
        GLMprm.gal.zPerRun = 1; % zscore net value per run
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            GLMprm.chosen.(Epm_nm).RPpool = 0;
            for iRP = 1:length(RP_conds)
                RP_nm = RP_conds{iRP};
                GLMprm.chosen.(Epm_nm).(RP_nm).E.money_varOption = 1;
                GLMprm.chosen.(Epm_nm).(RP_nm).E.E_varOption = 1;
                GLMprm.chosen.(Epm_nm).(RP_nm).E.trialN = 1;
                GLMprm.chosen.(Epm_nm).(RP_nm).E.RT = 1;
            end
            % effort perf
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 30 % non-default option R/E/trial number during effort performance
        % general parameters
        GLMprm.gal.orth_vars = 0; % no orthogonalization (nothing to orthogonalize anyway)
        GLMprm.gal.zPerRun = 1; % zscore net value per run
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            GLMprm.Eperf.(Epm_nm).RPpool = 0;
            for iRP = 1:length(RP_conds)
                RP_nm = RP_conds{iRP};
                GLMprm.Eperf.(Epm_nm).(RP_nm).E.money_varOption = 1;
                GLMprm.Eperf.(Epm_nm).(RP_nm).E.E_varOption = 1;
                GLMprm.Eperf.(Epm_nm).(RP_nm).E.trialN = 1;
            end
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 31 % NV non-default option/RT during choice
        % general parameters
        GLMprm.gal.orth_vars = 0; % no orthogonalization (nothing to orthogonalize anyway)
        GLMprm.gal.zPerRun = 1; % zscore net value per run
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RPpool = 0;
            for iRP = 1:length(RP_conds)
                RP_nm = RP_conds{iRP};
                GLMprm.choice.(Epm_nm).(RP_nm).E.NV_varOption = 1;
                GLMprm.choice.(Epm_nm).(RP_nm).E.NV_mdl = 'mdl_4';
                GLMprm.choice.(Epm_nm).(RP_nm).E.RT = 1;
            end
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 32 % NV non-default option/RT during chosen
        % general parameters
        GLMprm.gal.orth_vars = 0; % no orthogonalization (nothing to orthogonalize anyway)
        GLMprm.gal.zPerRun = 1; % zscore net value per run
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            GLMprm.chosen.(Epm_nm).RPpool = 0;
            for iRP = 1:length(RP_conds)
                RP_nm = RP_conds{iRP};
                GLMprm.chosen.(Epm_nm).(RP_nm).E.NV_varOption = 1;
                GLMprm.chosen.(Epm_nm).(RP_nm).E.NV_mdl = 'mdl_4';
                GLMprm.chosen.(Epm_nm).(RP_nm).E.RT = 1;
            end
            % effort perf
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 33 % NV non-default option during effort performance
        % general parameters
        GLMprm.gal.orth_vars = 0; % no orthogonalization (nothing to orthogonalize anyway)
        GLMprm.gal.zPerRun = 1; % zscore net value per run
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            GLMprm.Eperf.(Epm_nm).RPpool = 0;
            for iRP = 1:length(RP_conds)
                RP_nm = RP_conds{iRP};
                GLMprm.Eperf.(Epm_nm).(RP_nm).E.NV_varOption = 1;
                GLMprm.Eperf.(Epm_nm).(RP_nm).E.NV_mdl = 'mdl_4';
            end
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 34 % RT during choice/chosen option R/E/trial number during chosen and pool RP trials
        % general parameters
        GLMprm.gal.orth_vars = 0; % no orthogonalization (nothing to orthogonalize anyway)
        GLMprm.gal.zPerRun = 1; % zscore net value per run
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RPpool = 1;
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            GLMprm.chosen.(Epm_nm).RPpool = 1;
            GLMprm.chosen.(Epm_nm).RP.E.money_chosen = 1;
            GLMprm.chosen.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.chosen.(Epm_nm).RP.E.trialN = 1;
            % effort perf
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 35 % same as GLM 34 R/E/trial number during choice + RT during choice
        % general parameters
        GLMprm.gal.orth_vars = 0; % no orthogonalization (nothing to orthogonalize anyway)
        GLMprm.gal.zPerRun = 1; % zscore net value per run
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RPpool = 1;
            GLMprm.choice.(Epm_nm).RP.E.money_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.trialN = 1;
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 36 % RT during choice/NV chosen during chosen
        % general parameters
        GLMprm.gal.orth_vars = 0; % no orthogonalization (nothing to orthogonalize anyway)
        GLMprm.gal.zPerRun = 1; % zscore net value per run
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RPpool = 1;
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            GLMprm.chosen.(Epm_nm).RPpool = 1;
            GLMprm.chosen.(Epm_nm).RP.E.NV_chosen = 1;
            GLMprm.chosen.(Epm_nm).RP.E.NV_mdl = 'mdl_4';
            % effort perf
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 37 % choice: NV variable option/confidence inferred/RT; chosen: NV chosen
        % general parameters
        GLMprm.gal.orth_vars = 0; % no orthogonalization (nothing to orthogonalize anyway)
        GLMprm.gal.zPerRun = 1; % zscore net value per run
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RPpool = 1;
            GLMprm.choice.(Epm_nm).RP.E.NV_varOption = 1;
            GLMprm.choice.(Epm_nm).RP.E.NV_mdl = 'mdl_4';
            GLMprm.choice.(Epm_nm).RP.E.confidence = 2;
            GLMprm.choice.(Epm_nm).RP.E.conf_mdl = 'mdl_4';
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            GLMprm.chosen.(Epm_nm).RPpool = 1;
            GLMprm.chosen.(Epm_nm).RP.E.NV_chosen = 1;
            GLMprm.chosen.(Epm_nm).RP.E.NV_mdl = 'mdl_4';
            % effort perf
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
            
        end % physical/mental loop
    case 38 % dACC special= choice: confidence/RT; chosen+effort: incentive + effort
        % general parameters
        GLMprm.gal.orth_vars = 1; % orthogonalization to remove uncertainty from RT regressor
        GLMprm.gal.zPerRun = 1; % zscore net value per run
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RPpool = 1;
            GLMprm.choice.(Epm_nm).RP.E.confidence = 2;
            GLMprm.choice.(Epm_nm).RP.E.conf_mdl = 'mdl_4';
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'boxcar_bis';
            GLMprm.chosen.(Epm_nm).RPpool = 1;
            GLMprm.chosen.(Epm_nm).RP.E.money_chosen = 3;
            GLMprm.chosen.(Epm_nm).RP.E.E_chosen = 1;
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 39 % dACC special bis = choice: confidence/RT; effort cross+effort: incentive + effort
        % general parameters
        GLMprm.gal.orth_vars = 1; % orthogonalization to remove uncertainty from RT regressor
        GLMprm.gal.zPerRun = 1; % zscore net value per run
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % initial cross
            GLMprm.model_onset.(Epm_nm).preChoiceCross = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RPpool = 1;
            GLMprm.choice.(Epm_nm).RP.E.confidence = 2;
            GLMprm.choice.(Epm_nm).RP.E.conf_mdl = 'mdl_4';
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % pre-effort cross (effort preparation)
            GLMprm.model_onset.(Epm_nm).preEffortCross = 'boxcar_bis';
            GLMprm.preEffortCross.(Epm_nm).RPpool = 1;
            GLMprm.preEffortCross.(Epm_nm).RP.E.money_chosen = 3;
            GLMprm.preEffortCross.(Epm_nm).RP.E.E_chosen = 1;
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 40 % dACC special bis = choice: confidence/RT; effort cross and effort: incentive + effort
        % general parameters
        GLMprm.gal.orth_vars = 1; % orthogonalization to remove uncertainty from RT regressor
        GLMprm.gal.zPerRun = 1; % zscore net value per run
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % initial cross
            GLMprm.model_onset.(Epm_nm).preChoiceCross = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RPpool = 1;
            GLMprm.choice.(Epm_nm).RP.E.confidence = 2;
            GLMprm.choice.(Epm_nm).RP.E.conf_mdl = 'mdl_4';
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % pre-effort cross (effort preparation)
            GLMprm.model_onset.(Epm_nm).preEffortCross = 'stick';
            GLMprm.preEffortCross.(Epm_nm).RPpool = 1;
            GLMprm.preEffortCross.(Epm_nm).RP.E.money_chosen = 3;
            GLMprm.preEffortCross.(Epm_nm).RP.E.E_chosen = 1;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            GLMprm.Eperf.(Epm_nm).RPpool = 1;
            GLMprm.Eperf.(Epm_nm).RP.E.money_chosen = 3;
            GLMprm.Eperf.(Epm_nm).RP.E.E_chosen = 1;
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 41 % dACC special bis = choice: confidence/RT; effort cross and effort: incentive + effort
        % general parameters
        GLMprm.gal.orth_vars = 1; % orthogonalization to remove uncertainty from RT regressor
        GLMprm.gal.zPerRun = 1; % zscore net value per run
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % initial cross
            GLMprm.model_onset.(Epm_nm).preChoiceCross = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RPpool = 1;
            GLMprm.choice.(Epm_nm).RP.E.confidence = 2;
            GLMprm.choice.(Epm_nm).RP.E.conf_mdl = 'mdl_4';
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            GLMprm.chosen.(Epm_nm).RPpool = 1;
            GLMprm.chosen.(Epm_nm).RP.E.money_chosen = 3;
            GLMprm.chosen.(Epm_nm).RP.E.E_chosen = 1;
            % pre-effort cross (effort preparation)
            GLMprm.model_onset.(Epm_nm).preEffortCross = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            GLMprm.Eperf.(Epm_nm).RPpool = 1;
            GLMprm.Eperf.(Epm_nm).RP.E.money_chosen = 3;
            GLMprm.Eperf.(Epm_nm).RP.E.E_chosen = 1;
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 42
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % initial cross
            GLMprm.model_onset.(Epm_nm).preChoiceCross = 'boxcar';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'boxcar';
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'boxcar';
            GLMprm.chosen.(Epm_nm).RP.E.confidence = 2;
            GLMprm.chosen.(Epm_nm).RP.E.conf_mdl = 'mdl_4';
            % pre-effort cross
            GLMprm.model_onset.(Epm_nm).preEffortCross = 'boxcar';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'boxcar';
            GLMprm.Eperf.(Epm_nm).RPpool = 0;
            for iRP = 1:length(RP_conds)
                RP_nm = RP_conds{iRP};
                GLMprm.Eperf.(Epm_nm).(RP_nm).E.money_chosen = 1;
                GLMprm.Eperf.(Epm_nm).(RP_nm).E.E_chosen = 1;
            end
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'boxcar';
            GLMprm.fbk.(Epm_nm).RP.E.win_vs_loss = 1;
        end % physical/mental  loop
    case 43 % modulation of effort execution by effort level to see if dmPFC
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % initial cross
            GLMprm.model_onset.(Epm_nm).preChoiceCross = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            GLMprm.chosen.(Epm_nm).RP.E.confidence = 2;
            GLMprm.chosen.(Epm_nm).RP.E.conf_mdl = 'mdl_4';
            % pre-effort cross
            GLMprm.model_onset.(Epm_nm).preEffortCross = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            GLMprm.Eperf.(Epm_nm).RP.E.money_chosen = 1;
            GLMprm.Eperf.(Epm_nm).RP.E.E_chosen = 1;
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental  loop
    case 44 % same as GLM43 but to look for effort preparation during chosen option instead of effort execution
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % initial cross
            GLMprm.model_onset.(Epm_nm).preChoiceCross = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.confidence = 2;
            GLMprm.choice.(Epm_nm).RP.E.conf_mdl = 'mdl_4';
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            GLMprm.chosen.(Epm_nm).RP.E.money_chosen = 1;
            GLMprm.chosen.(Epm_nm).RP.E.E_chosen = 1;
            % pre-effort cross
            GLMprm.model_onset.(Epm_nm).preEffortCross = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental  loop
    case 45 % check chosen - default frame
        % general parameters
        GLMprm.gal.orth_vars = 1;
        GLMprm.gal.zPerRun = 1;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % initial cross
            GLMprm.model_onset.(Epm_nm).preChoiceCross = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.money_ch_min_fixOption = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_ch_min_fixOption = 1;
            GLMprm.choice.(Epm_nm).RP.E.confidence = 2;
            GLMprm.choice.(Epm_nm).RP.E.conf_mdl = 'mdl_4';
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % pre-effort cross
            GLMprm.model_onset.(Epm_nm).preEffortCross = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental  loop
    case 46 % check non-default option during choice
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % initial cross
            GLMprm.model_onset.(Epm_nm).preChoiceCross = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RPpool = 0;
            for iRP = 1:length(RP_conds)
                RP_nm = RP_conds{iRP};
                GLMprm.choice.(Epm_nm).(RP_nm).E.money_varOption = 1;
                GLMprm.choice.(Epm_nm).(RP_nm).E.E_varOption = 1;
                GLMprm.choice.(Epm_nm).(RP_nm).E.confidence = 2;
                GLMprm.choice.(Epm_nm).(RP_nm).E.conf_mdl = 'mdl_4';
                GLMprm.choice.(Epm_nm).(RP_nm).E.RT = 1;
            end
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % pre-effort cross
            GLMprm.model_onset.(Epm_nm).preEffortCross = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental  loop
    case 47 % same as GLM 46 but pooling R and P trials now to increase power
        % also removed the zscoring of the variables
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % initial cross
            GLMprm.model_onset.(Epm_nm).preChoiceCross = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.money_varOption = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_varOption = 1;
            GLMprm.choice.(Epm_nm).RP.E.confidence = 2;
            GLMprm.choice.(Epm_nm).RP.E.conf_mdl = 'mdl_4';
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % pre-effort cross
            GLMprm.model_onset.(Epm_nm).preEffortCross = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental  loop
    case 48 % same as GLM 47 but orthogonalizing variables
        % general parameters
        GLMprm.gal.orth_vars = 1;
        GLMprm.gal.zPerRun = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % initial cross
            GLMprm.model_onset.(Epm_nm).preChoiceCross = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.money_varOption = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_varOption = 1;
            GLMprm.choice.(Epm_nm).RP.E.confidence = 2;
            GLMprm.choice.(Epm_nm).RP.E.conf_mdl = 'mdl_4';
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % pre-effort cross
            GLMprm.model_onset.(Epm_nm).preEffortCross = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental  loop
    case 49 % same as GLM 48 but zscoring variables
        % general parameters
        GLMprm.gal.orth_vars = 1;
        GLMprm.gal.zPerRun = 1;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % initial cross
            GLMprm.model_onset.(Epm_nm).preChoiceCross = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.money_varOption = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_varOption = 1;
            GLMprm.choice.(Epm_nm).RP.E.confidence = 2;
            GLMprm.choice.(Epm_nm).RP.E.conf_mdl = 'mdl_4';
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % pre-effort cross
            GLMprm.model_onset.(Epm_nm).preEffortCross = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental  loop
    case 50 % NVnd/Conf/DT during choice/ Ech anticipation during chosen/Money obtained during feedback
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % initial cross
            GLMprm.model_onset.(Epm_nm).preChoiceCross = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.NV_varOption = 1;
            GLMprm.choice.(Epm_nm).RP.E.NV_mdl = 'mdl_4';
            GLMprm.choice.(Epm_nm).RP.E.confidence = 2;
            GLMprm.choice.(Epm_nm).RP.E.conf_mdl = 'mdl_4';
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            GLMprm.chosen.(Epm_nm).RP.E.E_chosen = 1;
            % pre-effort cross
            GLMprm.model_onset.(Epm_nm).preEffortCross = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
            GLMprm.fbk.(Epm_nm).RP.E.money_obtained = 1;
        end % physical/mental  loop
    case 51 % same as GLM 50 but decomposing net value in Rnd/End/trialN*End
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % initial cross
            GLMprm.model_onset.(Epm_nm).preChoiceCross = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.money_varOption = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_varOption = 1;
            GLMprm.choice.(Epm_nm).RP.E.trialN = 4;
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            GLMprm.chosen.(Epm_nm).RP.E.E_chosen = 1;
            % pre-effort cross
            GLMprm.model_onset.(Epm_nm).preEffortCross = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
            GLMprm.fbk.(Epm_nm).RP.E.money_obtained = 1;
        end % physical/mental  loop
    case 52 % same as GLM 50 but with behavioral model 3 instead of model 4
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % initial cross
            GLMprm.model_onset.(Epm_nm).preChoiceCross = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.NV_varOption = 1;
            GLMprm.choice.(Epm_nm).RP.E.NV_mdl = 'mdl_3';
            GLMprm.choice.(Epm_nm).RP.E.confidence = 2;
            GLMprm.choice.(Epm_nm).RP.E.conf_mdl = 'mdl_3';
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            GLMprm.chosen.(Epm_nm).RP.E.E_chosen = 1;
            % pre-effort cross
            GLMprm.model_onset.(Epm_nm).preEffortCross = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
            GLMprm.fbk.(Epm_nm).RP.E.money_obtained = 1;
        end % physical/mental  loop
    case 53 % chosen option during choice
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % initial cross
            GLMprm.model_onset.(Epm_nm).preChoiceCross = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.money_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.confidence = 2;
            GLMprm.choice.(Epm_nm).RP.E.conf_mdl = 'mdl_3';
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % pre-effort cross
            GLMprm.model_onset.(Epm_nm).preEffortCross = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental  loop
    case 54 % like GLM45 but with DT as first regressor
        % general parameters
        GLMprm.gal.orth_vars = 1;
        GLMprm.gal.zPerRun = 1;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % initial cross
            GLMprm.model_onset.(Epm_nm).preChoiceCross = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.RT = 4;
            GLMprm.choice.(Epm_nm).RP.E.money_ch_min_fixOption = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_ch_min_fixOption = 1;
            GLMprm.choice.(Epm_nm).RP.E.confidence = 2;
            GLMprm.choice.(Epm_nm).RP.E.conf_mdl = 'mdl_4';
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % pre-effort cross
            GLMprm.model_onset.(Epm_nm).preEffortCross = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental  loop
    case 55 % like GLM45 but replacing confidence model4 by model3
        % general parameters
        GLMprm.gal.orth_vars = 1;
        GLMprm.gal.zPerRun = 1;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % initial cross
            GLMprm.model_onset.(Epm_nm).preChoiceCross = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.money_ch_min_fixOption = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_ch_min_fixOption = 1;
            GLMprm.choice.(Epm_nm).RP.E.confidence = 2;
            GLMprm.choice.(Epm_nm).RP.E.conf_mdl = 'mdl_3';
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % pre-effort cross
            GLMprm.model_onset.(Epm_nm).preEffortCross = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental  loop
    case 56 % like GLM54 but replacing confidence model4 by model3
        % (in other words, like GLM55 but with DT as first regressor)
        % general parameters
        GLMprm.gal.orth_vars = 1;
        GLMprm.gal.zPerRun = 1;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % initial cross
            GLMprm.model_onset.(Epm_nm).preChoiceCross = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.RT = 4;
            GLMprm.choice.(Epm_nm).RP.E.money_ch_min_fixOption = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_ch_min_fixOption = 1;
            GLMprm.choice.(Epm_nm).RP.E.confidence = 2;
            GLMprm.choice.(Epm_nm).RP.E.conf_mdl = 'mdl_3';
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % pre-effort cross
            GLMprm.model_onset.(Epm_nm).preEffortCross = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental  loop
    case 57 % like GLM55 but no orthogonalization
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % initial cross
            GLMprm.model_onset.(Epm_nm).preChoiceCross = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.money_ch_min_fixOption = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_ch_min_fixOption = 1;
            GLMprm.choice.(Epm_nm).RP.E.confidence = 2;
            GLMprm.choice.(Epm_nm).RP.E.conf_mdl = 'mdl_3';
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % pre-effort cross
            GLMprm.model_onset.(Epm_nm).preEffortCross = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental  loop
    case 58 % model both choice and E perf
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % initial cross
            GLMprm.model_onset.(Epm_nm).preChoiceCross = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.choiceHighE = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % pre-effort cross
            GLMprm.model_onset.(Epm_nm).preEffortCross = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            GLMprm.Eperf.(Epm_nm).RP.E.E_chosen = 1;
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental  loop
    case 59 % model low vs high effort chosen separately during choice
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % initial cross
            GLMprm.model_onset.(Epm_nm).preChoiceCross = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).splitPerE = 3;
            GLMprm.choice.(Epm_nm).RP.lEch.money_varOption = 1;
            GLMprm.choice.(Epm_nm).RP.lEch.E_varOption = 1;
            GLMprm.choice.(Epm_nm).RP.hEch.money_varOption = 1;
            GLMprm.choice.(Epm_nm).RP.hEch.E_varOption = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'boxcar';
            % pre-effort cross
            GLMprm.model_onset.(Epm_nm).preEffortCross = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            GLMprm.Eperf.(Epm_nm).RP.E.E_chosen = 1;
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'boxcar';
        end % physical/mental  loop
    case 60 % model high effort choice during fixation cross and during choice
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % initial cross
            GLMprm.model_onset.(Epm_nm).preChoiceCross = 'stick';
            GLMprm.preChoiceCross.(Epm_nm).choiceHighE = 1;
            GLMprm.preChoiceCross.(Epm_nm).RT = 1;
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.choiceHighE = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'boxcar';
            % pre-effort cross
            GLMprm.model_onset.(Epm_nm).preEffortCross = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            GLMprm.Eperf.(Epm_nm).RP.E.E_chosen = 1;
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'boxcar';
        end % physical/mental  loop
    case 61 % same as GLM 60 but modeling choice high effort only during pre-choice cross and adding RT during choice
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % initial cross
            GLMprm.model_onset.(Epm_nm).preChoiceCross = 'stick';
            GLMprm.preChoiceCross.(Epm_nm).choiceHighE = 1;
            GLMprm.preChoiceCross.(Epm_nm).RT = 1;
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'boxcar';
            % pre-effort cross
            GLMprm.model_onset.(Epm_nm).preEffortCross = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            GLMprm.Eperf.(Epm_nm).RP.E.E_chosen = 1;
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'boxcar';
        end % physical/mental  loop
    case 62 % model effort chosen during choice and during performance
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % initial cross
            GLMprm.model_onset.(Epm_nm).preChoiceCross = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % pre-effort cross
            GLMprm.model_onset.(Epm_nm).preEffortCross = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            GLMprm.Eperf.(Epm_nm).RP.E.E_chosen = 1;
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental  loop
    case 63 % same as GLM 62 with slight differences (boxcar/stick)
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % initial cross
            GLMprm.model_onset.(Epm_nm).preChoiceCross = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.E_varOption = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'boxcar';
            % pre-effort cross
            GLMprm.model_onset.(Epm_nm).preEffortCross = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            GLMprm.Eperf.(Epm_nm).RP.E.E_chosen = 1;
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'boxcar';
        end % physical/mental  loop
    case 64 % onsets-only GLM
        GLMprm.gal.onsets_only = 1;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % initial cross
            GLMprm.model_onset.(Epm_nm).preChoiceCross = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'boxcar';
            % pre-effort cross
            GLMprm.model_onset.(Epm_nm).preEffortCross = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'boxcar';
        end % physical/mental  loop
    case 65 % onsets-only GLM but with less events modelled to avoid overlaps
        GLMprm.gal.onsets_only = 1;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'boxcar';
        end % physical/mental  loop
    case 66 % GLM with Ech during choice and during Eperf
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.money_ch_min_fixOption = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.confidence = 2;
            GLMprm.choice.(Epm_nm).RP.E.conf_mdl = 'mdl_3';
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            GLMprm.Eperf.(Epm_nm).RP.E.E_chosen = 1;
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'boxcar';
        end % physical/mental  loop
end % GLM number
%% warnings: check compatibility of the GLM parameters entered
isGLMokCheck(GLMprm);

end % function