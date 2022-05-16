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
%   .choice/chosen/preEffortCross/Eperf/fbk: for each event, for each task (Ep/Em) and 
%   for each condition (R/P/RP) indicate if a given regressor should be 
%   included or not.
%   .Ep/Em: physical (Ep) or mental (Em) effort task
%   .R/P/RP: reward only (R), punishment only (P) or reward and punishment
%   trials mixed (RP)
%
%       .(choice/chosen).Ep/Em.RPpool:
%       (0) split rewards and punishments as separate events
%       (1) pool reward and punishment trials
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).R_vs_P:
%       (1) 0 when punishment trial, 1 when reward trial
%
%       .choice.(Ep/Em).(R/P/RP).money_left
%       (1) money amount associated to left option
%
%       .choice.(Ep/Em).(R/P/RP).money_right
%       (1) money amount associated to right option
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).money_chosen:
%       (1) money chosen amount
%       (2) |money chosen amount|
%       (3) money levels (-4/-3/-2/-1/1/2/3/4) chosen from highest loss
%       until highest gain
%       (4) |money levels (1/2/3/4) chosen|
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).money_unchosen:
%       (1) money unchosen amount
%       (2) |money unchosen amount|
%       (3) money levels (-4/-3/-2/-1/1/2/3/4) unchosen
%       (4) |money levels (1/2/3/4) unchosen|
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).money_varOption:
%       (1) money amount non-default option
%       (2) |money amount non-default option|
%       (3) money levels  (-3/-2/-1/1/2/3) non-default option
%       (4) |money levels (1/2/3) non-default option|
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).money_ch_min_unch:
%       (1) money chosen - money unchosen amount
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).money_sum:
%       (1) money default + money non-default amount
%
%       .choice.(Ep/Em).(R/P/RP).E_left
%       (1) effort level (0/1/2/3) associated to left option
%       (2) effort difficulty associated to left option (Ep: duration to hold; Em: nb answers to give)
%
%       .choice.(Ep/Em).(R/P/RP).E_right
%       (1) effort level (0/1/2/3) associated to right option
%       (2) effort difficulty associated to right option (Ep: duration to hold; Em: nb answers to give)
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).E_chosen
%       (1) effort level (0/1/2/3) associated to chosen option
%       (2) effort difficulty associated to chosen option (Ep: duration to hold; Em: nb answers to give)
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).E_unchosen
%       (1) effort level (0/1/2/3) associated to unchosen option
%       (2) effort difficulty associated to unchosen option (Ep: duration to hold; Em: nb answers to give)
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).E_varOption
%       (1) effort level (0/1/2/3) associated to the non-default option
%       (2) effort difficulty associated to the non-default option (Ep: duration to hold; Em: nb answers to give)
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).E_ch_min_unch
%       (1) effort level (0/1/2/3) associated to the chosen option minus
%       the effort level (0/1/2/3) associated to the unchosen option
%       (2) effort difficulty associated to the chosen option (Ep: duration to hold; Em: nb answers to give)
%       minus the effort difficulty associated to the unchosen option (Ep: duration to hold; Em: nb answers to give)
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).E_sum
%       (1) sum of the effort levels (0/1/2/3) associated to both options
%       (2) sum of the effort difficulties (Ep: duration to hold; Em: nb answers to give) associated to both options
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).NV_chosen
%       (1) net value of the chosen option based on the model defined in
%       .(choice/chosen).(Ep.Em).(R/P/RP).NV_mdl (='mdl_X' or 'bayesianModel_X')
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).NV_varOption
%       (1) net value of the chosen option based on the model defined in
%       .(choice/chosen).(Ep.Em).(R/P/RP).NV_mdl (='mdl_X' or 'bayesianModel_X')
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).trialN
%       (1) trial number
%       (2) trial number*(E chosen - E non-chosen option)
%       (3) trial number*(E non-default - E default option)
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).confidence
%       (1) confidence level (0/1) given by the subject for each choice
%       (2) confidence inferred by the model (p(choice)-0.5)² for the model
%       defined in .(choice/chosen).(Ep.Em).(R/P/RP).conf_mdl (='mdl_X' or 'bayesianModel_X')
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).RT: reaction time for choice
%       (1) raw reaction time
%       (2) reaction time zscored per run
%       (3) reaction time zscored per subject across all runs
%
%
%       .Eperf.(Ep/Em).(R/P/RP).RPpool
%       (0) split rewards and punishments as separate events
%       (1) pool reward and punishment trials
%
%       .Eperf.(Ep/Em).(R/P/RP).money_chosen
%       (1) money chosen amount
%       (2) |money chosen amount|
%       (3) money levels (-4/-3/-2/-1/1/2/3/4) chosen
%       (4) |money levels (1/2/3/4) chosen|
%
%       .Eperf.(Ep/Em).(R/P/RP).E_chosen
%       (1) effort level (0/1/2/3) associated to chosen option
%       (2) effort difficulty associated to chosen option (Ep: duration to hold; Em: nb answers to give)
%
%       .Eperf.Ep.(R/P/RP).F_peak: physical effort only: force peak
%       (1) force peak in newtons
%
%       .Eperf.Ep.(R/P/RP).F_integral: physical effort only: force integral
%       (1) integral of effort performed during the effort period (sum of 
%       efforts in newtons)
%
%       .Eperf.Em.(R/P/RP).RT_avg: mental effort only: average reaction
%       time for answering N-back task
%       (1) average reaction time for answering to questions (in seconds)
%
%       .Eperf.Em.(R/P/RP).n_errors: mental effort only: number of errors
%       made per trial
%       (1) number of errors made
%
%       .Eperf.(Ep/Em).(R/P/RP).NV_chosen
%       (1) net value of the chosen option based on the model defined in
%       .Eperf.(Ep.Em).(R/P/RP).NV_mdl (='mdl_X' or 'bayesianModel_X')
%
%       .Eperf.(Ep/Em).(R/P/RP).NV_varOption
%       (1) net value of the chosen option based on the model defined in
%       .Eperf.(Ep.Em).(R/P/RP).NV_mdl (='mdl_X' or 'bayesianModel_X')
%
%       .Eperf.(Ep/Em).(R/P/RP).RT_1stAnswer
%       (1) raw reaction time for first answer (force above threshold for Ep 
%       and first answer to first digit for Em)
%
%       .Eperf.(Ep/Em).(R/P/RP).trialN
%       (1) trial number
%       (2) trial number*(E chosen - E non-chosen option)
%       (3) trial number*(E non-default - E default option)
%
%       .Eperf.(Ep/Em).(R/P/RP).confidence
%       (1) confidence level (0/1) given by the subject for each choice
%       (2) confidence inferred by the model (p(choice)-0.5)² for the model
%       defined in .Eperf.(Ep.Em).(R/P/RP).conf_mdl (='mdl_X' or 'bayesianModel_X')
%
%
%       .fbk.(Ep/Em).(R/P/RP).win_vs_loss
%       (1) win (1) - loss (0) trials
%
%       .fbk.(Ep/Em).(R/P/RP).money_obtained
%       (1) money level obtained at the end of the trial
%
%       .fbk.(Ep/Em).(R/P/RP).E_made
%       (1) level of effort performed during effort performance
%
%       .fbk.(Ep/Em).(R/P/RP).trialN
%       (1) trial number
%       (2) trial number*(E chosen - E non-chosen option)
%       (3) trial number*(E non-default - E default option)
%
%       .fbk.(Ep/Em).(R/P/RP).confidence
%       (1) confidence level (0/1) given by the subject for each choice
%       (2) confidence inferred by the model (p(choice)-0.5)² for the model
%       defined in .fbk.(Ep.Em).(R/P/RP).conf_mdl (='mdl_X' or 'bayesianModel_X')
%
% See also GLM_details.m
%
% Nicolas Clairis - 2021/2022


%% initialize all variables

% general parameters: derivative, grey matter filtering
[GLMprm.gal.add_drv,...
    GLMprm.gal.grey_mask,...
    GLMprm.gal.zPerRun,...
    GLMprm.gal.orth_vars] = deal(0);

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
    RPconditions = {'R','P','RP'};
    for iRP = 1:length(RPconditions)
        RP_nm = RPconditions{iRP};
        
        % pool reward and punishment together (default)
        GLMprm.choice.(EpEm_nm).RPpool = 1;
        [GLMprm.choice.(EpEm_nm).(RP_nm).R_vs_P,...
            GLMprm.choice.(EpEm_nm).(RP_nm).money_left,...
            GLMprm.choice.(EpEm_nm).(RP_nm).money_right,...
            GLMprm.choice.(EpEm_nm).(RP_nm).money_chosen,...
            GLMprm.choice.(EpEm_nm).(RP_nm).money_unchosen,...
            GLMprm.choice.(EpEm_nm).(RP_nm).money_varOption,...
            GLMprm.choice.(EpEm_nm).(RP_nm).money_ch_min_unch,...
            GLMprm.choice.(EpEm_nm).(RP_nm).money_sum,...
            GLMprm.choice.(EpEm_nm).(RP_nm).E_left,...
            GLMprm.choice.(EpEm_nm).(RP_nm).E_right,...
            GLMprm.choice.(EpEm_nm).(RP_nm).E_chosen,...
            GLMprm.choice.(EpEm_nm).(RP_nm).E_unchosen,...
            GLMprm.choice.(EpEm_nm).(RP_nm).E_varOption,...
            GLMprm.choice.(EpEm_nm).(RP_nm).E_ch_min_unch,...
            GLMprm.choice.(EpEm_nm).(RP_nm).E_sum,...
            GLMprm.choice.(EpEm_nm).(RP_nm).NV_chosen,...
            GLMprm.choice.(EpEm_nm).(RP_nm).NV_varOption,...
            GLMprm.choice.(EpEm_nm).(RP_nm).confidence,...
            GLMprm.choice.(EpEm_nm).(RP_nm).RT,...
            GLMprm.choice.(EpEm_nm).(RP_nm).trialN] = deal(0);
        [GLMprm.choice.(EpEm_nm).(RP_nm).NV_mdl,...
            GLMprm.choice.(EpEm_nm).(RP_nm).conf_mdl] = deal('');
        
        % chosen option display
        % pool reward and punishment together (default)
        GLMprm.chosen.(EpEm_nm).RPpool = 1;
        [GLMprm.chosen.(EpEm_nm).(RP_nm).R_vs_P,...
            GLMprm.chosen.(EpEm_nm).(RP_nm).money_chosen,...
            GLMprm.chosen.(EpEm_nm).(RP_nm).money_unchosen,...
            GLMprm.chosen.(EpEm_nm).(RP_nm).money_varOption,...
            GLMprm.chosen.(EpEm_nm).(RP_nm).money_ch_min_unch,...
            GLMprm.chosen.(EpEm_nm).(RP_nm).money_sum,...
            GLMprm.chosen.(EpEm_nm).(RP_nm).E_chosen,...
            GLMprm.chosen.(EpEm_nm).(RP_nm).E_unchosen,...
            GLMprm.chosen.(EpEm_nm).(RP_nm).E_varOption,...
            GLMprm.chosen.(EpEm_nm).(RP_nm).E_ch_min_unch,...
            GLMprm.chosen.(EpEm_nm).(RP_nm).E_sum,...
            GLMprm.chosen.(EpEm_nm).(RP_nm).NV_chosen,...
            GLMprm.chosen.(EpEm_nm).(RP_nm).NV_varOption,...
            GLMprm.chosen.(EpEm_nm).(RP_nm).confidence,...
            GLMprm.chosen.(EpEm_nm).(RP_nm).RT,...
            GLMprm.chosen.(EpEm_nm).(RP_nm).trialN] = deal(0);
        [GLMprm.chosen.(EpEm_nm).(RP_nm).NV_mdl,...
            GLMprm.chosen.(EpEm_nm).(RP_nm).conf_mdl] = deal('');
        
        % pre-effort black cross
        % pool reward and punishment together (default)
        GLMprm.preEffortCross.(EpEm_nm).RPpool = 1;
        [GLMprm.preEffortCross.(EpEm_nm).(RP_nm).money_chosen,...
            GLMprm.preEffortCross.(EpEm_nm).(RP_nm).E_chosen,...
            GLMprm.preEffortCross.(EpEm_nm).(RP_nm).NV_chosen,...
            GLMprm.preEffortCross.(EpEm_nm).(RP_nm).NV_varOption,...
            GLMprm.preEffortCross.(EpEm_nm).(RP_nm).RT_1stAnswer,...
            GLMprm.preEffortCross.(EpEm_nm).(RP_nm).trialN,...
            GLMprm.preEffortCross.(EpEm_nm).(RP_nm).confidence] = deal(0);
        % specific variables for each effort type
        switch EpEm_nm
            case 'Ep'
                % effort performance
                [GLMprm.preEffortCross.Ep.(RP_nm).F_peak,...
                    GLMprm.preEffortCross.Ep.(RP_nm).F_integral] = deal(0);
            case 'Em'
                % effort performance
                [GLMprm.preEffortCross.Em.(RP_nm).RT_avg,...
                    GLMprm.preEffortCross.Em.(RP_nm).n_errors] = deal(0);
        end
        [GLMprm.preEffortCross.(EpEm_nm).(RP_nm).NV_mdl,...
            GLMprm.preEffortCross.(EpEm_nm).(RP_nm).conf_mdl] = deal('');
        
        % effort performance
        % pool reward and punishment together (default)
        GLMprm.Eperf.(EpEm_nm).RPpool = 1;
        [GLMprm.Eperf.(EpEm_nm).(RP_nm).money_chosen,...
            GLMprm.Eperf.(EpEm_nm).(RP_nm).E_chosen,...
            GLMprm.Eperf.(EpEm_nm).(RP_nm).NV_chosen,...
            GLMprm.Eperf.(EpEm_nm).(RP_nm).NV_varOption,...
            GLMprm.Eperf.(EpEm_nm).(RP_nm).RT_1stAnswer,...
            GLMprm.Eperf.(EpEm_nm).(RP_nm).trialN,...
            GLMprm.Eperf.(EpEm_nm).(RP_nm).confidence] = deal(0);
        % specific variables for each effort type
        switch EpEm_nm
            case 'Ep'
                % effort performance
                [GLMprm.Eperf.Ep.(RP_nm).F_peak,...
                    GLMprm.Eperf.Ep.(RP_nm).F_integral] = deal(0);
            case 'Em'
                % effort performance
                [GLMprm.Eperf.Em.(RP_nm).RT_avg,...
                    GLMprm.Eperf.Em.(RP_nm).n_errors] = deal(0);
        end
        [GLMprm.Eperf.(EpEm_nm).(RP_nm).NV_mdl,...
            GLMprm.Eperf.(EpEm_nm).(RP_nm).conf_mdl] = deal('');
        
        % feedback
        % pool reward and punishment together (default)
        GLMprm.fbk.(EpEm_nm).RPpool = 1;
        [GLMprm.fbk.(EpEm_nm).(RP_nm).money_obtained,...
            GLMprm.fbk.(EpEm_nm).(RP_nm).win_vs_loss,...
            GLMprm.fbk.(EpEm_nm).(RP_nm).E_made,...
            GLMprm.fbk.(EpEm_nm).(RP_nm).confidence,...
            GLMprm.fbk.(EpEm_nm).(RP_nm).trialN] = deal(0);
        GLMprm.fbk.(EpEm_nm).(RP_nm).conf_mdl = '';
    end % RP
end % effort type
    
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
        GLMprm.choice.Ep.RP.money_sum = 1;
        GLMprm.choice.Em.RP.money_sum = 1;
        GLMprm.choice.Ep.RP.E_sum = 1;
        GLMprm.choice.Em.RP.E_sum = 1;
        GLMprm.choice.Ep.RP.RT = 1;
        GLMprm.choice.Em.RP.RT = 1;
        % disp chosen
        GLMprm.model_onset.Ep.chosen = 'stick';
        GLMprm.model_onset.Em.chosen = 'stick';
        GLMprm.chosen.Ep.RP.money_chosen = 1;
        GLMprm.chosen.Em.RP.money_chosen = 1;
        GLMprm.chosen.Ep.RP.E_chosen = 1;
        GLMprm.chosen.Em.RP.E_chosen = 1;
        % effort perf
        GLMprm.model_onset.Ep.Eperf = 'stick';
        GLMprm.model_onset.Em.Eperf = 'stick';
        % feedback
        GLMprm.model_onset.Ep.fbk = 'stick';
        GLMprm.model_onset.Em.fbk = 'stick';
        GLMprm.fbk.Ep.RP.money_obtained = 1;
        GLMprm.fbk.Em.RP.money_obtained = 1;
    case 2
        % general parameters
        GLMprm.gal.orth_vars = 1;
        % cross
        GLMprm.model_onset.Ep.preChoiceCross = 'stick';
        GLMprm.model_onset.Em.preChoiceCross = 'stick';
        % choice
        GLMprm.model_onset.Ep.choice = 'stick';
        GLMprm.model_onset.Em.choice = 'stick';
        GLMprm.choice.Ep.RP.RT = 1;
        GLMprm.choice.Em.RP.RT = 1;
        % disp chosen
        GLMprm.model_onset.Ep.chosen = 'stick';
        GLMprm.model_onset.Em.chosen = 'stick';
        GLMprm.chosen.Ep.RP.money_chosen = 1;
        GLMprm.chosen.Em.RP.money_chosen = 1;
        GLMprm.chosen.Ep.RP.E_chosen = 1;
        GLMprm.chosen.Em.RP.E_chosen = 1;
        % effort perf
        GLMprm.model_onset.Ep.Eperf = 'stick';
        GLMprm.model_onset.Em.Eperf = 'stick';
        % feedback
        GLMprm.model_onset.Ep.fbk = 'stick';
        GLMprm.model_onset.Em.fbk = 'stick';
        GLMprm.fbk.Ep.RP.win_vs_loss = 1;
        GLMprm.fbk.Em.RP.win_vs_loss = 1;
    case 3
        % general parameters
        GLMprm.gal.orth_vars = 1;
        % cross
        GLMprm.model_onset.Ep.preChoiceCross = 'stick';
        GLMprm.model_onset.Em.preChoiceCross = 'stick';
        % choice
        GLMprm.model_onset.Ep.choice = 'stick';
        GLMprm.model_onset.Em.choice = 'stick';
        GLMprm.choice.Ep.RP.R_vs_P = 1;
        GLMprm.choice.Em.RP.R_vs_P = 1;
        GLMprm.choice.Ep.RP.RT = 1;
        GLMprm.choice.Em.RP.RT = 1;
        % disp chosen
        GLMprm.model_onset.Ep.chosen = 'stick';
        GLMprm.model_onset.Em.chosen = 'stick';
        GLMprm.chosen.Ep.RP.money_chosen = 1;
        GLMprm.chosen.Em.RP.money_chosen = 1;
        GLMprm.chosen.Ep.RP.E_chosen = 1;
        GLMprm.chosen.Em.RP.E_chosen = 1;
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
            GLMprm.choice.(Epm_nm).RP.R_vs_P = 1;
            GLMprm.choice.(Epm_nm).RP.money_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E_chosen = 1;
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
            GLMprm.choice.(Epm_nm).RP.R_vs_P = 1;
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
            GLMprm.choice.(Epm_nm).RP.money_ch_min_unch = 1;
            GLMprm.choice.(Epm_nm).RP.E_ch_min_unch = 1;
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
                GLMprm.choice.(Epm_nm).(RP_nm).money_chosen = 3;
                GLMprm.choice.(Epm_nm).(RP_nm).money_unchosen = 3;
                GLMprm.choice.(Epm_nm).(RP_nm).E_chosen = 1;
                GLMprm.choice.(Epm_nm).(RP_nm).E_unchosen = 1;
                GLMprm.choice.(Epm_nm).(RP_nm).RT = 1;
            end
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort performance
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback - split R/P
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
            GLMprm.fbk.(Epm_nm).RP.win_vs_loss = 1;
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
            GLMprm.choice.(Epm_nm).RP.RT = 1;
            % chosen - split R/P and use R/P and E levels
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            GLMprm.chosen.(Epm_nm).RPpool = 0;
            for iRP = 1:length(RP_conds)
                RP_nm = RP_conds{iRP};
                GLMprm.chosen.(Epm_nm).(RP_nm).money_chosen = 3;
                GLMprm.chosen.(Epm_nm).(RP_nm).money_unchosen = 3;
                GLMprm.chosen.(Epm_nm).(RP_nm).E_chosen = 1;
                GLMprm.chosen.(Epm_nm).(RP_nm).E_unchosen = 1;
            end
            % effort performance
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback - split R/P
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
            GLMprm.fbk.(Epm_nm).RP.win_vs_loss = 1;
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
                GLMprm.choice.(Epm_nm).(RP_nm).money_varOption = 3;
                GLMprm.choice.(Epm_nm).(RP_nm).E_varOption = 1;
                GLMprm.choice.(Epm_nm).(RP_nm).RT = 1;
            end
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort performance
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback - split R/P
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
            GLMprm.fbk.(Epm_nm).RP.win_vs_loss = 1;
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
                GLMprm.Eperf.(Epm_nm).(RP_nm).money_chosen = 3;
                GLMprm.Eperf.(Epm_nm).(RP_nm).E_chosen = 1;
            end
            % feedback - split R/P
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
            GLMprm.fbk.(Epm_nm).RP.win_vs_loss = 1;
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
            GLMprm.choice.(Epm_nm).RP.RT = 1;
            % chosen - split R/P and use R/P and E levels
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            GLMprm.chosen.(Epm_nm).RPpool = 0;
            for iRP = 1:length(RP_conds)
                RP_nm = RP_conds{iRP};
                GLMprm.chosen.(Epm_nm).(RP_nm).money_chosen = 3;
                GLMprm.chosen.(Epm_nm).(RP_nm).E_chosen = 1;
            end
            % effort performance
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback - split R/P
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
            GLMprm.fbk.(Epm_nm).RP.win_vs_loss = 1;
        end
    case 12 % Vch/R-P/VE during dispChosen option
        GLMprm.gal.orth_vars = 1;
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.RT = 1;
            % chosen - split R/P and use R/P and E levels
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            GLMprm.chosen.(Epm_nm).RPpool = 1;
            GLMprm.chosen.(Epm_nm).RP.money_chosen = 1;
            GLMprm.chosen.(Epm_nm).RP.E_chosen = 1;
            GLMprm.chosen.(Epm_nm).RP.R_vs_P = 1;
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
            GLMprm.choice.(Epm_nm).RP.RT = 1;
            % chosen - split R/P and use R/P and E levels
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            GLMprm.chosen.(Epm_nm).RPpool = 0;
            for iRP = 1:length(RP_conds)
                RP_nm = RP_conds{iRP};
                GLMprm.chosen.(Epm_nm).(RP_nm).money_chosen = 3;
                GLMprm.chosen.(Epm_nm).(RP_nm).E_chosen = 1;
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
            GLMprm.choice.(Epm_nm).RP.RT = 1;
            % chosen - split R/P and use R/P and E levels
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            GLMprm.chosen.(Epm_nm).RPpool = 0;
            for iRP = 1:length(RP_conds)
                RP_nm = RP_conds{iRP};
                GLMprm.chosen.(Epm_nm).(RP_nm).money_chosen = 3;
                GLMprm.chosen.(Epm_nm).(RP_nm).E_chosen = 1;
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
            GLMprm.choice.(Epm_nm).RP.R_vs_P = 1;
            GLMprm.choice.(Epm_nm).RP.money_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E_chosen = 1;
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
            GLMprm.choice.(Epm_nm).RP.NV_varOption = 1;
            GLMprm.choice.(Epm_nm).RP.NV_mdl = 'mdl_2';
            GLMprm.choice.(Epm_nm).RP.RT = 1;
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
                GLMprm.choice.(Epm_nm).(RP_nm).NV_varOption = 1;
                GLMprm.choice.(Epm_nm).(RP_nm).NV_mdl = 'mdl_2';
                GLMprm.choice.(Epm_nm).(RP_nm).RT = 1;
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
                GLMprm.chosen.(Epm_nm).(RP_nm).NV_chosen = 1;
                GLMprm.chosen.(Epm_nm).(RP_nm).NV_mdl = 'mdl_2';
                GLMprm.chosen.(Epm_nm).(RP_nm).RT = 1;
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
                GLMprm.choice.(Epm_nm).(RP_nm).money_varOption = 3;
                GLMprm.choice.(Epm_nm).(RP_nm).E_varOption = 1;
                GLMprm.choice.(Epm_nm).(RP_nm).trialN = 3;
                GLMprm.choice.(Epm_nm).(RP_nm).RT = 1;
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
                GLMprm.chosen.(Epm_nm).(RP_nm).money_varOption = 3;
                GLMprm.chosen.(Epm_nm).(RP_nm).E_varOption = 1;
                GLMprm.chosen.(Epm_nm).(RP_nm).trialN = 3;
                GLMprm.chosen.(Epm_nm).(RP_nm).RT = 1;
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
                GLMprm.Eperf.(Epm_nm).(RP_nm).NV_chosen = 1;
                GLMprm.Eperf.(Epm_nm).(RP_nm).NV_mdl = 'mdl_4';
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
                GLMprm.choice.(Epm_nm).(RP_nm).money_chosen = 1;
                GLMprm.choice.(Epm_nm).(RP_nm).E_chosen = 1;
                GLMprm.choice.(Epm_nm).(RP_nm).trialN = 1;
                GLMprm.choice.(Epm_nm).(RP_nm).RT = 1;
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
                GLMprm.chosen.(Epm_nm).(RP_nm).money_chosen = 1;
                GLMprm.chosen.(Epm_nm).(RP_nm).E_chosen = 1;
                GLMprm.chosen.(Epm_nm).(RP_nm).trialN = 1;
                GLMprm.chosen.(Epm_nm).(RP_nm).RT = 1;
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
                GLMprm.Eperf.(Epm_nm).(RP_nm).money_chosen = 1;
                GLMprm.Eperf.(Epm_nm).(RP_nm).E_chosen = 1;
                GLMprm.Eperf.(Epm_nm).(RP_nm).trialN = 1;
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
                GLMprm.choice.(Epm_nm).(RP_nm).NV_chosen = 1;
                GLMprm.choice.(Epm_nm).(RP_nm).NV_mdl = 'mdl_4';
                GLMprm.choice.(Epm_nm).(RP_nm).RT = 1;
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
                GLMprm.chosen.(Epm_nm).(RP_nm).NV_chosen = 1;
                GLMprm.chosen.(Epm_nm).(RP_nm).NV_mdl = 'mdl_4';
                GLMprm.chosen.(Epm_nm).(RP_nm).RT = 1;
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
                GLMprm.Eperf.(Epm_nm).(RP_nm).NV_chosen = 1;
                GLMprm.Eperf.(Epm_nm).(RP_nm).NV_mdl = 'mdl_4';
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
                GLMprm.choice.(Epm_nm).(RP_nm).money_varOption = 1;
                GLMprm.choice.(Epm_nm).(RP_nm).E_varOption = 1;
                GLMprm.choice.(Epm_nm).(RP_nm).trialN = 1;
                GLMprm.choice.(Epm_nm).(RP_nm).RT = 1;
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
                GLMprm.chosen.(Epm_nm).(RP_nm).money_varOption = 1;
                GLMprm.chosen.(Epm_nm).(RP_nm).E_varOption = 1;
                GLMprm.chosen.(Epm_nm).(RP_nm).trialN = 1;
                GLMprm.chosen.(Epm_nm).(RP_nm).RT = 1;
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
                GLMprm.Eperf.(Epm_nm).(RP_nm).money_varOption = 1;
                GLMprm.Eperf.(Epm_nm).(RP_nm).E_varOption = 1;
                GLMprm.Eperf.(Epm_nm).(RP_nm).trialN = 1;
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
                GLMprm.choice.(Epm_nm).(RP_nm).NV_varOption = 1;
                GLMprm.choice.(Epm_nm).(RP_nm).NV_mdl = 'mdl_4';
                GLMprm.choice.(Epm_nm).(RP_nm).RT = 1;
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
                GLMprm.chosen.(Epm_nm).(RP_nm).NV_varOption = 1;
                GLMprm.chosen.(Epm_nm).(RP_nm).NV_mdl = 'mdl_4';
                GLMprm.chosen.(Epm_nm).(RP_nm).RT = 1;
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
                GLMprm.Eperf.(Epm_nm).(RP_nm).NV_varOption = 1;
                GLMprm.Eperf.(Epm_nm).(RP_nm).NV_mdl = 'mdl_4';
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
            GLMprm.choice.(Epm_nm).RP.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            GLMprm.chosen.(Epm_nm).RPpool = 1;
            GLMprm.chosen.(Epm_nm).RP.money_chosen = 1;
            GLMprm.chosen.(Epm_nm).RP.E_chosen = 1;
            GLMprm.chosen.(Epm_nm).RP.trialN = 1;
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
            GLMprm.choice.(Epm_nm).RP.money_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.trialN = 1;
            GLMprm.choice.(Epm_nm).RP.RT = 1;
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
            GLMprm.choice.(Epm_nm).RP.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            GLMprm.chosen.(Epm_nm).RPpool = 1;
            GLMprm.chosen.(Epm_nm).RP.NV_chosen = 1;
            GLMprm.chosen.(Epm_nm).RP.NV_mdl = 'mdl_4';
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
            GLMprm.choice.(Epm_nm).RP.NV_varOption = 1;
            GLMprm.choice.(Epm_nm).RP.NV_mdl = 'mdl_4';
            GLMprm.choice.(Epm_nm).RP.confidence = 2;
            GLMprm.choice.(Epm_nm).RP.conf_mdl = 'mdl_4';
            GLMprm.choice.(Epm_nm).RP.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            GLMprm.chosen.(Epm_nm).RPpool = 1;
            GLMprm.chosen.(Epm_nm).RP.NV_chosen = 1;
            GLMprm.chosen.(Epm_nm).RP.NV_mdl = 'mdl_4';
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
            GLMprm.choice.(Epm_nm).RP.confidence = 2;
            GLMprm.choice.(Epm_nm).RP.conf_mdl = 'mdl_4';
            GLMprm.choice.(Epm_nm).RP.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'boxcar_bis';
            GLMprm.chosen.(Epm_nm).RPpool = 1;
            GLMprm.chosen.(Epm_nm).RP.money_chosen = 3;
            GLMprm.chosen.(Epm_nm).RP.E_chosen = 1;
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
            GLMprm.choice.(Epm_nm).RP.confidence = 2;
            GLMprm.choice.(Epm_nm).RP.conf_mdl = 'mdl_4';
            GLMprm.choice.(Epm_nm).RP.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % pre-effort cross (effort preparation)
            GLMprm.model_onset.(Epm_nm).preEffortCross = 'boxcar_bis';
            GLMprm.preEffortCross.(Epm_nm).RPpool = 1;
            GLMprm.preEffortCross.(Epm_nm).RP.money_chosen = 3;
            GLMprm.preEffortCross.(Epm_nm).RP.E_chosen = 1;
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
end % GLM number
%% warnings: check compatibility of the GLM parameters entered
isGLMokCheck(GLMprm);

end % function