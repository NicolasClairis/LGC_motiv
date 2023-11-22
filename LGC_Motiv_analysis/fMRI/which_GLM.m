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
%       (0) no manual mask, only SPM implicit masking
%       (1) grey-matter filter based on individual grey matter
%       (2) use SPM template grey matter mask
%       (3) grey-matter filter based on group average
%       (4) SPM implicit mask but with a lower threshold (0.5 instead of 0.8)
%       (5) SPM implicit mask but with a lower threshold (0.3 instead of 0.8)
%       (6) SPM implicit mask but with a lower threshold (0.1 instead of 0.8)
%       (7) grey + white matter filter based on group average
%
%       .mask_probaThreshold:
%       value indicating probability to be in grey(+white) matter to
%       consider in case .grey_mask = 1, 2, 3 or 7
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
%       .autocorrel: which method to use for 1st level temporal
%       autocorrelation of fMRI time series
%       (0) default AR algorithm for temporal autocorrelation
%       (1) FAST algorithm. Better suited for short TR and generally said
%       to bring more accurate results (see (Olszowy et al, 2019)) but
%       first test (GLM182) not super conclusive. Some maps look the same
%       and others just look uglier. Probably better suited for short TR.
%
%   .model_onset: indicate for each task (Ep/Em: physical/mental) for each
%   event (preChoiceCross/choice/chosen/preEffortCross/Eperf/fbk) if it should be modelled as a
%   stick ('stick') as a boxcar ('boxcar') or not included in the GLM
%   ('none'), for some cases, 'boxcar_bis' and 'boxcar_ter' correspond to 
%   situations where the boxcar also entails following periods of the task
%       .allCrosses: pool of both pre-choice and pre-effort fixation crosses
%       periods (no regressors included)
%       .preChoiceCross: (white) fixation cross before choice period
%           - 'boxcar_bis': entails fixation cross and display of the two
%           options until the choice is done
%       .choice: choice period (when options are displayed on screen)
%           - 'boxcar_bis': entails choice + chosen periods
%       .chosen: moment when the chosen option is displayed on screen
%           - 'boxcar_bis': entails chosen option display until the end of
%           the effort exertion period (effort preparation and execution period)
%           - 'boxcar_ter': entails chosen option display until the end of
%           the fixation cross related to effort preparation (effort
%           preparation period)
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
%       (2) (-1) when low effort/default option is selected and (+1) when high
%       effort/non-default option is selected
%
%   .choice/chosen/preEffortCross/Eperf/fbk: for each event, for each task (Ep/Em) and 
%   for each condition (R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch) indicate if a given regressor should be 
%   included or not.
%       .Ep/Em: physical (Ep) or mental (Em) effort task
%           .R/P/RP: reward only (R), punishment only (P) or reward and punishment
%           trials mixed (RP)
%           .splitPerE/splitPerEch: split per effort proposed (E1/E2/E3) (splitPerE=1)
%           or per effort chosen (Ech0/Ech/Ech1/Ech2/Ech3) (splitPerE=2)
%           or per choice (low vs high E choice) made (lEch/hEch) (splitPerE=3)
%               
%% choice/chosen period regressors
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
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).R_vs_P:
%       (1) 0 when punishment trial, 1 when reward trial
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).choiceHighE
%       (1) 0 when low effort/default option is selected and 1 when high
%       effort/non-default option is selected
%       (2) (-1) when low effort/default option is selected and (+1) when high
%       effort/non-default option is selected
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).R_varOption
%       (1) reward amount for the high effort option (zeros for punishment)
%       (2) reward level for the high effort option (zeros for punishment)
%       (3) zscored reward amount for the high effort option (zeros for punishment)
%       (4) zscored reward level for the high effort option (zeros for punishment)
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).R_chosen
%       (1) reward amount for the chosen option (zeros for punishment)
%       (2) reward level for the chosen option (zeros for punishment)
%       (3) zscored reward amount for the chosen option (zeros for punishment)
%       (4) zscored reward level for the chosen option (zeros for punishment)
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).P_varOption
%       (1) punishment amount for the high effort option (zeros for rewards)
%       (2) punishment level for the high effort option (zeros for rewards)
%       (3) zscored punishment amount for the high effort option (zeros for rewards)
%       (4) zscored punishment level for the high effort option (zeros for rewards)
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).P_chosen
%       (1) punishment amount for the chosen option (zeros for rewards)
%       (2) punishment level for the chosen option (zeros for rewards)
%       (3) zscored punishment amount for the chosen option (zeros for rewards)
%       (4) zscored punishment level for the chosen option (zeros for rewards)
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).money_left
%       (1) money amount associated to left option
%       (2) |money amount| associated to left option
%       (3) money level associated to left option
%       (4) |money level| associated to left option
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).money_right
%       (1) money amount associated to right option
%       (2) |money amount| associated to right option
%       (3) money level associated to right option
%       (4) |money level| associated to right option
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).money_chosen:
%       (1) money chosen amount
%       (2) |money chosen amount|
%       (3) money levels (-4/-3/-2/-1/1/2/3/4) chosen from highest loss
%       until highest gain
%       (4) |money levels (1/2/3/4) chosen|
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).money_unchosen:
%       (1) money unchosen amount
%       (2) |money unchosen amount|
%       (3) money levels (-4/-3/-2/-1/1/2/3/4) unchosen
%       (4) |money levels (1/2/3/4) unchosen|
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).money_varOption:
%       (1) money amount non-default option
%       (2) |money amount non-default option|
%       (3) money levels  (-3/-2/-1/1/2/3) non-default option
%       (4) |money levels (1/2/3) non-default option|
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).money_ch_min_unch:
%       (1) money chosen - money unchosen amount
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).money_ch_min_fixOption:
%       (1) money chosen - money default option - amount
%       (2) money chosen - money default option - levels (1/2/3/4 - 1)
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).money_sum:
%       (1) money default + money non-default amount
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).E_left
%       (1) effort level (0/1/2/3) associated to left option
%       (2) effort difficulty associated to left option (Ep: duration to hold; Em: nb answers to give)
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).E_right
%       (1) effort level (0/1/2/3) associated to right option
%       (2) effort difficulty associated to right option (Ep: duration to hold; Em: nb answers to give)
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).E_chosen
%       (1) effort level (0/1/2/3) associated to chosen option
%       (2) effort difficulty associated to chosen option (Ep: duration to hold; Em: nb answers to give)
%       (3) effort level of the high effort option (1/2/3)*choice (-1/+1)
%       => Ech between (-3/-2/-1/1/2/3)
%       (4) zscored effort level (0/1/2/3) associated to chosen option
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).E_unchosen
%       (1) effort level (0/1/2/3) associated to unchosen option
%       (2) effort difficulty associated to unchosen option (Ep: duration to hold; Em: nb answers to give)
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).E_varOption
%       (1) effort level (1/2/3) associated to the non-default option
%       (2) effort difficulty associated to the non-default option (Ep: duration to hold; Em: nb answers to give)
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).E_ch_min_unch
%       (1) effort level (0/1/2/3) associated to the chosen option minus
%       the effort level (0/1/2/3) associated to the unchosen option
%       (2) effort difficulty associated to the chosen option (Ep: duration to hold; Em: nb answers to give)
%       minus the effort difficulty associated to the unchosen option (Ep: duration to hold; Em: nb answers to give)
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).E_ch_min_fixOption:
%       (1) E chosen - E default option - levels (0/1/2/3)
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).E_sum
%       (1) sum of the effort levels (0/1/2/3) associated to both options
%       (2) sum of the effort difficulties (Ep: duration to hold; Em: nb answers to give) associated to both options
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).money_level_x_E_varOption
%       (1) (P from 1 to 4 and R from 5 to 8)*(E from 1 to 3) for high effort option
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).money_level_x_E_chosen
%       (1) (P from 1 to 4 and R from 5 to 8)*(E from 1 to 3) for chosen option
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).R_level_x_E_varOption
%       (1) (R 1 to 4)*(E 1 to 4) for high E option (0 for punishment trials)
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).R_level_x_E_chosen
%       (1) (R 1 to 4)*(E 1 to 4) for chosen option (0 for punishment trials)
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).P_level_x_E_varOption
%       (1) (P 1 to 4)*(E 1 to 4) for high E option (0 for reward trials)
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).P_level_x_E_chosen
%       (1) (P 1 to 4)*(E 1 to 4) for chosen option (0 for reward trials)
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_chosen
%       (1) net value of the chosen-unchosen option based on the model defined in
%       .(choice/chosen).(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_mdl (='mdl_X' or 'bayesianModel_X')
%       (2) p(choice = hE) = probability of choosing the chosen option
%       based on the model defined in .(choice/chosen).(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_mdl (='mdl_X' or 'bayesianModel_X')
%       This is relatively similar to (1) but will include the motivational
%       bias and the sigmoid transformation
%       (3) net value of the chosen-unchosen option based on the model defined in
%       .(choice/chosen).(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_mdl (='mdl_X' or 'bayesianModel_X')
%       Same as (1) but will also include the bias
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_varOption
%       (1) difference between net value of the high effort option and low effort options based on the model defined in
%       .(choice/chosen).(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_mdl (='mdl_X' or 'bayesianModel_X')
%       (2) absolute difference between net value of the high effort option
%       and low effort option based on the model defined in .(choice/chosen).(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_mdl (='mdl_X' or 'bayesianModel_X')
%       (3) p(choice = hE) = probability of choosing the high effort option
%       based on the model defined in .(choice/chosen).(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_mdl (='mdl_X' or 'bayesianModel_X')
%       This is relatively similar to (2) but will include the motivational
%       bias and the sigmoid transformation
%       (4) difference between net value of the high effort option
%       and low effort option based on the model defined in .(choice/chosen).(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_mdl (='mdl_X' or 'bayesianModel_X')
%       This is basically like (3) but before the sigmoid transformation
%       (5) absolute difference between net value of the high effort option
%       and low effort option based on the model defined in .(choice/chosen).(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_mdl (='mdl_X' or 'bayesianModel_X')
%       This is basically like (2) but including the bias.
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_varOption_bis
%       (1) difference between net value of the high effort option and low effort options based on the model defined in
%       .(choice/chosen).(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_mdl (='mdl_X' or 'bayesianModel_X')
%       (2) absolute difference between net value of the high effort option
%       and low effort option based on the model defined in .(choice/chosen).(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_mdl (='mdl_X' or 'bayesianModel_X')
%       (3) p(choice = hE) = probability of choosing the high effort option
%       based on the model defined in .(choice/chosen).(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_mdl (='mdl_X' or 'bayesianModel_X')
%       This is relatively similar to (2) but will include the motivational
%       bias and the sigmoid transformation
%       (4) difference between net value of the high effort option
%       and low effort option based on the model defined in .(choice/chosen).(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_mdl (='mdl_X' or 'bayesianModel_X')
%       This is basically like (3) but before the sigmoid transformation
%       (5) absolute difference between net value of the high effort option
%       and low effort option based on the model defined in .(choice/chosen).(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_mdl (='mdl_X' or 'bayesianModel_X')
%       This is basically like (2) but including the bias.
%
%       .(choice/chosen).Ep.(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).F_integral: physical effort only: 
%       force integral for the performance of the current trial
%       (=anticipated force produced)
%       (1) integral of effort performed during the effort period
%       (2) integral of effort performed during the effort period
%       considering only the overshoot (force produced above the red
%       threshold)
%       (3) integral of effort performed during the effort period with
%       force in newtons
%       (4) integral of effort performed during the effort period
%       considering only the overshoot (force produced above the red
%       threshold) with force in newtons
%       (5-8) same as (1-4) but zscored
%
%       .(choice/chosen).Em.(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).efficacy:
%       mental effort only: efficacy for the performance of the current trial
%       (=anticipated effort produced) defined as:
%       (1) (nb correct answers - nb errors made)/total time of effort
%       period
%       (2) (nb correct answers - nb errors made)/time spent between answer
%       of second useless number and end of the trial
%       (3) nb correct answers/total time of effort period
%       (4) nb correct answers/time spent between answer of second useless 
%       number and end of the trial
%       (5-8) same as (1-4) but zscored
%
%       .(choice/chosen).Ep.(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).fatigue:
%       physical effort only:
%       (1) sum of previous AUC of force produced
%       (2) same as (1) but zscored)
%
%       .(choice/chosen).Em.(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).prevEfficacy:
%       mental effort only:
%       (1) efficacy of previous trial computed as (nb correct answers - nb errors made)/total time of effort
%       period
%       (2) efficacy of previous trial computed  as (nb correct answers - nb errors made)/time spent between answer
%       of second useless number and end of the trial
%       (3) efficacy of previous trial computed as (nb correct
%       answers)/total time of effort period (no errors)
%       (4) efficacy of previous trial computed as (nb correct
%       answers)/total time of effort period (no errors) between answer
%       of second useless number and end of the trial
%       (5-8) same as (1-4) but zscored
%
%       .(choice/chosen).Ep.(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).Ech_x_fatigue:
%       physical effort only:
%       (1) (Effort chosen) x (sum of previous AUC of force produced)
%
%       .(choice/chosen).Em.(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).Ech_x_prevEfficacy:
%       mental effort only:
%       (1) (Effort chosen) x (efficacy of previous trial computed as (nb correct answers - nb errors made)/total time of effort
%       period)
%       (2) (Effort chosen) x (efficacy of previous trial computed  as (nb correct answers - nb errors made)/time spent between answer
%       of second useless number and end of the trial)
%       (3) (Effort chosen) x (efficacy of previous trial computed as (nb correct
%       answers)/total time of effort period (no errors))
%       (4) (Effort chosen) x (efficacy of previous trial computed as (nb correct
%       answers)/total time of effort period (no errors) between answer
%       of second useless number and end of the trial)
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).trialN
%       (1) trial number
%       (2) (trial number)*(E chosen - E non-chosen option)
%       (3) (trial number)*(E non-default - E default option)
%       (4) (trial number)*(E non-default)
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).confidence
%       (1) confidence level (0/1) given by the subject for each choice
%       (2) confidence inferred by the model (p(high E)-0.5)² for the model
%       defined in .(choice/chosen).(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).conf_mdl (='mdl_X' or 'bayesianModel_X')
%       (3) confidence inferred by the model (p(left)-0.5)² for the model
%       defined in .(choice/chosen).(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).conf_mdl (='mdl_X' or 'bayesianModel_X')
%       (4) confidence inferred by the model (p(chosen)-0.5)² for the model
%       defined in .(choice/chosen).(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).conf_mdl (='mdl_X' or 'bayesianModel_X')
%
%       .(choice/chosen).(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).RT: reaction time for choice
%       (1) raw reaction time
%       (2) reaction time zscored per run
%       (3) reaction time zscored per subject across all runs
%       (4-6) same as (1-3) but with RT as first regressor instead of last
%
%% pre-effort cross period regressors
%       .preEffortCross.(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).choiceHighE
%       (1) 0 when low effort/default option is selected and 1 when high
%       effort/non-default option is selected
%       (2) (-1) when low effort/default option is selected and (+1) when high
%       effort/non-default option is selected
%
%       .preEffortCross.(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).money_chosen:
%       (1) money chosen amount
%       (2) |money chosen amount|
%       (3) money levels (-4/-3/-2/-1/1/2/3/4) chosen from highest loss
%       until highest gain
%       (4) |money levels (1/2/3/4) chosen|
%
%       .preEffortCross.(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).E_chosen
%       (1) effort level (0/1/2/3) associated to chosen option
%       (2) effort difficulty associated to chosen option (Ep: duration to hold; Em: nb answers to give)
%       (3) effort level of the high effort option (1/2/3)*choice (-1/+1)
%       => Ech between (-3/-2/-1/1/2/3)
%       (4) zscored effort level (0/1/2/3) associated to chosen option
%
%       .preEffortCross.Ep.(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).F_peak: physical effort only: force peak
%       (1) force peak in voltage
%       (2) force peak in newtons
%       (3-4) like (1-2) but zscored
%
%       .preEffortCross.Ep.(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).F_integral: physical effort only: force integral
%       (1) integral of effort performed during the effort period
%       (2) integral of effort performed during the effort period
%       considering only the overshoot (force produced above the red
%       threshold)
%       (3) integral of effort performed during the effort period with
%       force in newtons
%       (4) integral of effort performed during the effort period
%       considering only the overshoot (force produced above the red
%       threshold) with force in newtons
%       (5-8) same as (1-4) but zscored
%
%       .preEffortCross.Em.(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).RT_avg: mental effort only: average reaction
%       time for answering N-back task
%       (1) average reaction time for answering to questions (in seconds)
%
%       .preEffortCross.Em.(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).n_correct:
%       mental effort only: number of correct answers made per trial
%       (1) number of correct answers provided (ignoring the two first
%       useless digits)
%
%       .preEffortCross.Em.(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).n_errors: mental effort only: number of errors
%       made per trial
%       (1) number of errors made
%
%       .preEffortCross.(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_chosen
%       (1) net value of the chosen option based on the model defined in
%       .(choice/chosen).(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_mdl (='mdl_X' or 'bayesianModel_X')
%       (2) p(choice = hE) = probability of choosing the chosen option
%       based on the model defined in .(choice/chosen).(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_mdl (='mdl_X' or 'bayesianModel_X')
%       This is relatively similar to (1) but will include the motivational
%       bias and the sigmoid transformation
%       (3) net value of the chosen-unchosen option based on the model defined in
%       .(preEffortCross).(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_mdl (='mdl_X' or 'bayesianModel_X')
%       Same as (1) but will also include the bias
%
%       .preEffortCross.(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_varOption
%       (1) delta of the net value of the high E - low E based on the model defined in
%       .Eperf.(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_mdl (='mdl_X' or 'bayesianModel_X')
%       (2) absolute difference between net value of the high effort option
%       and low effort option based on the model defined in .(choice/chosen).(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_mdl (='mdl_X' or 'bayesianModel_X')
%       (3) p(choice = hE) = probability of choosing the high effort option
%       based on the model defined in .Eperf.(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_mdl (='mdl_X' or 'bayesianModel_X')
%       This is relatively similar to (2) but will include the motivational
%       bias and the sigmoid transformation
%       (4) difference between net value of the high effort option
%       and low effort option based on the model defined in .(choice/chosen).(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_mdl (='mdl_X' or 'bayesianModel_X')
%       This is basically like (3) but before the sigmoid transformation
%       (5) absolute difference between net value of the high effort option
%       and low effort option based on the model defined in .(choice/chosen).(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_mdl (='mdl_X' or 'bayesianModel_X')
%       This is basically like (2) but including the bias.
%
%       .preEffortCross.(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_varOption_bis
%       (1) difference between net value of the high effort option and low effort options based on the model defined in
%       .(choice/chosen).(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_mdl (='mdl_X' or 'bayesianModel_X')
%       (2) absolute difference between net value of the high effort option
%       and low effort option based on the model defined in .(choice/chosen).(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_mdl (='mdl_X' or 'bayesianModel_X')
%       (3) p(choice = hE) = probability of choosing the high effort option
%       based on the model defined in .(choice/chosen).(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_mdl (='mdl_X' or 'bayesianModel_X')
%       This is relatively similar to (2) but will include the motivational
%       bias and the sigmoid transformation
%       (4) difference between net value of the high effort option
%       and low effort option based on the model defined in .(choice/chosen).(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_mdl (='mdl_X' or 'bayesianModel_X')
%       This is basically like (3) but before the sigmoid transformation
%       (5) absolute difference between net value of the high effort option
%       and low effort option based on the model defined in .(choice/chosen).(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_mdl (='mdl_X' or 'bayesianModel_X')
%       This is basically like (2) but including the bias.
%
%       .preEffortCross.Eperf.(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).RT_1stAnswer
%       (1) raw reaction time for first answer (force above threshold for Ep 
%       and first answer to first digit for Em)
%
%       .preEffortCross.(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).trialN
%       (1) trial number
%       (2) (trial number)*(E chosen - E non-chosen option)
%       (3) (trial number)*(E non-default - E default option)
%       (4) (trial number)*(E non-default)
%
%% Effort period regressors
%       .Eperf.(Ep/Em).RPpool
%       (0) split rewards and punishments as separate events
%       (1) pool reward and punishment trials
%
%       .Eperf.(Ep/Em).splitPerE:
%       (0) pool all effort together (by default)
%       (1) split per level of effort proposed (E1/E2/E3)
%       (2) split per level of effort chosen (Ech0/Ech1/Ech2/Ech3)
%       (3) split per option chosen (low/default vs high/non-default effort chosen)
%
%       .Eperf.(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).choiceHighE
%       (1) 0 when low effort/default option is selected and 1 when high
%       effort/non-default option is selected
%       (2) (-1) when low effort/default option is selected and (+1) when high
%       effort/non-default option is selected
%
%       .Eperf.(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).money_chosen
%       (1) money chosen amount
%       (2) |money chosen amount|
%       (3) money levels (-4/-3/-2/-1/1/2/3/4) chosen
%       (4) |money levels (1/2/3/4) chosen|
%
%       .Eperf.(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).E_chosen
%       (1) effort level (0/1/2/3) associated to chosen option
%       (2) effort difficulty associated to chosen option (Ep: duration to hold; Em: nb answers to give)
%       (3) effort level of the high effort option (1/2/3)*choice (-1/+1)
%       => Ech between (-3/-2/-1/1/2/3)
%       (4) zscored effort level (0/1/2/3) associated to chosen option
%
%       .Eperf.Ep.(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).F_peak: physical effort only: force peak
%       (1) force peak in voltage
%       (2) force peak in newtons
%       (3-4) same as (1-2) but zscored
%
%       .Eperf.Ep.(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).F_integral: physical effort only: force integral
%       (1) integral of effort performed during the effort period
%       (2) integral of effort performed during the effort period
%       considering only the overshoot (force produced above the red
%       threshold)
%       (3) integral of effort performed during the effort period with
%       force in newtons
%       (4) integral of effort performed during the effort period
%       considering only the overshoot (force produced above the red
%       threshold) with force in newtons
%       (5-8) same as (1-4) but zscored
%
%       .Eperf.Em.(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).efficacy:
%       mental effort only: efficacy defined as:
%       (1) (nb correct answers - nb errors made)/total time of effort
%       period
%       (2) (nb correct answers - nb errors made)/time spent between answer
%       of second useless number and end of the trial
%       (3) nb correct answers/total time of effort period
%       (4) nb correct answers/time spent between answer of second useless 
%       number and end of the trial
%       (5-8) same as (1-4) but zscored
%
%       .Eperf.Em.(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).RT_avg: mental effort only: average reaction
%       time for answering N-back task
%       (1) average reaction time for answering to questions (in seconds)
%
%       .Eperf.Em.(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).n_correct:
%       mental effort only: number of correct answers made per trial
%       (1) number of correct answers provided (ignoring the two first
%       useless digits)
%
%       .Eperf.Em.(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).n_errors: mental effort only: number of errors
%       made per trial
%       (1) number of errors made
%
%       .Eperf.(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_chosen
%       (1) net value of the chosen option based on the model defined in
%       .Eperf.(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_mdl (='mdl_X' or 'bayesianModel_X')
%       (2) p(choice = hE) = probability of choosing the chosen option
%       based on the model defined in .Eperf.(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_mdl (='mdl_X' or 'bayesianModel_X')
%       This is relatively similar to (1) but will include the motivational
%       bias and the sigmoid transformation
%       (3) net value of the chosen-unchosen option based on the model defined in
%       .(preEffortCross).(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_mdl (='mdl_X' or 'bayesianModel_X')
%       Same as (1) but will also include the bias
%
%       .Eperf.(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_varOption
%       (1) delta of the net value of the high E - low E based on the model defined in
%       .Eperf.(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_mdl (='mdl_X' or 'bayesianModel_X')
%       (2) absolute difference between net value of the high effort option
%       and low effort option based on the model defined in .(choice/chosen).(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_mdl (='mdl_X' or 'bayesianModel_X')
%       (3) p(choice = hE) = probability of choosing the high effort option
%       based on the model defined in .Eperf.(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_mdl (='mdl_X' or 'bayesianModel_X')
%       This is relatively similar to (2) but will include the motivational
%       bias and the sigmoid transformation
%       (4) difference between net value of the high effort option
%       and low effort option based on the model defined in .(choice/chosen).(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_mdl (='mdl_X' or 'bayesianModel_X')
%       This is basically like (3) but before the sigmoid transformation
%       (5) absolute difference between net value of the high effort option
%       and low effort option based on the model defined in .(choice/chosen).(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_mdl (='mdl_X' or 'bayesianModel_X')
%       This is basically like (2) but including the bias.
%
%       .Eperf.(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_varOption_bis
%       (1) difference between net value of the high effort option and low effort options based on the model defined in
%       .(choice/chosen).(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_mdl (='mdl_X' or 'bayesianModel_X')
%       (2) absolute difference between net value of the high effort option
%       and low effort option based on the model defined in .(choice/chosen).(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_mdl (='mdl_X' or 'bayesianModel_X')
%       (3) p(choice = hE) = probability of choosing the high effort option
%       based on the model defined in .(choice/chosen).(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_mdl (='mdl_X' or 'bayesianModel_X')
%       This is relatively similar to (2) but will include the motivational
%       bias and the sigmoid transformation
%       (4) difference between net value of the high effort option
%       and low effort option based on the model defined in .(choice/chosen).(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_mdl (='mdl_X' or 'bayesianModel_X')
%       This is basically like (3) but before the sigmoid transformation
%       (5) absolute difference between net value of the high effort option
%       and low effort option based on the model defined in .(choice/chosen).(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).NV_mdl (='mdl_X' or 'bayesianModel_X')
%       This is basically like (2) but including the bias.
%
%       .Eperf.(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).RT_1stAnswer
%       (1) raw reaction time for first answer (force above threshold for Ep 
%       and first answer to first digit for Em)
%
%       .Eperf.Ep.(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).fatigue:
%       physical effort only:
%       (1) sum of previous AUC of force produced
%       (2) same as (1) but zscored
%
%       .Eperf.Em.(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).prevEfficacy:
%       mental effort only:
%       (1) efficacy of previous trial computed as (nb correct answers - nb errors made)/total time of effort
%       period
%       (2) efficacy of previous trial computed  as (nb correct answers - nb errors made)/time spent between answer
%       of second useless number and end of the trial
%       (3) efficacy of previous trial computed as nb correct answers/total time of effort
%       period
%       (4) efficacy of previous trial computed  as nb correct answers/time spent between answer
%       of second useless number and end of the trial
%       (5-8) same as (1-4) but zscored
%
%       .Eperf.Ep.(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).Ech_x_fatigue:
%       physical effort only:
%       (1) (Effort chosen) x (sum of previous AUC of force produced)
%
%       .Eperf.Em.(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).Ech_x_prevEfficacy:
%       mental effort only:
%       (1) (Effort chosen) x (efficacy of previous trial computed as (nb correct answers - nb errors made)/total time of effort
%       period)
%       (2) (Effort chosen) x (efficacy of previous trial computed  as (nb correct answers - nb errors made)/time spent between answer
%       of second useless number and end of the trial)
%       (3) (Effort chosen) x (efficacy of previous trial computed as (nb correct
%       answers)/total time of effort period (no errors))
%       (4) (Effort chosen) x (efficacy of previous trial computed as (nb correct
%       answers)/total time of effort period (no errors) between answer
%       of second useless number and end of the trial)
%
%       .Eperf.(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).trialN
%       (1) trial number
%       (2) (trial number)*(E chosen - E non-chosen option)
%       (3) (trial number)*(E non-default - E default option)
%       (4) (trial number)*(E non-default)
%
%       .Eperf.(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).confidence
%       (1) confidence level (0/1) given by the subject for each choice
%       (2) confidence inferred by the model (p(high E)-0.5)² for the model
%       defined in .Eperf.(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).conf_mdl (='mdl_X' or 'bayesianModel_X')
%       (3) confidence inferred by the model (p(left)-0.5)² for the model
%       defined in .Eperf.(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).conf_mdl (='mdl_X' or 'bayesianModel_X')
%       (4) confidence inferred by the model (p(chosen)-0.5)² for the model
%       defined in .(choice/chosen).(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).conf_mdl (='mdl_X' or 'bayesianModel_X')
%
%% feedback period regressors
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
%       .fbk.(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).win_vs_loss
%       (1) win (1) - loss (0) trials
%
%       .fbk.(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).choiceHighE
%       (1) 0 when low effort/default option is selected and 1 when high
%       effort/non-default option is selected
%       (2) (-1) when low effort/default option is selected and (+1) when high
%       effort/non-default option is selected
%
%       .fbk.(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).money_obtained
%       (1) money amount obtained at the end of the trial
%       (1) |money amount| obtained at the end of the trial (~saliency)
%
%       .fbk.(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).E_made
%       (1) level of effort performed during effort performance
%
%       .fbk.(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).trialN
%       (1) trial number
%       (2) (trial number)*(E chosen - E non-chosen option)
%       (3) (trial number)*(E non-default - E default option)
%       (4) (trial number)*(E non-default)
%
%       .fbk.(Ep/Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).confidence
%       (1) confidence level (0/1) given by the subject for each choice
%       (2) confidence inferred by the model (p(high E)-0.5)² for the model
%       defined in .fbk.(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).conf_mdl (='mdl_X' or 'bayesianModel_X')
%       (3) confidence inferred by the model (p(left)-0.5)² for the model
%       defined in .fbk.(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).conf_mdl (='mdl_X' or 'bayesianModel_X')
%       (4) confidence inferred by the model (p(chosen)-0.5)² for the model
%       defined in .(choice/chosen).(Ep.Em).(R/P/RP).(E/E1/E2/E3/Ech0/Ech1/Ech2/Ech3/lEch/hEch).conf_mdl (='mdl_X' or 'bayesianModel_X')
%
% See also GLM_details.m
%
% Nicolas Clairis - 2021/2022


%% initialize all variables

% general parameters: derivative, grey matter filtering
[GLMprm.gal.add_drv,...
    GLMprm.gal.grey_mask,...
    GLMprm.gal.mask_probaThreshold,...
    GLMprm.gal.zPerRun,...
    GLMprm.gal.orth_vars,...
    GLMprm.gal.onsets_only,...
    GLMprm.gal.autocorrel] = deal(0);

% onsets: not modelled (none), modelled as stick function (stick) or as
% boxcar function (boxcar)
[GLMprm.model_onset.Ep.allCrosses,...
    GLMprm.model_onset.Ep.preChoiceCross,...
    GLMprm.model_onset.Ep.preEffortCross,...
    GLMprm.model_onset.Ep.choice,...
    GLMprm.model_onset.Ep.chosen,...
    GLMprm.model_onset.Ep.Eperf,...
    GLMprm.model_onset.Ep.fbk,...
    GLMprm.model_onset.Em.allCrosses,...
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
    GLMprm.preChoiceCross.(EpEm_nm).E_chosen = 0;
    
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
                GLMprm.choice.(EpEm_nm).(RP_nm).(Econd_nm).R_varOption,...
                GLMprm.choice.(EpEm_nm).(RP_nm).(Econd_nm).R_chosen,...
                GLMprm.choice.(EpEm_nm).(RP_nm).(Econd_nm).P_varOption,...
                GLMprm.choice.(EpEm_nm).(RP_nm).(Econd_nm).P_chosen,...
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
                GLMprm.choice.(EpEm_nm).(RP_nm).(Econd_nm).money_level_x_E_varOption,...
                GLMprm.choice.(EpEm_nm).(RP_nm).(Econd_nm).money_level_x_E_chosen,...
                GLMprm.choice.(EpEm_nm).(RP_nm).(Econd_nm).R_level_x_E_varOption,...
                GLMprm.choice.(EpEm_nm).(RP_nm).(Econd_nm).R_level_x_E_chosen,...
                GLMprm.choice.(EpEm_nm).(RP_nm).(Econd_nm).P_level_x_E_varOption,...
                GLMprm.choice.(EpEm_nm).(RP_nm).(Econd_nm).P_level_x_E_chosen,...
                GLMprm.choice.(EpEm_nm).(RP_nm).(Econd_nm).NV_chosen,...
                GLMprm.choice.(EpEm_nm).(RP_nm).(Econd_nm).NV_varOption,...
                GLMprm.choice.(EpEm_nm).(RP_nm).(Econd_nm).NV_varOption_bis,...
                GLMprm.choice.Ep.(RP_nm).(Econd_nm).F_integral,...
                GLMprm.choice.Em.(RP_nm).(Econd_nm).efficacy,...
                GLMprm.choice.Ep.(RP_nm).(Econd_nm).fatigue,...
                GLMprm.choice.Em.(RP_nm).(Econd_nm).prevEfficacy,...
                GLMprm.choice.Ep.(RP_nm).(Econd_nm).Ech_x_fatigue,...
                GLMprm.choice.Em.(RP_nm).(Econd_nm).Ech_x_prevEfficacy,...
                GLMprm.choice.(EpEm_nm).(RP_nm).(Econd_nm).confidence,...
                GLMprm.choice.(EpEm_nm).(RP_nm).(Econd_nm).RT,...
                GLMprm.choice.(EpEm_nm).(RP_nm).(Econd_nm).trialN] = deal(0);
            [GLMprm.choice.(EpEm_nm).(RP_nm).(Econd_nm).NV_mdl,...
                GLMprm.choice.(EpEm_nm).(RP_nm).(Econd_nm).conf_mdl] = deal('');
            
            % chosen option display
            % by default all regressors are not included
            [GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).R_vs_P,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).choiceHighE,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).R_varOption,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).R_chosen,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).P_varOption,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).P_chosen,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).money_left,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).money_right,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).money_chosen,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).money_unchosen,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).money_varOption,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).money_ch_min_unch,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).money_ch_min_fixOption,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).money_sum,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).E_left,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).E_right,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).E_chosen,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).E_unchosen,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).E_varOption,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).E_ch_min_unch,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).E_ch_min_fixOption,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).E_sum,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).money_level_x_E_varOption,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).money_level_x_E_chosen,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).R_level_x_E_varOption,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).R_level_x_E_chosen,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).P_level_x_E_varOption,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).P_level_x_E_chosen,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).NV_chosen,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).NV_varOption,...
                GLMprm.chosen.(EpEm_nm).(RP_nm).(Econd_nm).NV_varOption_bis,...
                GLMprm.chosen.Ep.(RP_nm).(Econd_nm).F_integral,...
                GLMprm.chosen.Em.(RP_nm).(Econd_nm).efficacy,...
                GLMprm.chosen.Ep.(RP_nm).(Econd_nm).fatigue,...
                GLMprm.chosen.Em.(RP_nm).(Econd_nm).prevEfficacy,...
                GLMprm.chosen.Ep.(RP_nm).(Econd_nm).Ech_x_fatigue,...
                GLMprm.chosen.Em.(RP_nm).(Econd_nm).Ech_x_prevEfficacy,...
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
                GLMprm.preEffortCross.(EpEm_nm).(RP_nm).(Econd_nm).NV_varOption_bis,...
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
                        GLMprm.preEffortCross.Em.(RP_nm).(Econd_nm).n_correct,...
                        GLMprm.preEffortCross.Em.(RP_nm).(Econd_nm).n_errors] = deal(0);
            end
            [GLMprm.preEffortCross.(EpEm_nm).(RP_nm).(Econd_nm).NV_mdl,...
                GLMprm.preEffortCross.(EpEm_nm).(RP_nm).(Econd_nm).conf_mdl] = deal('');
            
            % effort performance
            % by default all regressors are not included
            [GLMprm.Eperf.(EpEm_nm).(RP_nm).(Econd_nm).choiceHighE,...
                GLMprm.Eperf.(EpEm_nm).(RP_nm).(Econd_nm).money_chosen,...
                GLMprm.Eperf.(EpEm_nm).(RP_nm).(Econd_nm).E_chosen,...
                GLMprm.Eperf.(EpEm_nm).(RP_nm).(Econd_nm).R_varOption,...
                GLMprm.Eperf.(EpEm_nm).(RP_nm).(Econd_nm).P_varOption,...
                GLMprm.Eperf.(EpEm_nm).(RP_nm).(Econd_nm).R_chosen,...
                GLMprm.Eperf.(EpEm_nm).(RP_nm).(Econd_nm).P_chosen,...
                GLMprm.Eperf.(EpEm_nm).(RP_nm).(Econd_nm).R_level_x_E_varOption,...
                GLMprm.Eperf.(EpEm_nm).(RP_nm).(Econd_nm).P_level_x_E_varOption,...
                GLMprm.Eperf.(EpEm_nm).(RP_nm).(Econd_nm).R_level_x_E_chosen,...
                GLMprm.Eperf.(EpEm_nm).(RP_nm).(Econd_nm).P_level_x_E_chosen,...
                GLMprm.Eperf.(EpEm_nm).(RP_nm).(Econd_nm).NV_chosen,...
                GLMprm.Eperf.(EpEm_nm).(RP_nm).(Econd_nm).NV_varOption,...
                GLMprm.Eperf.(EpEm_nm).(RP_nm).(Econd_nm).NV_varOption_bis,...
                GLMprm.Eperf.(EpEm_nm).(RP_nm).(Econd_nm).RT_1stAnswer,...
                GLMprm.Eperf.(EpEm_nm).(RP_nm).(Econd_nm).trialN,...
                GLMprm.Eperf.(EpEm_nm).(RP_nm).(Econd_nm).confidence,...
                GLMprm.Eperf.Em.(RP_nm).(Econd_nm).RT_avg,...
                GLMprm.Eperf.Ep.(RP_nm).(Econd_nm).F_peak,...
                GLMprm.Eperf.Ep.(RP_nm).(Econd_nm).F_integral,...
                GLMprm.Eperf.Em.(RP_nm).(Econd_nm).efficacy,...
                GLMprm.Eperf.Em.(RP_nm).(Econd_nm).n_correct,...
                GLMprm.Eperf.Em.(RP_nm).(Econd_nm).n_errors,...
                GLMprm.Eperf.Ep.(RP_nm).(Econd_nm).fatigue,...
                GLMprm.Eperf.Em.(RP_nm).(Econd_nm).prevEfficacy,...
                GLMprm.Eperf.Ep.(RP_nm).(Econd_nm).Ech_x_fatigue,...
                GLMprm.Eperf.Em.(RP_nm).(Econd_nm).Ech_x_prevEfficacy] = deal(0);
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
Elvl_conds = {'E1','E2','E3'};
Ech_conds = {'Ech0','Ech1','Ech2','Ech3'};
Ech_conds_bis = {'lEch','hEch'};
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
        end % physical/mental loop
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
        end % physical/mental loop
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
        end % physical/mental loop
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
        end % physical/mental loop
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
        end % physical/mental loop
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
        end % physical/mental loop
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
        end % physical/mental loop
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
        end % physical/mental loop
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
        end % physical/mental loop
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
        end % physical/mental loop
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
        end % physical/mental loop
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
        end % physical/mental loop
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
        end % physical/mental loop
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
        end % physical/mental loop
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
        end % physical/mental loop
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
        end % physical/mental loop
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
        end % physical/mental loop
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
        end % physical/mental loop
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
        end % physical/mental loop
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
        end % physical/mental loop
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
        end % physical/mental loop
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
        end % physical/mental loop
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
        end % physical/mental loop
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
        end % physical/mental loop
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
        end % physical/mental loop
    case 67 % same as GLM 66 but with varying effort instead of chosen effort
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.money_varOption = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_varOption = 1;
            GLMprm.choice.(Epm_nm).RP.E.confidence = 2;
            GLMprm.choice.(Epm_nm).RP.E.conf_mdl = 'mdl_3';
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            GLMprm.Eperf.(Epm_nm).RP.E.E_chosen = 1;
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'boxcar';
        end % physical/mental loop
    case 68 % choice = high effort during choice
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.choiceHighE = 1;
            GLMprm.choice.(Epm_nm).RP.E.confidence = 2;
            GLMprm.choice.(Epm_nm).RP.E.conf_mdl = 'mdl_3';
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            GLMprm.Eperf.(Epm_nm).RP.E.E_chosen = 1;
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'boxcar';
        end % physical/mental loop
    case 69 % same as GLM 66 but orthogonalizaing variables
        % general parameters
        GLMprm.gal.orth_vars = 1;
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
        end % physical/mental loop
    case 70 % onsets-only GLM (like GLM64 but all events are modeled but with a stick, while in GLM 64 some have a boxcar)
        GLMprm.gal.onsets_only = 1;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % initial cross
            GLMprm.model_onset.(Epm_nm).preChoiceCross = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % pre-effort cross
            GLMprm.model_onset.(Epm_nm).preEffortCross = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 71 % same as GLM 66 but adding R vs P split of trials
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_vs_P = 1;
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
        end % physical/mental loop
    case 72 % GLM to test activation of VS during effort performance
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen option
            GLMprm.model_onset.(Epm_nm).chosen = 'boxcar';
            % effort performance (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            GLMprm.Eperf.(Epm_nm).RPpool = 0; % split R and P events
            for iRP = 1:length(RP_conds) % loop through R/P conditions
                RP_nm = RP_conds{iRP};
                GLMprm.Eperf.(Epm_nm).(RP_nm).E.money_chosen = 4;
                GLMprm.Eperf.Ep.(RP_nm).E.F_integral = 1;
                GLMprm.Eperf.Em.(RP_nm).E.efficacy = 2;
            end
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'boxcar';
        end % physical/mental loop
    case 73 % fatigue/learning impact on choice? + performance measure
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.Ep.RP.E.fatigue = 1;
            GLMprm.choice.Em.RP.E.prevEfficacy = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen option
            GLMprm.model_onset.(Epm_nm).chosen = 'boxcar';
            % effort performance (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            GLMprm.Eperf.(Epm_nm).RPpool = 0; % split R and P events
            for iRP = 1:length(RP_conds) % loop through R/P conditions
                RP_nm = RP_conds{iRP};
                GLMprm.Eperf.Ep.(RP_nm).E.F_integral = 1;
                GLMprm.Eperf.Em.(RP_nm).E.efficacy = 1;
            end
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'boxcar';
        end % physical/mental loop
    case 74 % look at fatigue during performance + net value high effort option during choice
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.NV_mdl = 'bayesianModel_3';
            GLMprm.choice.(Epm_nm).RP.E.NV_varOption = 1;
            % chosen option
            GLMprm.model_onset.(Epm_nm).chosen = 'boxcar';
            % effort performance (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'boxcar';
            GLMprm.Eperf.Ep.RP.E.F_integral = 1;
            GLMprm.Eperf.Em.RP.E.efficacy = 1;
            GLMprm.Eperf.Ep.RP.E.fatigue = 1;
            GLMprm.Eperf.Em.RP.E.prevEfficacy = 1;
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'boxcar';
        end % physical/mental loop
    case 75 % look for incentive effect during performance and for confidence
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.conf_mdl = 'bayesianModel_3';
            GLMprm.choice.(Epm_nm).RP.E.confidence = 2;
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen option
            GLMprm.model_onset.(Epm_nm).chosen = 'boxcar';
            % effort performance (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'boxcar';
            GLMprm.Eperf.(Epm_nm).RPpool = 0;
            for iRP = 1:length(RP_conds) % loop through R/P conditions
                RP_nm = RP_conds{iRP};
                GLMprm.Eperf.(Epm_nm).(RP_nm).E.money_chosen = 2;
                GLMprm.Eperf.Ep.(RP_nm).E.F_integral = 1;
                GLMprm.Eperf.Em.(RP_nm).E.efficacy = 1;
            end
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'boxcar';
        end % physical/mental loop
    case 76 % look at net value chosen option during choice
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.NV_mdl = 'bayesianModel_3';
            GLMprm.choice.(Epm_nm).RP.E.NV_chosen = 1;
            % chosen option
            GLMprm.model_onset.(Epm_nm).chosen = 'boxcar';
            % effort performance (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'boxcar';
            GLMprm.Eperf.Ep.RP.E.F_integral = 1;
            GLMprm.Eperf.Em.RP.E.efficacy = 1;
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'boxcar';
        end % physical/mental loop
    case 77 % model low vs high effort chosen separately during choice (similar to GLM59)
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).splitPerE = 3;
            GLMprm.choice.(Epm_nm).RP.lEch.money_varOption = 1;
            GLMprm.choice.(Epm_nm).RP.lEch.E_varOption = 1;
            GLMprm.choice.(Epm_nm).RP.hEch.money_varOption = 1;
            GLMprm.choice.(Epm_nm).RP.hEch.E_varOption = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'boxcar';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'boxcar';
        end % physical/mental loop
    case 78 % focus on R vs P during choice and during performance
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RPpool = 0;
            for iRP = 1:length(RP_conds)
                RP_nm = RP_conds{iRP};
                GLMprm.choice.(Epm_nm).(RP_nm).E.money_varOption = 2;
                GLMprm.choice.(Epm_nm).(RP_nm).E.E_varOption = 1;
            end
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'boxcar';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            GLMprm.Eperf.(Epm_nm).RPpool = 0;
            for iRP = 1:length(RP_conds)
                RP_nm = RP_conds{iRP};
                GLMprm.Eperf.(Epm_nm).(RP_nm).E.money_chosen = 2;
            end
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'boxcar';
        end % physical/mental loop
    case 79 % variant of GLM55/66, no zscore
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.money_ch_min_fixOption = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 80 % add derivative and pool chosen + cross together for effort preparation
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.add_drv = 1;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.money_ch_min_fixOption = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % effort preparation
            GLMprm.model_onset.(Epm_nm).chosen = 'boxcar_ter';
            GLMprm.chosen.(Epm_nm).RP.E.E_chosen = 1;
            % effort execution
            GLMprm.model_onset.(Epm_nm).Eperf = 'boxcar';
            GLMprm.Eperf.(Epm_nm).RP.E.E_chosen = 1;
            % feedback
            GLMprm.model_onset.(Epm_nm).Eperf = 'boxcar';
        end % physical/mental loop
    case 81 % similar to GLM 80 but without modeling effort during effort preparation
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.add_drv = 1;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.money_ch_min_fixOption = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % effort preparation
            GLMprm.model_onset.(Epm_nm).chosen = 'boxcar';
            % effort execution
            GLMprm.model_onset.(Epm_nm).Eperf = 'boxcar';
            GLMprm.Eperf.(Epm_nm).RP.E.E_chosen = 1;
            % feedback
            GLMprm.model_onset.(Epm_nm).Eperf = 'boxcar';
        end % physical/mental loop
    case 82 % like GLM 79 but adding temporal derivative
        % general parameters
        GLMprm.gal.add_drv = 1;
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.money_ch_min_fixOption = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 83 % look at correlates of |delta net value| for identifying hard trials
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.NV_varOption = 2;
            GLMprm.choice.(Epm_nm).RP.E.NV_mdl = 'bayesianModel_3';
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'boxcar';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'boxcar';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'boxcar';
        end % physical/mental loop
    case 84 % same as GLM 83, including Effort
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.NV_varOption = 2;
            GLMprm.choice.(Epm_nm).RP.E.NV_mdl = 'bayesianModel_3';
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'boxcar';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'boxcar';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'boxcar';
        end % physical/mental loop
    case 85 % focus on the slope for high efforts chosen
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            % split low/high effort choice and look at the slope for the
            % high effort choice
            GLMprm.choice.(Epm_nm).splitPerE = 3;
            GLMprm.choice.(Epm_nm).RP.hEch.money_ch_min_fixOption = 1;
            GLMprm.choice.(Epm_nm).RP.hEch.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.hEch.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 86 % GLM only with Effort chosen during choice
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 87 % GLM similar to 86 but including several confounding factors (NV, RT, confidence)
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.NV_varOption = 3;
            GLMprm.choice.(Epm_nm).RP.E.NV_mdl = 'bayesianModel_3';
            GLMprm.choice.(Epm_nm).RP.E.confidence = 2;
            GLMprm.choice.(Epm_nm).RP.E.conf_mdl = 'bayesianModel_3';
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 88 % GLM with Effort, Reward and Punishment chosen during choice
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 89 % same as GLM 88 but including several confounding factors (NV, RT, confidence)
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.NV_varOption = 3;
            GLMprm.choice.(Epm_nm).RP.E.NV_mdl = 'bayesianModel_3';
            GLMprm.choice.(Epm_nm).RP.E.confidence = 2;
            GLMprm.choice.(Epm_nm).RP.E.conf_mdl = 'bayesianModel_3';
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 90 % onsets-only GLM identical to GLM88 for the periods modeled but with no parametric modulation
        % general parameters
        GLMprm.gal.onsets_only = 1;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 91 % same as GLM 88 but with boxcar instead of stick
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'boxcar';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'boxcar';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'boxcar';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'boxcar';
        end % physical/mental loop
    case 92 % same as GLM 89 but with boxcar instead of stick
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'boxcar';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.NV_varOption = 3;
            GLMprm.choice.(Epm_nm).RP.E.NV_mdl = 'bayesianModel_3';
            GLMprm.choice.(Epm_nm).RP.E.confidence = 2;
            GLMprm.choice.(Epm_nm).RP.E.conf_mdl = 'bayesianModel_3';
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'boxcar';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'boxcar';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'boxcar';
        end % physical/mental loop
    case 93 % same as GLM 90 but with boxcar instead of stick
        % general parameters
        GLMprm.gal.onsets_only = 1;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'boxcar';
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'boxcar';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'boxcar';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'boxcar';
        end % physical/mental loop
    case 94 % same as GLM 90 but no chosen period
        % general parameters
        GLMprm.gal.onsets_only = 1;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 95 % same as GLM 93 but no chosen period
        % general parameters
        GLMprm.gal.onsets_only = 1;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'boxcar';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'boxcar';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'boxcar';
        end % physical/mental loop
    case 96 % like GLM 88 but no chosen period
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 97 % like GLM 89 but no chosen period
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.NV_varOption = 3;
            GLMprm.choice.(Epm_nm).RP.E.NV_mdl = 'bayesianModel_3';
            GLMprm.choice.(Epm_nm).RP.E.confidence = 2;
            GLMprm.choice.(Epm_nm).RP.E.conf_mdl = 'bayesianModel_3';
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 98 % test whether BOLD baseline will predict choice made
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % pre-choice fixation cross
            GLMprm.model_onset.(Epm_nm).preChoiceCross = 'boxcar_bis';
            GLMprm.preChoiceCross.(Epm_nm).choiceHighE = 1;
            % chosen option display
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 99 % test whether BOLD baseline will predict Effort chosen
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % pre-choice fixation cross
            GLMprm.model_onset.(Epm_nm).preChoiceCross = 'boxcar_bis';
            GLMprm.preChoiceCross.(Epm_nm).E_chosen = 1;
            % chosen option display
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 100 % like GLM 96 but boxcar
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'boxcar';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'boxcar';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'boxcar';
        end % physical/mental loop
    case 101 % like GLM 97 but boxcar
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'boxcar';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.NV_varOption = 3;
            GLMprm.choice.(Epm_nm).RP.E.NV_mdl = 'bayesianModel_3';
            GLMprm.choice.(Epm_nm).RP.E.confidence = 2;
            GLMprm.choice.(Epm_nm).RP.E.conf_mdl = 'bayesianModel_3';
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'boxcar';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'boxcar';
        end % physical/mental loop
    case 102 % extract average dmPFC per effort level proposed
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).splitPerE = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 103 % extract average dmPFC per period
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 104 % like GLM 96 but adding time component to the GLM (Ep fatigue and Em prev. efficacy)
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 105 % like GLM 97 but adding time component to the GLM (Ep fatigue and Em prev. efficacy)
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            GLMprm.choice.(Epm_nm).RP.E.NV_varOption = 3;
            GLMprm.choice.(Epm_nm).RP.E.NV_mdl = 'bayesianModel_3';
            GLMprm.choice.(Epm_nm).RP.E.confidence = 2;
            GLMprm.choice.(Epm_nm).RP.E.conf_mdl = 'bayesianModel_3';
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 106 % like GLM 104 but Ech*fatigue instead of fatigue
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.Ech_x_fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.Ech_x_prevEfficacy = 3;
            end
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 107 % like GLM 105 but Ech*fatigue instead of fatigue
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.Ech_x_fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.Ech_x_prevEfficacy = 3;
            end
            GLMprm.choice.(Epm_nm).RP.E.NV_varOption = 3;
            GLMprm.choice.(Epm_nm).RP.E.NV_mdl = 'bayesianModel_3';
            GLMprm.choice.(Epm_nm).RP.E.confidence = 2;
            GLMprm.choice.(Epm_nm).RP.E.conf_mdl = 'bayesianModel_3';
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 108 % look at slope with effort, but only for when high effort option is chosen
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).splitPerE = 3; % split trials according to choice made
            GLMprm.choice.(Epm_nm).RP.hEch.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.hEch.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.hEch.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.hEch.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.hEch.prevEfficacy = 3;
            end
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 109 % look at SV correlates during performance
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            GLMprm.Eperf.(Epm_nm).RP.E.NV_chosen = 2;
            GLMprm.Eperf.(Epm_nm).RP.E.NV_mdl = 'bayesianModel_3';
            switch Epm_nm
                case 'Ep'
                    GLMprm.Eperf.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.Eperf.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 110 % similar to GLM109 but only E chosen and fatigue modulator during performance
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            GLMprm.Eperf.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.Eperf.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.Eperf.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 111 % looking at efficacy (in Em), force integral (in Ep) and latency
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            switch Epm_nm
                case 'Ep'
                    GLMprm.Eperf.Ep.RP.E.F_integral = 3;
                case 'Em'
                    GLMprm.Eperf.Em.RP.E.efficacy = 3;
            end
            GLMprm.Eperf.(Epm_nm).RP.E.RT_1stAnswer = 1;
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 112 % looking at nb correct (in Em) and force integral (in Ep)
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            switch Epm_nm
                case 'Ep'
                    GLMprm.Eperf.Ep.RP.E.F_integral = 3;
                case 'Em'
                    GLMprm.Eperf.Em.RP.E.n_correct = 1;
            end
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 113 % look at performance like GLM 111 but separately for each effort level
        error('need debugging');
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            GLMprm.Eperf.(Epm_nm).splitPerE = 2;
            for iEch = 1:length(Ech_conds)
                Ech_nm = Ech_conds{iEch};
                switch Epm_nm
                    case 'Ep'
                        GLMprm.Eperf.Ep.RP.(Ech_nm).F_integral = 3;
                    case 'Em'
                        GLMprm.Eperf.Em.RP.(Ech_nm).efficacy = 3;
                end
                GLMprm.Eperf.(Epm_nm).RP.(Ech_nm).RT_1stAnswer = 1;
            end
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 114 % look at performance like GLM 112 but separately for each effort level
        error('need debugging');
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            switch Epm_nm
                case 'Ep'
                    GLMprm.Eperf.Ep.RP.E.F_integral = 3;
                case 'Em'
                    GLMprm.Eperf.Em.RP.E.n_correct = 1;
            end
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 115 % like GLM 104 but focus on high effort option instead of chosen option
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_varOption = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_varOption = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_varOption = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 116 % like GLM 104 and 115 but include both high effort option and chosen effort and see which is best
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_varOption = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_varOption = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_varOption = 1;
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 117 % like GLM 104 but include both high effort option and chosen effort and see which is best
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.E_varOption = 1;
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 118 % like GLM 104 but adding RT and confidence in the model
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            GLMprm.choice.(Epm_nm).RP.E.confidence = 2;
            GLMprm.choice.(Epm_nm).RP.E.conf_mdl = 'bayesianModel_3';
            GLMprm.choice.(Epm_nm).RP.E.RT = 2;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 119 % like GLM 104 but adding RT, p(choice) and confidence in the model
        % +compared to GLM 118, the confidence is defined slightly
        % differently
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            GLMprm.choice.(Epm_nm).RP.E.NV_varOption = 3;
            GLMprm.choice.(Epm_nm).RP.E.NV_mdl = 'bayesianModel_3';
            GLMprm.choice.(Epm_nm).RP.E.confidence = 3;
            GLMprm.choice.(Epm_nm).RP.E.conf_mdl = 'bayesianModel_3';
            GLMprm.choice.(Epm_nm).RP.E.RT = 2;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 120 % like GLM 104 but with grey-matter filter instead of SPM's default
        % with the objective to get back some signal in the striatum
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.grey_mask = 3;
        GLMprm.gal.zPerRun = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 121
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.grey_mask = 3;
        GLMprm.gal.zPerRun = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_varOption = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_varOption = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_varOption = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
        %% GLM to check % signal change dmPFC during choice vs % high E
    case 122 % like GLM 104 without any parametric regulator
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        %%
        
    case 123 % like GLM 104 but more lenient SPM threshold at 1st level to include striatum
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 4;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 124 % like GLM 104 and 123 but even more lenient SPM threshold at 1st level to include striatum
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 125 % like GLM 104, 123 and 124 but even more lenient SPM threshold at 1st level to include striatum
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 6;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 126 % like GLM 104, 123, 124 and 125 but using grey + white matter mask
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 127 % like GLM 104, 123, 124, 125 and 126 but using grey + white matter mask and more lenient threshold for mask
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 128 % same as GLM 94 but with grey+white matter mask
        % general parameters
        GLMprm.gal.onsets_only = 1;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 129 % same as GLM 94 but with SPM implicit mask more lenient at 0.1
        % general parameters
        GLMprm.gal.onsets_only = 1;
        GLMprm.gal.grey_mask = 6;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 130 % like GLM 127 but only regressor = choice
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.choiceHighE = 1;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 131 % like GLM 130 but p(high E) instead of actual choice
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.NV_varOption = 3;
            GLMprm.choice.(Epm_nm).RP.E.NV_mdl = 'bayesianModel_3';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 132 % split onsets per effort level and modulate by incentive
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).splitPerE = 1;
            for iE = 1:length(Elvl_conds)
                E_lvl_nm = Elvl_conds{iE};
                GLMprm.choice.(Epm_nm).RP.(E_lvl_nm).R_varOption = 2;
                GLMprm.choice.(Epm_nm).RP.(E_lvl_nm).P_varOption = 2;
            end % effort loop
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 133 % check interaction R*E and P*E
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_level_x_E_varOption = 1;
            GLMprm.choice.(Epm_nm).RP.E.P_level_x_E_varOption = 1;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 134 % like GLM 127 but adding R>P and RT as regressors
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_vs_P = 1;
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            GLMprm.choice.(Epm_nm).RP.E.RT = 2;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 135 % like GLM 127 but adding Conf and z(RT) as regressors
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            GLMprm.choice.(Epm_nm).RP.E.confidence = 3;
            GLMprm.choice.(Epm_nm).RP.E.conf_mdl = 'bayesianModel_3';
            GLMprm.choice.(Epm_nm).RP.E.RT = 2;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 136 % like GLM 127 but adding z(RT) as regressors
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            GLMprm.choice.(Epm_nm).RP.E.RT = 2;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 137 % like GLM 127 but adding Conf and raw RT as regressors
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            GLMprm.choice.(Epm_nm).RP.E.confidence = 3;
            GLMprm.choice.(Epm_nm).RP.E.conf_mdl = 'bayesianModel_3';
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 138 % like GLM 127 but adding raw RT as regressor
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 0;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 139 % like GLM 118 but with grey+white matter mask at 20%
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 20;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            GLMprm.choice.(Epm_nm).RP.E.confidence = 3;
            GLMprm.choice.(Epm_nm).RP.E.conf_mdl = 'bayesianModel_3';
            GLMprm.choice.(Epm_nm).RP.E.RT = 2;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 140 % like GLM 139 but with variables orthogonalized to RT
        % general parameters
        GLMprm.gal.orth_vars = 1;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 20;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.RT = 5;
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            GLMprm.choice.(Epm_nm).RP.E.confidence = 3;
            GLMprm.choice.(Epm_nm).RP.E.conf_mdl = 'bayesianModel_3';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
        %% new series of GLMS with temporal derivative + mask at 20%
    case 141 % like GLM 139 but with temporal derivative
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.add_drv = 1;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 20;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            GLMprm.choice.(Epm_nm).RP.E.confidence = 3;
            GLMprm.choice.(Epm_nm).RP.E.conf_mdl = 'bayesianModel_3';
            GLMprm.choice.(Epm_nm).RP.E.RT = 2;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 142 % same as GLM 127 but threshold at 20% for mask + adding temporal derivative
        % general parameters
        GLMprm.gal.add_drv = 1;
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 20;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 143 % same as GLM 142 but adding z(RT)
        % general parameters
        GLMprm.gal.add_drv = 1;
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 20;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            GLMprm.choice.(Epm_nm).RT.E.RT = 2;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
        %% new series of GLM with mask at 20% but no temporal derivative
    case 144 % same as GLM 127 but threshold at 20% for mask
        % and same as GLM 142 but without temporal derivative
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 20;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 145 % same as GLM 144 but only keep Ech regressor
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 20;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
        %% 20% mask + temporal derivative
    case 146 % same as GLM 142 but E chosen as the only regressor
        % general parameters
        GLMprm.gal.add_drv = 1;
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 20;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 147 % same as GLM 142 but no time-varying regressor (physical fatigue + mental facilitation)
        % general parameters
        GLMprm.gal.add_drv = 1;
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 20;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
        %% new series of GLM back to more lenient mask to avoid problems close to midline
    case 148 % (Kurniawan et al, 2021)-like GLM with 5% mask
        % R/P/E for high E option + Ech + RT + choice
        % general parameters
        GLMprm.gal.add_drv = 1;
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.choiceHighE = 1;
            GLMprm.choice.(Epm_nm).RP.E.R_varOption = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_varOption = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_varOption = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.RT = 2;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 149 % (Kurniawan et al, 2021)-like GLM without temporal derivative
        % R/P/E for high E option + Ech + RT + choice
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.choiceHighE = 1;
            GLMprm.choice.(Epm_nm).RP.E.R_varOption = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_varOption = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_varOption = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.RT = 2;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
        
    case 150 % basic GLM with Rch/Pch/Ech
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 151 % GLM 150 + Fp/Fm
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 152 % GLM 151 + z(RT)
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            GLMprm.choice.(Epm_nm).RP.E.RT = 2;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 153 % GLM 152 + model-derived confidence
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            GLMprm.choice.(Epm_nm).RP.E.RT = 2;
            GLMprm.choice.(Epm_nm).RP.E.confidence = 2;
            GLMprm.choice.(Epm_nm).RP.E.conf_mdl = 'bayesianModel_3';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 154 % p(high E) as only regressor
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.NV_varOption = 3;
            GLMprm.choice.(Epm_nm).RP.E.NV_mdl = 'bayesianModel_3';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 155 % dR/dP/dE/Ech
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_varOption = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_varOption = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_varOption = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 156 % like GLM 150 but split R and P trials
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RPpool = 0;
            % R  trials
            GLMprm.choice.(Epm_nm).R.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).R.E.E_chosen = 1;
            % P trials
            GLMprm.choice.(Epm_nm).P.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).P.E.E_chosen = 1;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 157 % like GLM 156 but adding chosen period
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RPpool = 0;
            % R  trials
            GLMprm.choice.(Epm_nm).R.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).R.E.E_chosen = 1;
            % P trials
            GLMprm.choice.(Epm_nm).P.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).P.E.E_chosen = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        %% onsets-only
    case 158 % same as GLM 128 but threshold at 5%
        % general parameters
        GLMprm.gal.onsets_only = 1;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        %%
        
    case 159 % like GLM 150 but adding chosen period
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 160 % p(chosen)
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.NV_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.NV_mdl = 'bayesianModel_3';
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 161 % competition between -p(chosen) and Echosen
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.NV_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.NV_mdl = 'bayesianModel_3';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 162 % competition between -p(chosen), Echosen and RT
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.NV_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.NV_mdl = 'bayesianModel_3';
            GLMprm.choice.(Epm_nm).RP.E.RT = 2;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        %% GLM 163-164 testing interaction
    case 163 % R*E chosen and P*E chosen
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_level_x_E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.P_level_x_E_chosen = 1;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 164 % R*E chosen and P*E chosen + control regressors (Rch/Pch/Ech alone and RT)
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.R_level_x_E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.P_level_x_E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.RT = 2;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 165 % dR/dP/dE; like GLM 155 without E chosen
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_varOption = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_varOption = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_varOption = 1;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 166 % Rch/Pch/Ech/Fp-Fm/Conf/RT with zscore
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            GLMprm.choice.(Epm_nm).RP.E.confidence = 2;
            GLMprm.choice.(Epm_nm).RP.E.conf_mdl = 'bayesianModel_3';
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen option display
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 167 % same as GLM 166 but with SPM implicit mask
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            GLMprm.choice.(Epm_nm).RP.E.confidence = 2;
            GLMprm.choice.(Epm_nm).RP.E.conf_mdl = 'bayesianModel_3';
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen option display
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 168 % look at effort invested but during incentive/RT
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.F_integral = 3;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.efficacy = 1;
            end
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen option display
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 169 % like GLM168 but adding Echosen
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.F_integral = 3;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.efficacy = 1;
            end
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen option display
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 170 % like GLM169 but orthogonalizing variables + RT at the beginning
        % general parameters
        GLMprm.gal.orth_vars = 1;
        GLMprm.gal.zPerRun = 1;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.F_integral = 3;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.efficacy = 1;
            end
            GLMprm.choice.(Epm_nm).RP.E.RT = 4;
            % chosen option display
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 171 % like GLM 166 but no zscore (Rch/Pch/Ech/Fp-Fm/Conf/RT)
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            GLMprm.choice.(Epm_nm).RP.E.confidence = 2;
            GLMprm.choice.(Epm_nm).RP.E.conf_mdl = 'bayesianModel_3';
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen option display
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        %% new series of GLM
    case 172 % Rch/Pch/Ech/Fp/Fm during choice
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            % chosen period
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 173 % Rch/Pch/Ech/Fp/Fm during chosen
        % like GLM 172 but modulating chosen instead of choice period
        
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            % chosen period
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            GLMprm.chosen.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.chosen.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.chosen.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.chosen.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.chosen.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 174 % Rch/Pch/Ech/Fp/Fm during choice
        % like GLM 172 but with boxcar everywhere
        
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'boxcar_bis';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'boxcar';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'boxcar';
        end % physical/mental loop
    case 175 % Rch/Pch/Ech/Fp/Fm during choice
        % like GLM 174 but with boxcar only for choice+chosen
        
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'boxcar_bis';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
        %% new series of GLM
    case 176 % Rch/Pch/Ech during choice
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            % chosen period
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 177 % Rch/Pch/Ech during chosen
        % like GLM 176 but modulating chosen instead of choice period
        
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            % chosen period
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            GLMprm.chosen.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.chosen.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.chosen.(Epm_nm).RP.E.E_chosen = 1;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 178 % Rch/Pch/Ech during choice+chosen
        
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'boxcar_bis';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
        %% new series of GLM with boxcars
    case 179 % Rch/Pch/Ech with boxcars + zscore variables + performance during Eperf
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'boxcar';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'boxcar';
            switch Epm_nm
                case 'Ep'
                    GLMprm.Eperf.(Epm_nm).RP.E.F_integral = 3;
                    GLMprm.Eperf.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.Eperf.(Epm_nm).RP.E.efficacy = 1;
                    GLMprm.Eperf.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'boxcar';
        end % physical/mental loop
        
    case 180 % Rch/Pch/Ech with boxcars + zscore variables
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'boxcar';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'boxcar';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'boxcar';
        end % physical/mental loop
    case 181 % Rch/Pch/Ech with stick and boxcars
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'boxcar';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 182 % Rch/Pch/Ech with stick and boxcars + FAST
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        GLMprm.gal.autocorrel = 1;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'boxcar';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop

    case 183 % like GLM 151 but with chosen period modelled
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            % chosen
            GLMprm.chosen.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 184 % variant of GLM 168 including R/Pchosen and choice
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.choiceHighE = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.F_integral = 3;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.efficacy = 1;
            end
            % chosen option display
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 185 % variant of GLM 168 including R/Pchosen and choice
        % same as GLM 184 but without chosen period
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.choiceHighE = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.F_integral = 3;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.efficacy = 1;
            end
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 186 % same as GLM 158 but adding chosen period
        % general parameters
        GLMprm.gal.onsets_only = 1;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    
        %% new series of GLMs: back to pooling R and P together
    case 187 % money chosen/E chosen + zscore
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.money_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 188 % money chosen/E chosen/Fp-Fm + zscore
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.money_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 189 % money chosen/E chosen/Fp-Fm/RT + zscore
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.money_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 190 % money chosen/E chosen/RT + zscore
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.money_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 191 % basic GLM with Rch/Pch/Ech (like GLM 150 but with orthogonalization of variables)
        % general parameters
        GLMprm.gal.orth_vars = 1;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 192 % GLM trying to model all the periods of the trial
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % initial cross
            GLMprm.model_onset.(Epm_nm).preChoiceCross = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.money_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % pre-effort cross
            GLMprm.model_onset.(Epm_nm).preEffortCross = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 193 % basic GLM with Money ch/Ech/RT
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.money_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 194 % GLM trying to model all the periods of the trial
        % same as GLM 192 but split money into R and P regressors
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % initial cross
            GLMprm.model_onset.(Epm_nm).preChoiceCross = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % pre-effort cross
            GLMprm.model_onset.(Epm_nm).preEffortCross = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 195 % GLM trying to model all the periods of the trial
        % similar to GLM 192 but using boxcars instead of sticks
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % initial cross
            GLMprm.model_onset.(Epm_nm).preChoiceCross = 'boxcar';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'boxcar';
            GLMprm.choice.(Epm_nm).RP.E.money_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'boxcar';
            % pre-effort cross
            GLMprm.model_onset.(Epm_nm).preEffortCross = 'boxcar';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'boxcar';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'boxcar';
        end % physical/mental loop
    case 196 % similar to GLM 195 but adding time regressor
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % initial cross
            GLMprm.model_onset.(Epm_nm).preChoiceCross = 'boxcar';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'boxcar';
            GLMprm.choice.(Epm_nm).RP.E.money_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'boxcar';
            % pre-effort cross
            GLMprm.model_onset.(Epm_nm).preEffortCross = 'boxcar';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'boxcar';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'boxcar';
        end % physical/mental loop
    case 197 % basic GLM with Rch/Pch/Ech (like GLM150) but modeling all trial periods + using boxcars
        % + one single regressor for all fixation crosses
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 7;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % fixation crosses
            GLMprm.model_onset.(Epm_nm).allCrosses = 'boxcar';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'boxcar';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            % chosen
            GLMprm.chosen.(Epm_nm).RP.E.chosen = 'boxcar';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'boxcar';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'boxcar';
        end % physical/mental loop
        
        %% from here and on using grey matter mask (5%) instead of grey + white matter
        %% GLM with stick + all trial periods
    case 198 % Rch/Pch/Ech/RT
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        GLMprm.gal.grey_mask = 3;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % fixation crosses
            GLMprm.model_onset.(Epm_nm).allCrosses = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 199 % same as GLM 198 but adding Fp-Fm
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        GLMprm.gal.grey_mask = 3;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % fixation crosses
            GLMprm.model_onset.(Epm_nm).allCrosses = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
        %%
    case 200 % same as GLM 192 but removing RT + using grey-matter mask instead of grey+white matter
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        GLMprm.gal.grey_mask = 3;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % initial cross
            GLMprm.model_onset.(Epm_nm).preChoiceCross = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.money_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % pre-effort cross
            GLMprm.model_onset.(Epm_nm).preEffortCross = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 201 % same as GLM 186 but grey matter mask only (not grey+white)
        % general parameters
        GLMprm.gal.onsets_only = 1;
        GLMprm.gal.grey_mask = 3;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 202 % Rch/Pch/Ech/RT
        % variant on GLM 198 but including chosen period + using R-P levels
        % instead of amounts + using zscore only for RT
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 3;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % fixation crosses
            GLMprm.model_onset.(Epm_nm).allCrosses = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.RT = 2;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 203 % same as GLM 198 but adding back 'chosen' period removed by mistake
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        GLMprm.gal.grey_mask = 3;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % fixation crosses
            GLMprm.model_onset.(Epm_nm).allCrosses = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 204 % same as GLM 200 but adding back RT
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        GLMprm.gal.grey_mask = 3;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % initial cross
            GLMprm.model_onset.(Epm_nm).preChoiceCross = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.money_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % pre-effort cross
            GLMprm.model_onset.(Epm_nm).preEffortCross = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 205 % same as GLM 203 but without RT
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        GLMprm.gal.grey_mask = 3;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % fixation crosses
            GLMprm.model_onset.(Epm_nm).allCrosses = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 206 % Rch/Pch/Ech all periods (but cross) modeled + boxcar everywhere
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.grey_mask = 3;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'boxcar';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'boxcar';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'boxcar';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'boxcar';
        end % physical/mental loop
        
    case 207 % like GLM 206 but stick instead of boxcars
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.grey_mask = 3;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 208 % like GLM 150 but with grey matter
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 3;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 209 % like GLM 208 but adding RT
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 3;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.RT = 2;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 210 % GLM including all possible confounding variables
        % Val(dSV+bias)/Saliency(|dSV+bias|)/Ech/Fp-Fm/RT/Uncertainty
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        GLMprm.gal.grey_mask = 3;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % fixation crosses
            GLMprm.model_onset.(Epm_nm).allCrosses = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.NV_mdl = 'bayesianModel_3';
            GLMprm.choice.(Epm_nm).RP.E.NV_varOption = 4;
            GLMprm.choice.(Epm_nm).RP.E.NV_varOption_bis = 5;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.confidence = 3;
            GLMprm.choice.(Epm_nm).RP.E.conf_mdl = 'bayesianModel_3';
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 211 % GLM including all possible confounding variables
        % Val(dSV+bias)/Saliency(|dSV+bias|)/Ech/Fp-Fm/RT/Uncertainty
        % same as GLM 210 but using p(chosen) for Confidence instead of
        % p(left)
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        GLMprm.gal.grey_mask = 3;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % fixation crosses
            GLMprm.model_onset.(Epm_nm).allCrosses = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.NV_mdl = 'bayesianModel_3';
            GLMprm.choice.(Epm_nm).RP.E.NV_varOption = 4;
            GLMprm.choice.(Epm_nm).RP.E.NV_varOption_bis = 5;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.confidence = 4;
            GLMprm.choice.(Epm_nm).RP.E.conf_mdl = 'bayesianModel_3';
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 212 % reduced version of GLM 211: Val/Conf/RT/Fp-Fm/Ech (removed saliency)
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        GLMprm.gal.grey_mask = 3;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % fixation crosses
            GLMprm.model_onset.(Epm_nm).allCrosses = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.NV_mdl = 'bayesianModel_3';
            GLMprm.choice.(Epm_nm).RP.E.NV_varOption = 4;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.confidence = 4;
            GLMprm.choice.(Epm_nm).RP.E.conf_mdl = 'bayesianModel_3';
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 213 % variant of GLM 210-212: Val(chosen)/Saliency/RT/Fp-Fm/Ech (no confidence)
        % bug for some subjects due to high correlation between Val(chosen)
        % and saliency
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        GLMprm.gal.grey_mask = 3;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % fixation crosses
            GLMprm.model_onset.(Epm_nm).allCrosses = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.NV_mdl = 'bayesianModel_3';
            GLMprm.choice.(Epm_nm).RP.E.NV_chosen = 3;
            GLMprm.choice.(Epm_nm).RP.E.NV_varOption = 5;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 214 % same as GLM 213 but without saliency
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        GLMprm.gal.grey_mask = 3;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % fixation crosses
            GLMprm.model_onset.(Epm_nm).allCrosses = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.NV_mdl = 'bayesianModel_3';
            GLMprm.choice.(Epm_nm).RP.E.NV_chosen = 3;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 215 % same as GLM 203 but adding time effect (Fp/Fm)
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        GLMprm.gal.grey_mask = 3;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % fixation crosses
            GLMprm.model_onset.(Epm_nm).allCrosses = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 216 % same as GLM 213 but without value
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        GLMprm.gal.grey_mask = 3;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % fixation crosses
            GLMprm.model_onset.(Epm_nm).allCrosses = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.NV_mdl = 'bayesianModel_3';
            GLMprm.choice.(Epm_nm).RP.E.NV_varOption = 5;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 217 % same as GLM 215 but using R/P levels instead of amounts
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        GLMprm.gal.grey_mask = 3;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % fixation crosses
            GLMprm.model_onset.(Epm_nm).allCrosses = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 218 % same as GLM 212 but without NVhE-NVlE
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        GLMprm.gal.grey_mask = 3;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % fixation crosses
            GLMprm.model_onset.(Epm_nm).allCrosses = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.confidence = 4;
            GLMprm.choice.(Epm_nm).RP.E.conf_mdl = 'bayesianModel_3';
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 219 % same as GLM 217 but without Fp/Fm (ie like GLM 203 but with Rch/Pch levels instead of amounts)
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        GLMprm.gal.grey_mask = 3;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % fixation crosses
            GLMprm.model_onset.(Epm_nm).allCrosses = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 220  % same as GLM 217 but adding performance regressor during performance
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        GLMprm.gal.grey_mask = 3;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % fixation crosses
            GLMprm.model_onset.(Epm_nm).allCrosses = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            switch Epm_nm
                case 'Ep'
                    GLMprm.Eperf.Ep.RP.E.F_integral = 3;
                case 'Em'
                    GLMprm.Eperf.Em.RP.E.efficacy = 3;
            end
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
        %% new family of GLM: fixed R/P regressors with zscore
    case 221 % Rch/Pch/Ech/Fp-Fm/RT
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 3;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % fixation crosses
            GLMprm.model_onset.(Epm_nm).allCrosses = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 4;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 4;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.Ep.RP.E.fatigue = 2;
                case 'Em'
                    GLMprm.choice.Em.RP.E.prevEfficacy = 7;
            end
            GLMprm.choice.(Epm_nm).RP.E.RT = 2;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 222 % same as GLM 221 but adding perf during perf
        % Rch/Pch/Ech/Fp-Fm/RT during choice + perf during perf
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 3;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % fixation crosses
            GLMprm.model_onset.(Epm_nm).allCrosses = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 4;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 4;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.Ep.RP.E.fatigue = 2;
                case 'Em'
                    GLMprm.choice.Em.RP.E.prevEfficacy = 7;
            end
            GLMprm.choice.(Epm_nm).RP.E.RT = 2;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            switch Epm_nm
                case 'Ep'
                    GLMprm.Eperf.Ep.RP.E.F_integral = 7;
                case 'Em'
                    GLMprm.Eperf.Em.RP.E.efficacy = 7;
            end
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 223 % Rch/Pch/Ech/RT using levels and no zscore apart for RT
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 3;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % fixation crosses
            GLMprm.model_onset.(Epm_nm).allCrosses = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.RT = 2;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
    case 224 % Rch/Pch/Ech/Fp/Fm/RT using levels and zscore only Fp/Fm/RT variables
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 3;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % fixation crosses
            GLMprm.model_onset.(Epm_nm).allCrosses = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.Ep.RP.E.fatigue = 2;
                case 'Em'
                    GLMprm.choice.Em.RP.E.prevEfficacy = 7;
            end
            GLMprm.choice.(Epm_nm).RP.E.RT = 2;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 225 % same as GLM 223 but with boxcar + removing cross regressor
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 3;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'boxcar';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.RT = 2;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'boxcar';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'boxcar';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'boxcar';
        end % physical/mental loop
    case 226 % same as GLM 223 but removing cross regressor
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 3;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.R_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_chosen = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.RT = 2;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 227 % same as GLM 214 but without cross
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        GLMprm.gal.grey_mask = 3;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.NV_mdl = 'bayesianModel_3';
            GLMprm.choice.(Epm_nm).RP.E.NV_chosen = 3;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 228 % same as GLM 227 but without zscore
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        GLMprm.gal.grey_mask = 3;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.NV_mdl = 'bayesianModel_3';
            GLMprm.choice.(Epm_nm).RP.E.NV_chosen = 3;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
        %% class of GLMs only modeling choice and perf
        % + temporal derivative
    case 229
        % general parameters
        GLMprm.gal.add_drv = 1; % temporal derivative
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        GLMprm.gal.grey_mask = 3;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'boxcar';
            GLMprm.choice.(Epm_nm).RP.E.NV_mdl = 'bayesianModel_3';
            GLMprm.choice.(Epm_nm).RP.E.NV_chosen = 3;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'boxcar';
        end % physical/mental loop
        
    case 230 % like GLM 229 but without SV
        % general parameters
        GLMprm.gal.add_drv = 1; % temporal derivative
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        GLMprm.gal.grey_mask = 3;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'boxcar';
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'boxcar';
        end % physical/mental loop
        
    case 231 % same as GLM 229 but without boxcar
        % general parameters
        GLMprm.gal.add_drv = 1; % temporal derivative
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        GLMprm.gal.grey_mask = 3;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.NV_mdl = 'bayesianModel_3';
            GLMprm.choice.(Epm_nm).RP.E.NV_chosen = 3;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
        end % physical/mental loop
        
    case 232 % same as GLM 230 but without boxcar (or like GLM 231 without SV)
        % general parameters
        GLMprm.gal.add_drv = 1; % temporal derivative
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        GLMprm.gal.grey_mask = 3;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
        end % physical/mental loop
        
    case 233 % same as GLM 231 but without temporal derivative
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        GLMprm.gal.grey_mask = 3;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.NV_mdl = 'bayesianModel_3';
            GLMprm.choice.(Epm_nm).RP.E.NV_chosen = 3;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
        end % physical/mental loop
        
    case 234 % same as GLM 230 but without temporal derivative (or like GLM 233 without SV)
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        GLMprm.gal.grey_mask = 3;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
        end % physical/mental loop
        
    case 235 % same as GLM 214 but removing time effect Fp/Fm
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        GLMprm.gal.grey_mask = 3;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % fixation crosses
            GLMprm.model_onset.(Epm_nm).allCrosses = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.NV_mdl = 'bayesianModel_3';
            GLMprm.choice.(Epm_nm).RP.E.NV_chosen = 3;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 236 % choice/hR/hP/hE/Fp/Fm/Ech
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 0;
        GLMprm.gal.grey_mask = 3;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % fixation crosses
            GLMprm.model_onset.(Epm_nm).allCrosses = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.choiceHighE = 2;
            GLMprm.choice.(Epm_nm).RP.E.R_varOption = 2;
            GLMprm.choice.(Epm_nm).RP.E.P_varOption = 2;
            GLMprm.choice.(Epm_nm).RP.E.E_varOption = 1;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 3;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 2;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 7;
            end
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 237 % same as GLM 214 but split pre-choice and pre-effort cross
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        GLMprm.gal.grey_mask = 3;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % pre-choice fixation cross
            GLMprm.model_onset.(Epm_nm).preChoiceCross = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.NV_mdl = 'bayesianModel_3';
            GLMprm.choice.(Epm_nm).RP.E.NV_chosen = 3;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % pre-effort fixation cross
            GLMprm.model_onset.(Epm_nm).preEffortCross = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 238 % same as GLM 214 but E*choice instead of Echosen
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        GLMprm.gal.grey_mask = 3;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % fixation crosses
            GLMprm.model_onset.(Epm_nm).allCrosses = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.NV_mdl = 'bayesianModel_3';
            GLMprm.choice.(Epm_nm).RP.E.NV_chosen = 3;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 3;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            GLMprm.choice.(Epm_nm).RP.E.RT = 1;
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
        
    case 239 % same as GLM 214 but without RT
        % general parameters
        GLMprm.gal.orth_vars = 0;
        GLMprm.gal.zPerRun = 1;
        GLMprm.gal.grey_mask = 3;
        GLMprm.gal.mask_probaThreshold = 5;
        % loop per task
        for iEpm = 1:length(Epm)
            Epm_nm = Epm{iEpm};
            % fixation crosses
            GLMprm.model_onset.(Epm_nm).allCrosses = 'stick';
            % choice
            GLMprm.model_onset.(Epm_nm).choice = 'stick';
            GLMprm.choice.(Epm_nm).RP.E.NV_mdl = 'bayesianModel_3';
            GLMprm.choice.(Epm_nm).RP.E.NV_chosen = 3;
            GLMprm.choice.(Epm_nm).RP.E.E_chosen = 1;
            switch Epm_nm
                case 'Ep'
                    GLMprm.choice.(Epm_nm).RP.E.fatigue = 1;
                case 'Em'
                    GLMprm.choice.(Epm_nm).RP.E.prevEfficacy = 3;
            end
            % chosen
            GLMprm.model_onset.(Epm_nm).chosen = 'stick';
            % effort perf (effort execution)
            GLMprm.model_onset.(Epm_nm).Eperf = 'stick';
            % feedback
            GLMprm.model_onset.(Epm_nm).fbk = 'stick';
        end % physical/mental loop
end % GLM number
%% warnings: check compatibility of the GLM parameters entered
isGLMokCheck(GLMprm);

end % function