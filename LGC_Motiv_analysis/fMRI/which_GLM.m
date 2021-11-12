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
%       .add_drv: derivative (0: no derivative added; 1: temporal
%       derivative; 2: temporal AND spatial derivative)
%       .grey_mask: add a grey mask (0: include all voxels; 1: filter based
%       on individual grey matter; 2: use SPM template grey matter mask)
%       .orth_vars: orthogonalize variables of the GLM (1) or not (0)
%   .model_onset: indicate for each task (Ep/Em: physical/mental) for each
%   event (preChoiceCross/choice/chosen/preEffortCross/Eperf/fbk) if it should be modelled as a
%   stick ('stick') as a boxcar ('boxcar') or not included in the GLM
%   ('none')
%   .choice/chosen/Eperf/fbk: for each event and each task (Ep/Em) indicate
%   if a given regressor should be included or not
%

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
        [GLMprm.choice.(EpEm_nm).(RP_nm).money_left,...
            GLMprm.choice.(EpEm_nm).(RP_nm).money_right,...
            GLMprm.choice.(EpEm_nm).(RP_nm).money_chosen,...
            GLMprm.choice.(EpEm_nm).(RP_nm).money_unchosen,...
            GLMprm.choice.(EpEm_nm).(RP_nm).money_ch_min_unch,...
            GLMprm.choice.(EpEm_nm).(RP_nm).money_sum,...
            GLMprm.choice.(EpEm_nm).(RP_nm).money_varOption,...
            GLMprm.choice.(EpEm_nm).(RP_nm).E_left,...
            GLMprm.choice.(EpEm_nm).(RP_nm).E_right,...
            GLMprm.choice.(EpEm_nm).(RP_nm).E_chosen,...
            GLMprm.choice.(EpEm_nm).(RP_nm).E_unchosen,...
            GLMprm.choice.(EpEm_nm).(RP_nm).E_ch_min_unch,...
            GLMprm.choice.(EpEm_nm).(RP_nm).E_sum,...
            GLMprm.choice.(EpEm_nm).(RP_nm).E_varOption,...
            GLMprm.choice.(EpEm_nm).(RP_nm).confidence,...
            GLMprm.choice.(EpEm_nm).(RP_nm).RT,...
            GLMprm.choice.(EpEm_nm).(RP_nm).R_vs_P] = deal(0);
        
        % chosen option display
        % pool reward and punishment together (default)
        GLMprm.chosen.(EpEm_nm).RPpool = 1;
        [GLMprm.chosen.(EpEm_nm).(RP_nm).money_chosen,...
            GLMprm.chosen.(EpEm_nm).(RP_nm).E_chosen,...
            GLMprm.chosen.(EpEm_nm).(RP_nm).money_unchosen,...
            GLMprm.chosen.(EpEm_nm).(RP_nm).E_unchosen,...
            GLMprm.chosen.(EpEm_nm).(RP_nm).money_varOption,...
            GLMprm.chosen.(EpEm_nm).(RP_nm).E_varOption,...
            GLMprm.chosen.(EpEm_nm).(RP_nm).confidence] = deal(0);
        
        % effort performance
        % pool reward and punishment together (default)
        GLMprm.Eperf.(EpEm_nm).RPpool = 1;
        [GLMprm.Eperf.(EpEm_nm).(RP_nm).money_chosen,...
            GLMprm.Eperf.(EpEm_nm).(RP_nm).E_chosen,...
            GLMprm.Eperf.(EpEm_nm).(RP_nm).RT_1stAnswer] = deal(0);
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
        
        % feedback
        % pool reward and punishment together (default)
        GLMprm.fbk.(EpEm_nm).RPpool = 1;
        [GLMprm.fbk.(EpEm_nm).(RP_nm).money_obtained,...
            GLMprm.fbk.(EpEm_nm).(RP_nm).win_vs_loss,...
            GLMprm.fbk.(EpEm_nm).(RP_nm).E_made,...
            GLMprm.fbk.(EpEm_nm).(RP_nm).confidence] = deal(0);
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
        % choice
        GLMprm.model_onset.Ep.choice = 'stick';
        GLMprm.model_onset.Em.choice = 'stick';
        GLMprm.choice.Ep.RP.R_vs_P = 1;
        GLMprm.choice.Em.RP.R_vs_P = 1;
        GLMprm.choice.Ep.RP.money_chosen = 1;
        GLMprm.choice.Em.RP.money_chosen = 1;
        GLMprm.choice.Ep.RP.E_chosen = 1;
        GLMprm.choice.Em.RP.E_chosen = 1;
        % effort perf
        GLMprm.model_onset.Ep.Eperf = 'stick';
        GLMprm.model_onset.Em.Eperf = 'stick';
    case 5
        % general parameters
        GLMprm.gal.orth_vars = 1;
        % choice
        GLMprm.model_onset.Ep.choice = 'stick';
        GLMprm.model_onset.Em.choice = 'stick';
        GLMprm.choice.Ep.RP.R_vs_P = 1;
        GLMprm.choice.Em.RP.R_vs_P = 1;
        % effort perf
        GLMprm.model_onset.Ep.Eperf = 'stick';
        GLMprm.model_onset.Em.Eperf = 'stick';
    case 6 % testing Vchosen - Vunchosen during choice + adding temporal derivative
        % general parameters
        GLMprm.gal.orth_vars = 1;
        GLMprm.gal.add_drv = 1;
        % choice
        GLMprm.model_onset.Ep.choice = 'stick';
        GLMprm.model_onset.Em.choice = 'stick';
        GLMprm.choice.Ep.RP.money_ch_min_unch = 1;
        GLMprm.choice.Em.RP.money_ch_min_unch = 1;
        GLMprm.choice.Ep.RP.E_ch_min_unch = 1;
        GLMprm.choice.Em.RP.E_ch_min_unch = 1;
        % chosen
        GLMprm.model_onset.Ep.chosen = 'stick';
        GLMprm.model_onset.Em.chosen = 'stick';
        % effort perf
        GLMprm.model_onset.Ep.Eperf = 'stick';
        GLMprm.model_onset.Em.Eperf = 'stick';
        % feedback
        GLMprm.model_onset.Ep.fbk = 'stick';
        GLMprm.model_onset.Em.fbk = 'stick';
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
    case 10 % model money amounts during performance instead of choice periode (VS should be trigger by higher rewards)
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
end


end % function