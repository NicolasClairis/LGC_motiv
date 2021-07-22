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
%

%% initialize all variables

% general parameters: derivative, grey matter filtering
[GLMprm.gal.add_drv,...
    GLMprm.gal.grey_mask,...
    GLMprm.gal.orth_vars] = deal(0);

% onsets: not modelled (none), modelled as stick function (stick) or as
% boxcar function (boxcar)
[GLMprm.model_onset.Ep.cross,...
    GLMprm.model_onset.Ep.choice,...
    GLMprm.model_onset.Ep.chosen,...
    GLMprm.model_onset.Ep.Eperf,...
    GLMprm.model_onset.Ep.fbk,...
    GLMprm.model_onset.Em.cross,...
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
            GLMprm.choice.(EpEm_nm).(RP_nm).E_left,...
            GLMprm.choice.(EpEm_nm).(RP_nm).E_right,...
            GLMprm.choice.(EpEm_nm).(RP_nm).E_chosen,...
            GLMprm.choice.(EpEm_nm).(RP_nm).E_unchosen,...
            GLMprm.choice.(EpEm_nm).(RP_nm).E_ch_min_unch,...
            GLMprm.choice.(EpEm_nm).(RP_nm).E_sum,...
            GLMprm.choice.(EpEm_nm).(RP_nm).RT] = deal(0);
        
        % chosen option display
        % pool reward and punishment together (default)
        GLMprm.chosen.(EpEm_nm).RPpool = 1;
        [GLMprm.chosen.(EpEm_nm).(RP_nm).money_chosen,...
            GLMprm.chosen.(EpEm_nm).(RP_nm).E_chosen,...
            GLMprm.chosen.(EpEm_nm).(RP_nm).certainty] = deal(0);
        
        % effort performance
        % pool reward and punishment together (default)
        GLMprm.Eperf.(EpEm_nm).RPpool = 1;
        [GLMprm.Eperf.(EpEm_nm).(RP_nm).money,...
            GLMprm.Eperf.(EpEm_nm).(RP_nm).effort,...
            GLMprm.Eperf.(EpEm_nm).(RP_nm).RT_1stAnswer] = deal(0);
        % specific variables for each effort type
        switch EpEm_nm
            case 'Ep'
                % effort performance
                [GLMprm.Eperf.Ep.(RP_nm).Fpeak,...
                    GLMprm.Eperf.Ep.(RP_nm).Fintegral] = deal(0);
            case 'Em'
                % effort performance
                [GLMprm.Eperf.Em.(RP_nm).RT_avg,...
                    GLMprm.Eperf.Em.(RP_nm).n_errors] = deal(0);
        end
        
        % feedback
        % pool reward and punishment together (default)
        GLMprm.fbk.(EpEm_nm).RPpool = 1;
        [GLMprm.fbk.(EpEm_nm).(RP_nm).money_obtained,...
            GLMprm.fbk.(EpEm_nm).(RP_nm).E_made,...
            GLMprm.fbk.(EpEm_nm).(RP_nm).certainty] = deal(0);
    end % RP
end % effort type
    
%% define variables according to GLM number
switch GLM
    case 1
        % general parameters
        GLMprm.gal.orth_vars = 1;
        % cross
        GLMprm.model_onset.Ep.cross = 'stick';
        GLMprm.model_onset.Em.cross = 'stick';
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
end


end % function