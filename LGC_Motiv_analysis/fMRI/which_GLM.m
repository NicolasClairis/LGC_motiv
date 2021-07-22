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
[GLMprm.add_drv,...
    GLMprm.grey_mask] = deal(0);

% onsets: not modelled (none), modelled as stick function (stick) or as
% boxcar function (boxcar)
[GLMprm.model_onset.Ep.cross,...
    GLMprm.model_onset.Ep.choice,...
    GLMprm.model_onset.Ep.chosen,...
    GLMprm.model_onset.Ep.Eperf,...
    GLMprm.model_onset.Em.cross,...
    GLMprm.model_onset.Em.choice,...
    GLMprm.model_onset.Em.chosen,...
    GLMprm.model_onset.Em.Eperf] = deal('none');


% initialize all modulators of interest at 0
%
% use same conditions for physical and mental effort and pool reward and
% punishment trials or split reward and punishment trials
Ep_Em_Epm = {'Ep','Em','Epm'}; % apply same or different conditions to physical and mental effort
for iEpm = 1:length(EP_Em_Epm)
    EpEm_nm = Ep_Em_Epm{iEpm};
    RPconditions = {'R','P','RP'};
    GLMprm.choice.Epm.RPpool = 1; % pool reward and punishment together (default)
    for iRP = 1:length(RPconditions)
        RP_nm = RPconditions{iRP};
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
            GLMprm.choice.(EpEm_nm).(RP_nm).E_sum] = deal(0);
        
        % chosen option display
        [GLMprm.chosen.(EpEm_nm).money_chosen,...
            GLMprm.chosen.(EpEm_nm).E_chosen,...
            GLMprm.chosen.(EpEm_nm).certainty] = deal(0);
        
        % effort performance
        [GLMprm.Eperf.(EpEm_nm).money,...
            GLMprm.Eperf.(EpEm_nm).effort_level,...
            GLMprm.Eperf.(EpEm_nm).RT_1stAnswer] = deal(0);
        % specific variables for each effort type
        switch EpEm_nm
            case 'Ep'
                % effort performance
                [GLMprm.Eperf.Ep.Fpeak,...
                    GLMprm.Eperf.Ep.Fintegral] = deal(0);
            case 'Em'
                % effort performance
                [GLMprm.Eperf.Em.RT_avg,...
                    GLMprm.Eperf.Em.n_errors] = deal(0);
        end
    end % RP
end % effort type
    
%% define variables according to GLM number
switch GLM
    case 1
        GLMprm.model_onset.cross = 'stick';
        GLMprm.model_onset.choice = 'stick';
        GLMprm.choice.Epm.RP.money_sum = 1;
        GLMprm.choice.Epm.RP.E_sum = 1;
        GLMprm.model_onset.chosen = 'stick';
        GLMprm.chosen.Epm.RP.money_chosen = 1;
        GLMprm.chosen.(EpEm_nm).E_chosen = 1;
        GLMprm.model_onset.Eperf = 'stick';
end


end % function