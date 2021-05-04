function [ fx ] = LGCM_f_evolution( x, P, u, inF )
%[ fx ] = LGCM_f_evolution( x, P, u, inF )
%
% INPUTS
% x: hidden states = fatigue component
%
% P: parameters to estimate
%
% u: task parameters
%
% inF: indication about which model to use
%
% OUTPUTS
% fx: updated hidden state for the current trial

%% load model parameters
kR_tasksSplitOrPool = inF.mdlPrm.kR_tasksSplitOrPool;
fatigueType = inF.mdlPrm.fatigueType;

%% load hidden states
switch kR_tasksSplitOrPool
    case 'split'
        fatigue_lastTrial   = x(1);
    case 'pool' % both tasks are pooled but fatigue parameter should be estimated independently
        fatigue_Em_lastTrial  = x(1);
        fatigue_Ep_lastTrial  = x(2);
end

%% load parameters to estimate
switch kR_tasksSplitOrPool
    case 'split'
        iP = 1;
        kF = P(iP);
    case 'pool'
        iP = 1;
        kF_Em = P(iP);
        iP = iP + 1;
        kF_Ep = P(iP);
end

%% load task variables
switch kR_tasksSplitOrPool
    case 'split'
        trialN = u(1);
        % Rleft = u(2);
        % Rright = u(3);
        % Eleft = u(4);
        % Eright = u(5);
        % sumEprev = u(6);
        % intEprev = u(7);
        ElevelPrevTrial = u(8);
        EintegralPrevTrial = u(9);
    case 'pool'
        trialN = u(1);
        % Rleft = u(2);
        % Rright = u(3);
        % EleftPool_Em = u(4);
        % ErightPool_Em = u(5);
        % EleftPool_Ep = u(6);
        % ErightPool_Ep = u(7);
        % sumEprevPool_Em = u(8);
        % intEprevPool_Em = u(9);
        % sumEprevPool_Ep = u(10);
        % intEprevPool_Ep = u(11);
        effort_type = u(12);
        ElevelPrevTrialPool_Em = u(13);
        ElevelPrevTrialPool_Ep = u(14);
        EintegralPrevTrialPool_Em = u(15);
        EintegralPrevTrialPool_Ep = u(16);
end

% extract effort based on the current model
switch kR_tasksSplitOrPool
    case 'split'
        switch fatigueType
            case 'trial number'
                effortLastTrial = trialN;
            case 'sum of efforts'
                effortLastTrial = ElevelPrevTrial;
            case 'effort integral'
                effortLastTrial = EintegralPrevTrial;
        end
    case 'pool'
        switch fatigueType
            case 'trial number'
                effortLastTrial = trialN;
            case 'sum of efforts'
                switch effort_type
                    case 0 % mental
                        effortLastTrial = ElevelPrevTrialPool_Em;
                    case 1 % physical
                        effortLastTrial = ElevelPrevTrialPool_Ep;
                end
            case 'effort integral'
                switch effort_type
                    case 0 % mental
                        effortLastTrial = EintegralPrevTrialPool_Em;
                    case 1 % physical
                        effortLastTrial = EintegralPrevTrialPool_Ep;
                end
        end
end

%% update fatigue according to the model used
switch kR_tasksSplitOrPool
    case 'split'
        fatigue = fatigue_lastTrial + kF*effortLastTrial;
    case 'pool'
        switch effort_type
            case 0 % mental effort
                fatigue_Em = fatigue_Em_lastTrial + kF_Em*effortLastTrial;
                fatigue_Ep = 0;
            case 1 % physical effort
                fatigue_Em = 0;
                fatigue_Ep = fatigue_Ep_lastTrial + kF_Ep*effortLastTrial;
        end
end
%% output
switch kR_tasksSplitOrPool
    case 'split'
        fx = fatigue;
    case 'pool'
        fx = [fatigue_Em, fatigue_Ep];
end

end % function