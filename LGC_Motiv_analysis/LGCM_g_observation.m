function [gx] = LGCM_g_observation(x,P,u,inG)
% [gx] = LGCM_g_observation(x,P,u,inG)
% LGCM_g_observation= observation function for 
%
% INPUTS
% x: hidden states (could be fatigue if fatigue is entered as a hidden
% state)
%
% P: parameters to estimate (sensitivity to reward, effort, fatigue)
%
% u: task parameters
%
% inG: structure with additional variables of interest, in particular with
% model parameters stored in inG.mdlPrm
%
% OUTPUTS
% gx: output of the observation function, a binary variable corresponding
% to the choice (equal to (1) if the choice made was the left option, and
% to (0) if the choice made was the right option)

%% extract conditions of the model
fatigueEvolOrObs    = inG.mdlPrm.fatigueEvolOrObs;
kRP_tasksSplitOrPool = inG.mdlPrm.kRP_tasksSplitOrPool;
fatigue = inG.mdlPrm.fatigue;

%% load parameters to estimate
% reward
iP = 1;
kR = P(iP);
% punishment
iP = iP + 1;
kP = P(iP);
% effort and fatigue
switch fatigueEvolOrObs
    case 'split'
        % effort
        iP = iP + 1;
        kE = P(iP);
        
        % fatigue
        switch fatigueEvolOrObs
            case 'observation parameter'
                iP = iP + 1;
                kF = P(iP);
        end
        
    case 'pool'
        % effort mental
        iP = iP + 1;
        kEm = P(iP);
        
        % effort physical
        iP = iP + 1;
        kEp = P(iP);
        
        % fatigue
        switch fatigueEvolOrObs
            case 'observation parameter'
                % fatigue physical
                iP = iP + 1;
                kFEm = P(iP);
                % fatigue mental
                iP = iP + 1;
                kFEp = P(iP);
        end
end

% note: possible to pool mental and physical effort tasks together to
% estimate one single kR parameter for both tasks while estimating kE and
% kF separately for each

%% load task parameters
trialN      = u(1);
Rleft       = u(2);
Rright      = u(3);
switch kRP_tasksSplitOrPool
    case 'split'
        Eleft       = u(4);
        Eright      = u(5);
        sumEprev    = u(6);
        intEprev    = u(7);
        % ElevelPrevTrial = u(8);
        % EintegralPrevTrial = u(9);
        R_or_P = u(10);
    case 'pool'
        Eleft_Em     = u(4);
        Eright_Em    = u(5);
        Eleft_Ep     = u(6);
        Eright_Ep    = u(7);
        sumEprev_Em  = u(8);
        intEprev_Em  = u(9);
        sumEprev_Ep  = u(10);
        intEprev_Ep  = u(11);
        effort_type  = u(12);
        % ElevelPrevTrial_Em = u(13);
        % ElevelPrevTrial_Ep = u(14);
        % EintegralPrevTrial_Em = u(15);
        % EintegralPrevTrial_Ep = u(16);
        R_or_P = u(17);
end

%% compute variables of interest
DeltaR = Rleft - Rright;
switch R_or_P
    case 0 % punishment
        kRP = kP;
    case 1 % reward
        kRP = kR;
end
switch kRP_tasksSplitOrPool
    case 'pool'
        switch effort_type
            case 0 % mental
                kE = kEm;
                DeltaE = Eleft_Em - Eright_Em;
            case 1 % physical
                kE = kEp;
                DeltaE = Eleft_Ep - Eright_Ep;
        end
    case 'split'
        DeltaE = Eleft - Eright;
end
% define fatigue component to consider
switch fatigueEvolOrObs
    case 'hidden state'
        switch kRP_tasksSplitOrPool
            case 'split'
                FatigueComponent = x(1);
            case 'pool'
                FatigueComponent_Em = x(1);
                FatigueComponent_Ep = x(2);
                switch effort_type
                    case 0 %mental
                        FatigueComponent = FatigueComponent_Em;
                        kF = kFEm;
                    case 1 % physical
                        FatigueComponent = FatigueComponent_Ep;
                        kF = kFEp;
                end
        end
    case 'observation parameter'
        switch kRP_tasksSplitOrPool
            case 'split'
                switch fatigue
                    case 'trialN' % use trial number as as simple account of the cumulated level of fatigue
                        FatigueComponent = trialN;
                    case 'sumEffort' % use the sum of all the efforts done until now
                        FatigueComponent = sumEprev;
                    case 'integralEffort' % integral of the effort performed instead of sum of effort levels
                        % mental effort: sum of the number of correct answers provided
                        % physical effort: integral of the efforts performed until the
                        % current trial
                        FatigueComponent = intEprev;
                end
            case 'pool'
                switch fatigue
                    case 'trialN' % use trial number as as simple account of the cumulated level of fatigue
                        FatigueComponent = trialN;
                    case 'sumEffort' % use the sum of all the efforts done until now
                        switch effort_type
                            case 0 % mental
                                FatigueComponent = sumEprev_Em;
                            case 1 % physical
                                FatigueComponent = sumEprev_Ep;
                        end
                    case 'integralEffort' % integral of the effort performed instead of sum of effort levels
                        % mental effort: sum of the number of correct answers provided
                        % physical effort: integral of the efforts performed until the
                        % current trial
                        switch effort_type
                            case 0 % mental
                                FatigueComponent = intEprev_Em;
                            case 1 % physical
                                FatigueComponent = intEprev_Ep;
                        end
                end
        end
end

%% output
switch fatigueEvolOrObs
    case 'hidden state' % fatigue sensitivity computed in evolution function f already
        gx = 1./( 1 + exp( -( (kRP.*DeltaR) -kE.*DeltaE.*(1+FatigueComponent) ) ));
    case 'observation parameter'
        gx = 1./( 1 + exp( -( (kRP.*DeltaR) -kE.*DeltaE.*(1+kF.*FatigueComponent) ) ));
end

end % function