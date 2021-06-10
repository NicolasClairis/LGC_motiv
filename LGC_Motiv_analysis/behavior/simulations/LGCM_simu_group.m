function [] = LGCM_simu_group()


%% number of subjects for the simulation
NS_simu = 75;

%% define parameters of the simulation
n_RP_levels = 3;
n_E_levels = 3;
n_trials = 44;

mdlPrm.fatigueEvolOrObs = 'observation parameter'; % 'hidden state','observation parameter'
mdlPrm.fatigueType = 'trial number'; % 'trial number','sum of efforts','effort integral'
if ~strcmp(mdlPrm.fatigueType,'trial number')
    error('impossible yet to simulate other types of fatigue effects. Would need to include sum of previous efforts in the evolution function to make it work');
end
mdlPrm.kRP_tasksSplitOrPool = 'split'; % 'pool','split'

%% initialize parameters of interest
[choice_left] = deal(NaN(n_trials, NS_simu));

%% simulate the data
for iS = 1:NS_simu
    
    %% define priors
    switch fatigueEvolOrObs
        case 'hidden state'
            theta = ;
        case 'observation parameter'
            theta = []; % no evolution function
    end
    
    %% prepare inputs
    trialN = 1:n_trials;
    punishment_yn = 'yes';
    choiceOptions = choice_option_design(n_RP_levels, n_E_levels, punishment_yn, n_trials);
    Rleft = choiceOptions.R.left;
    Rright = choiceOptions.R.right;
    Eleft = choiceOptions.E.left;
    Eright = choiceOptions.E.right;
    R_or_P = choiceOptions.R_or_P;
    sumEprev = [];
    intEprev = [];
    ElevelPrevTrial = [];
    EintegralPrevTrial = [];
    switch kRP_tasksSplitOrPool
        case 'pool'
            %             u_t = [trialN;...
            %                 Rleft; Rright;...
            %                 EleftPool_Em; ErightPool_Em;...
            %                 EleftPool_Ep; ErightPool_Ep;...
            %                 sumEprevPool_Em; intEprevPool_Em;...
            %                 sumEprevPool_Ep; intEprevPool_Ep;...
            %                 effort_type;...
            %                 ElevelPrevTrialPool_Em; ElevelPrevTrialPool_Ep;...
            %                 EintegralPrevTrialPool_Em; EintegralPrevTrialPool_Ep;...
            %                 R_or_P];
            error('need to work a bit if you want to pool the data of the two tasks together');
        case 'split'
            u_t = [trialN;...
                Rleft; Rright;...
                Eleft; Eright;...
                sumEprev; intEprev;...
                ElevelPrevTrial,...
                EintegralPrevTrial;...
                R_or_P];
    end
    
    %% simulate the choices
    [choice_left(:,iS)] = LGCM_simu_generate_data(theta, phi, u_t,...
        n_trials, mdlPrm);
end % subject loop

end % function