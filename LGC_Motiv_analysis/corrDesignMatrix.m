function [m_R, sd_R, bestMatrix] = corrDesignMatrix()
% [m_R, sd_R, bestMatrix] = corrDesignMatrix()
% corrDesignMatrix will perform the correlation between the variables of
% the design matrix (delta reward and delta effort and also between
% rewards, efforts and trial number).
%
% OUTPUTS
% m_R: structure with mean correlation coefficient for each correlation
% tested
%
% sd_R: structure with standard deviation correlation coefficient for each
% correlation tested
%
% bestMatrix: matrix with the smaller sum of correlation coefficient for
% all possible correlations

%% initialize parameters of interest
n_RP_levels = 3;
punishment_yn = 'yes';
n_E_levels = 3;
n_trials = 44;
[R_money] = R_amounts(n_RP_levels, punishment_yn);

% how many design matrices do you want to test?
n_designMatrix = 1000;
trialN_vec = 1:n_trials;

[Rcorr_deltaR_vs_deltaE,...
    Rcorr_deltaP_vs_deltaE,...
    Rcorr_deltaR_vs_trialN,...
    Rcorr_deltaP_vs_trialN,...
    Rcorr_deltaE_vs_trialN_R,...
    Rcorr_deltaE_vs_trialN_P,...
    Rcorr_RPtrials_vs_trialN,...
    Rcorr_global] = deal( NaN(1, n_designMatrix) );

%% run the tests
for iTest = 1:n_designMatrix
    choice_opt_tmp = choice_option_design(n_RP_levels, n_E_levels, punishment_yn, n_trials, R_money);
    R_trials = strcmp(choice_opt_tmp.R_or_P,'R');
    P_trials = strcmp(choice_opt_tmp.R_or_P,'P');
    deltaR = choice_opt_tmp.R.left(R_trials)' - choice_opt_tmp.R.right(R_trials)';
    deltaP = choice_opt_tmp.R.left(P_trials)' - choice_opt_tmp.R.right(P_trials)';
    deltaE_R = choice_opt_tmp.E.left(R_trials)' - choice_opt_tmp.E.right(R_trials)';
    deltaE_P = choice_opt_tmp.E.left(P_trials)' - choice_opt_tmp.E.right(P_trials)';
    trialN_R = trialN_vec(R_trials)';
    trialN_P = trialN_vec(P_trials)';
    % store matrix
    choice_opt.(['matrix_',num2str(iTest)]) = choice_opt_tmp;
    
    % test correlations
    RE_matrix_tmp = corrcoef(deltaR, deltaE_R);
    Rcorr_deltaR_vs_deltaE(iTest)   = RE_matrix_tmp(1,2);
    PE_matrix_tmp = corrcoef(deltaP, deltaE_P);
    Rcorr_deltaP_vs_deltaE(iTest)   = PE_matrix_tmp(1,2);
    Rtrial_matrix_tmp = corrcoef(deltaR, trialN_R);
    Rcorr_deltaR_vs_trialN(iTest)   = Rtrial_matrix_tmp(1,2);
    Ptrial_matrix_tmp = corrcoef(deltaP, trialN_P);
    Rcorr_deltaP_vs_trialN(iTest)   = Ptrial_matrix_tmp(1,2);
    Etrial_matrix_tmp = corrcoef(deltaE_R, trialN_R);
    Rcorr_deltaE_vs_trialN_R(iTest) = Etrial_matrix_tmp(1,2);
    Etrial_matrix_tmp = corrcoef(deltaE_P, trialN_P);
    Rcorr_deltaE_vs_trialN_P(iTest) = Etrial_matrix_tmp(1,2);
    RPtrialN_matrix_tmp = corrcoef(R_trials, trialN_vec);
    Rcorr_RPtrials_vs_trialN(iTest) = RPtrialN_matrix_tmp(1,2);
    
    % check sum of correlations
    Rcorr_global(iTest) = Rcorr_deltaR_vs_deltaE(iTest) + Rcorr_deltaP_vs_deltaE(iTest) +...
        Rcorr_deltaR_vs_trialN(iTest) + Rcorr_deltaP_vs_trialN(iTest) +...
        Rcorr_deltaE_vs_trialN_R(iTest) + Rcorr_deltaE_vs_trialN_P(iTest) +...
        Rcorr_RPtrials_vs_trialN(iTest);
end % test loop

%% compute average and SD of correlation coefficient
% mean
m_R.Rcorr_deltaR_vs_deltaE     = mean(Rcorr_deltaR_vs_deltaE, 2);
m_R.Rcorr_deltaP_vs_deltaE     = mean(Rcorr_deltaP_vs_deltaE, 2);
m_R.Rcorr_deltaR_vs_trialN     = mean(Rcorr_deltaR_vs_trialN, 2);
m_R.Rcorr_deltaP_vs_trialN     = mean(Rcorr_deltaP_vs_trialN, 2);
m_R.Rcorr_deltaE_vs_trialN_R   = mean(Rcorr_deltaE_vs_trialN_R, 2);
m_R.Rcorr_deltaE_vs_trialN_P   = mean(Rcorr_deltaE_vs_trialN_P, 2);
m_R.Rcorr_RPtrials_vs_trialN   = mean(Rcorr_RPtrials_vs_trialN, 2);
% SD
sd_R.Rcorr_deltaR_vs_deltaE    = std(Rcorr_deltaR_vs_deltaE);
sd_R.Rcorr_deltaP_vs_deltaE    = std(Rcorr_deltaP_vs_deltaE);
sd_R.Rcorr_deltaR_vs_trialN    = std(Rcorr_deltaR_vs_trialN);
sd_R.Rcorr_deltaP_vs_trialN    = std(Rcorr_deltaP_vs_trialN);
sd_R.Rcorr_deltaE_vs_trialN_R  = std(Rcorr_deltaE_vs_trialN_R);
sd_R.Rcorr_deltaE_vs_trialN_P  = std(Rcorr_deltaE_vs_trialN_P);
sd_R.Rcorr_RPtrials_vs_trialN  = std(Rcorr_RPtrials_vs_trialN);

%% compute the matrix for which all correlations are minimized
bestMatrix_idx = find( min(Rcorr_global) );
bestMatrix = choice_opt.(['matrix_',num2str(bestMatrix_idx)]);
end % function