function [Rcorr_bestMatrix, m_R, sd_R, bestMatrix, Rcorr] = corrDesignMatrix_ter()
% [Rcorr_bestMatrix, m_R, sd_R, bestMatrix, Rcorr] = corrDesignMatrix_ter()
% corrDesignMatrix_bis will perform the correlation between the variables of
% the design matrix (delta reward and delta effort and also between
% rewards, efforts and the effort*trial number interaction).
% Note: corrDesignMatrix checks the correlation with trial number while corrDesignMatrix_bis
% will check the correlation with the interaction between deltaE and trial
% number.
%
% OUTPUTS
% Rcorr_bestMatrix: structure with correlation coefficients for the best
% matrix
%
% m_R: structure with mean correlation coefficient for each correlation
% tested
%
% sd_R: structure with standard deviation correlation coefficient for each
% correlation tested
%
% bestMatrix: matrix with the smaller sum of correlation coefficient for
% all possible correlations
%
% Rcorr: structure with all correlation coefficients

%% initialize parameters of interest
n_RP_levels = 3;
punishment_yn = 'yes';
n_E_levels = 3;
nTrials = 54;

% how many design matrices do you want to test?
n_designMatrix = 10000;
trialN_vec = 1:nTrials;

[Rcorr_deltaR_vs_deltaE,...
    Rcorr_deltaP_vs_deltaE,...
    Rcorr_deltaE_vs_deltaEtrialN_R,...
    Rcorr_deltaE_vs_deltatrialN_P,...
    Rcorr_RPtrials_vs_trialN,...
    Rcorr_global] = deal( NaN(1, n_designMatrix) );

%% run the tests
for iTest = 1:n_designMatrix
    choice_opt_tmp = design_trialOptions(n_RP_levels, n_E_levels, punishment_yn, nTrials);
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
    Etrial_R_matrix_tmp = corrcoef(deltaE_R, trialN_R.*deltaE_R);
    Rcorr_deltaE_vs_deltaEtrialN_R(iTest) = Etrial_R_matrix_tmp(1,2);
    Etrial_P_matrix_tmp = corrcoef(deltaE_P, trialN_P.*deltaE_P);
    Rcorr_deltaE_vs_deltatrialN_P(iTest) = Etrial_P_matrix_tmp(1,2);
    RPtrialN_matrix_tmp = corrcoef(R_trials, trialN_vec);
    Rcorr_RPtrials_vs_trialN(iTest) = RPtrialN_matrix_tmp(1,2);
    
    % check sum of correlations
    Rcorr_global(iTest) = Rcorr_deltaR_vs_deltaE(iTest)^2 + Rcorr_deltaP_vs_deltaE(iTest)^2 +...
        Rcorr_deltaE_vs_deltaEtrialN_R(iTest)^2 + Rcorr_deltaE_vs_deltatrialN_P(iTest)^2 +...
        Rcorr_RPtrials_vs_trialN(iTest)^2;
end % test loop

%% compute average and SD of correlation coefficient
% mean
m_R.Rcorr_deltaR_vs_deltaE     = mean(Rcorr_deltaR_vs_deltaE, 2);
m_R.Rcorr_deltaP_vs_deltaE     = mean(Rcorr_deltaP_vs_deltaE, 2);
m_R.Rcorr_deltaE_vs_deltaEtrialN_R   = mean(Rcorr_deltaE_vs_deltaEtrialN_R, 2);
m_R.Rcorr_deltaE_vs_deltatrialN_P   = mean(Rcorr_deltaE_vs_deltatrialN_P, 2);
m_R.Rcorr_RPtrials_vs_trialN   = mean(Rcorr_RPtrials_vs_trialN, 2);
% SD
sd_R.Rcorr_deltaR_vs_deltaE    = std(Rcorr_deltaR_vs_deltaE);
sd_R.Rcorr_deltaP_vs_deltaE    = std(Rcorr_deltaP_vs_deltaE);
sd_R.Rcorr_deltaE_vs_deltaEtrialN_R  = std(Rcorr_deltaE_vs_deltaEtrialN_R);
sd_R.Rcorr_deltaE_vs_deltatrialN_P  = std(Rcorr_deltaE_vs_deltatrialN_P);
sd_R.Rcorr_RPtrials_vs_trialN  = std(Rcorr_RPtrials_vs_trialN);

% store all correlation matrices
Rcorr.Rcorr_deltaR_vs_deltaE = Rcorr_deltaR_vs_deltaE;
Rcorr.deltaP_vs_deltaE      = Rcorr_deltaP_vs_deltaE;
Rcorr.deltaE_vs_deltaEtrialN_R    = Rcorr_deltaE_vs_deltaEtrialN_R;
Rcorr.deltaE_vs_deltatrialN_P    = Rcorr_deltaE_vs_deltatrialN_P;
Rcorr.RPtrials_vs_trialN    = Rcorr_RPtrials_vs_trialN;

%% compute the matrix for which all correlations are minimized
bestMatrix_idx = find( min(Rcorr_global) );
bestMatrix = choice_opt.(['matrix_',num2str(bestMatrix_idx)]);
Rfields = fieldnames(Rcorr);
for iField = 1:length(Rfields)
    Rcorr_bestMatrix.(Rfields{iField}) = Rcorr.(Rfields{iField})(bestMatrix_idx);
end
end % function