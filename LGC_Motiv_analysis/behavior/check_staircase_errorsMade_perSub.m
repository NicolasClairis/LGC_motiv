% script to check average proportion of errors made

% root = fullfile('C:','Users','Loco','Documents','GitHub','LGC_motiv','LGC_Motiv_results');
root = fullfile('C:','Users','clairis','Desktop','GitHub','LGC_motiv','LGC_Motiv_results');
groupFolder = 'pilots_v4_IP_Nback2_NOtaskSwitching'; % 'pilots_v3_IP_Nback2'
subFolder = [root, filesep, groupFolder, filesep];
sub = 1;
init = 'SN'; % 'DU'/'SN'/'AD'/

loadStruct = load([subFolder,'IP_pilot_data',init,'_sub_',num2str(sub),'.mat']);
n_trialsPerStaircase = 5;
n_efforts = 2;
n_repeats = 2;
n_sessions = 2;
n_trials = n_trialsPerStaircase*n_efforts*n_repeats*n_sessions;

% check average number of errors per condition
[n_errorsEffort.E1,...
    n_errorsEffort.E2,...
    n_errorsEffort.E3,...
    n_trialSuccess.E1,...
    n_trialSuccess.E2,...
    n_trialSuccess.E3] = deal( [] );
[nChoice.Em.E1, nChoice.Em.E2, nChoice.Em.E3] = deal(0);

for iEffort = 1:n_efforts
    for iSession = 1:n_sessions
        for iRepeat = 1:n_repeats
            Echosen_tmp = loadStruct.perfSummary.mental.(['repeat_nb',num2str(iRepeat)]).(['session_nb',num2str(iSession)]).(['Effort_lvl',num2str(iEffort)]).E_chosen;
            for iTrial = 1:n_trialsPerStaircase
                
                % update choice
                switch Echosen_tmp(iTrial)
                    case 1
                        nChoice.Em.E1 = nChoice.Em.E1 + 1;
                    case 2
                        nChoice.Em.E2 = nChoice.Em.E2 + 1;
                    case 3
                        nChoice.Em.E3 = nChoice.Em.E3 + 1;
                end
                
                % update proportion of errors and trial success
                switch loadStruct.perfSummary.mental.(['repeat_nb',num2str(iRepeat)]).(['session_nb',num2str(iSession)]).(['Effort_lvl',num2str(iEffort)]).perfSummary{1,iTrial}.n_max_to_reach
                    case 2
                        % check number of errors per trial
                        n_errorsEffort.E1 = [n_errorsEffort.E1, loadStruct.perfSummary.mental.(['repeat_nb',num2str(iRepeat)]).(['session_nb',num2str(iSession)]).(['Effort_lvl',num2str(iEffort)]).perfSummary{1,iTrial}.n_errorsMade];
                        % check proportion of failed trials (= more than 3 errors)
                        n_trialSuccess.E1 = [n_trialSuccess.E1, loadStruct.perfSummary.mental.(['repeat_nb',num2str(iRepeat)]).(['session_nb',num2str(iSession)]).(['Effort_lvl',num2str(iEffort)]).perfSummary{1,iTrial}.success];
                    case 4
                        % check number of errors per trial
                        n_errorsEffort.E2 = [n_errorsEffort.E2, loadStruct.perfSummary.mental.(['repeat_nb',num2str(iRepeat)]).(['session_nb',num2str(iSession)]).(['Effort_lvl',num2str(iEffort)]).perfSummary{1,iTrial}.n_errorsMade];
                        % check proportion of failed trials (= more than 3 errors)
                        n_trialSuccess.E2 = [n_trialSuccess.E2, loadStruct.perfSummary.mental.(['repeat_nb',num2str(iRepeat)]).(['session_nb',num2str(iSession)]).(['Effort_lvl',num2str(iEffort)]).perfSummary{1,iTrial}.success];
                    case 6
                        % check number of errors per trial
                        n_errorsEffort.E3 = [n_errorsEffort.E3, loadStruct.perfSummary.mental.(['repeat_nb',num2str(iRepeat)]).(['session_nb',num2str(iSession)]).(['Effort_lvl',num2str(iEffort)]).perfSummary{1,iTrial}.n_errorsMade];
                        % check proportion of failed trials (= more than 3 errors)
                        n_trialSuccess.E3 = [n_trialSuccess.E3, loadStruct.perfSummary.mental.(['repeat_nb',num2str(iRepeat)]).(['session_nb',num2str(iSession)]).(['Effort_lvl',num2str(iEffort)]).perfSummary{1,iTrial}.success];
                end
                
            end % trial
        end % repetition
    end % session
end % effort

n_AVGerrorsEffort.E1 = mean(n_errorsEffort.E1);
n_AVGerrorsEffort.E2 = mean(n_errorsEffort.E2);
n_AVGerrorsEffort.E3 = mean(n_errorsEffort.E3);
n_success.E1 = sum(n_trialSuccess.E1)/length(n_trialSuccess.E1);
n_success.E2 = sum(n_trialSuccess.E2)/length(n_trialSuccess.E2);
n_success.E3 = sum(n_trialSuccess.E3)/length(n_trialSuccess.E3);
