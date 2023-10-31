function[betas, y_var_fit] = glmfit_across_multiple_runs(run_cstt, x_var, y_var, n_runs)
% [betas, y_var_fit] = glmfit_across_multiple_runs(run_cstt, x_var, y_var, n_runs)
% glmfit_across_multiple_runs will perform glmfit (requires Statistics and
% Machine Learning Matlab Toolbox) on data containing multiple runs
% depending on which runs are to be included or not.
%
% INPUTS
% run_cstt: structure with one subfield for each run (named 'runX' for the
% run number X). Each subfield should be structured as a (nTrials*1) vector
% with zeros when it's not the current run number, ones for the current run
% number and NaN if the data has to be ignored. 
% Example: run_cstt.run1 = [1; 1; 1; 0; 0; 0; NaN; NaN; NaN] means that the
% 3 first trials are run 1, the 3 next are run 2 and the three last have to
% be removed. Accordingly, run_cstt.run3 = [0; 0; 0; 0; 0; 0; NaN; NaN;
% NaN]; meaning that sum(run_cstt.run3)=0 and therefore it will be
% automatically removed from the analysis.
%
% x_var: (nTrials*nRegressors) matrix containing all the regressors to be 
% included in the GLM but the run constant.
%
% y_var: (nTrials*1) vector containing the data to fit.
%
% n_runs: number indicating the total number of runs to check (ie how many
% subfields should be contained inside run_cstt)
%
% OUTPUTS
% betas: (nRegressors*1) vector containing the beta for each regressor
% starting with the run constant and then with the regressors inside x_var.
% Note that if a run number X is ignored, it will provide a beta(X) = NaN
%
% y_var_fit: corresponding fit for y_var.
%
% Designed by N. Clairis - 31/10/2023

%% check how many runs are ok
n_runs_ok = 0;
id_runs_ok = zeros(n_runs, 1);
for iR = 1:n_runs
    if sum(run_cstt.(['run',num2str(iR)]),'omitnan') > 0
        n_runs_ok = n_runs_ok + 1;
        id_runs_ok(iR) = 1;
    end
end % run loop

%% perform glmfit on correct sessions
x_var_with_cstt = [];
for iR = 1:n_runs
    if id_runs_ok(iR) == 1
        x_var_with_cstt = [x_var_with_cstt, run_cstt.(['run',num2str(iR)])];
    end
end
x_var_with_cstt = [x_var_with_cstt, x_var];
betas_tmp = glmfit(x_var_with_cstt, y_var, 'normal','Constant','Off');
y_var_fit = glmval(betas_tmp, x_var_with_cstt, 'identity','Constant','Off');

%% distribute data depending on which runs were removed
if n_runs_ok == n_runs
    betas = betas_tmp;
elseif n_runs_ok < n_runs
    betas = NaN(n_runs+size(x_var,2), 1);
    jR = 0;
    for iR = 1:n_runs
        if id_runs_ok(iR) == 1
            jR = jR + 1;
            betas(iR) = betas_tmp(jR);
        end
    end % loop through runs
    % upload the other regressors
    betas((n_runs+1):end) = betas_tmp((jR+1):end);
end % distribute betas differently in case where some runs were not included

end % function