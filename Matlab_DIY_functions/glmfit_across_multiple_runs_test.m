n_trials = 40;
n_runs = 2;

%% testing script
run_cstt.run1 = [ones(n_trials/2,1);zeros(n_trials/2,1)];
run_cstt.run2 = [zeros(n_trials/2,1);ones(n_trials/2,1)];
x_var = (1:n_trials)'+randn(n_trials,1).*10;
y_var = (run_cstt.run1.*35) + (x_var(:,1).*8) + ((randn(n_trials,1).*12).*15);
[betas, y_var_fit] = glmfit_across_multiple_runs(run_cstt, x_var, y_var, n_runs);
[x_var_sorted, x_var_sort_idx] = sort(x_var);
y_var_fit_sorted = y_var_fit(x_var_sort_idx);

figure;
scatter(x_var, y_var);
hold on;
plot(x_var_sorted, y_var_fit_sorted,'LineStyle','-','Color',[143 143 143]./255);

%% testing multiple X variables
x_var = [(1:n_trials)'+randn(n_trials,1).*10,...
    randn(n_trials,1).*10];
y_var = (run_cstt.run1.*35) + (x_var(:,1).*8)+ (x_var(:,2).*4) + ((randn(n_trials,1).*12).*15);
[betas, y_var_fit] = glmfit_across_multiple_runs(run_cstt, x_var, y_var, n_runs);

%% testing to remove one run
run_cstt.run1 = [ones(n_trials/2,1);NaN(n_trials/2,1)];
run_cstt.run2 = [zeros(n_trials/2,1);NaN(n_trials/2,1)];
x_var = (1:n_trials)'+randn(n_trials,1).*10;
y_var = (run_cstt.run1.*35) + (x_var(:,1).*8) + ((randn(n_trials,1).*12).*15);
[betas, y_var_fit] = glmfit_across_multiple_runs(run_cstt, x_var, y_var, n_runs);
[x_var_sorted, x_var_sort_idx] = sort(x_var);
y_var_fit_sorted = y_var_fit(x_var_sort_idx);

figure;
scatter(x_var, y_var);
hold on;
plot(x_var_sorted, y_var_fit_sorted,'LineStyle','-','Color',[143 143 143]./255);
