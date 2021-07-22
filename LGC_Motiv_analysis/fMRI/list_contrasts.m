function[list_con, con_vec, n_con] = list_contrasts(GLM, nb_runs)

GLMprm = which_GLM(GLM);

% extract number of regressors
n_con = 0;

n_con = n_con + 1;
list_con{n_con} = 'cross';
con_vec(n_con,:) = 

end % function