function[u] = convert_pval_to_u(p, dof, type)
% Successfully created by Jules Brochard in 2 seconds on the 29/11/2017 <3

if type == 't'
    u = icdf('t',p,dof);
elseif type == 'f'
    
    u = icdf('F',p,1,dof);
    warning('Values only work for a F test with a single regressor as design matrix')
    % error('not ready yet. Jules where are you? Snif');
end

end