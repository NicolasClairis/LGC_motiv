function [chi,p]=lratiotest(minval,orig_num,new_num)

% getting original ML value
orig_lhood=minval;

% getting new ML value
[new_pars,new_minval]=anneal(@model_est_revised,ones(1,new_num));
new_lhood=new_minval;

% calculating test statistic
dof=abs(orig_num-new_num);

if orig_num<new_num   
chi=-2*log(new_lhood/orig_lhood);
elseif orig_num>new_num      
chi=-2*log(orig_lhood/new_lhood);
else
error('Invalid model comparison: Identical or non-nested models are being compared');
end

% calculating associated p-value
p=1-chi2cd(chi,dof);

end