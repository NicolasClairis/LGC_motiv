function [chi2diff,p]=chi2difftest(minval,orig_num,new_num)

input_data=data_gen();
N=size(input_data,1);

% getting original ML value
chisq_old=(N-1)*minval;

% getting new ML value
[new_pars,new_minval]=anneal(@model_est_revised,ones(1,new_num));
chisq_new=(N-1)*new_minval;

dof=abs(orig_num-new_num);

% calculating test statistic
if orig_num<new_num
chi2diff=chisq_old-chisq_new;
elseif orig_num>new_num
chi2diff=chisq_new-chisq_old;
else
error('Invalid model comparison: Identical or non-nested models are being compared');
end

% calculating associated p-value
p=1-chi2cd(chi2diff,dof);

end