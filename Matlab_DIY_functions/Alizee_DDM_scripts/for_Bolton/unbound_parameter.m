function not_bound=unbound_parameter(bounded_parameter,lim1,lim2)

not_bound=-log(bounded_parameter-lim2)+log(lim1-bounded_parameter);
% a=lim1;
% b=lim2-lim1;
% not_bound=-log((a*b)/(bounded_parameter-a));
localEps = 1e-6;

% not_bound(bound_parameter(not_bound,lim1,lim2)>=lim2-localEps)=lim2;
% not_bound(bound_parameter(not_bound,lim1,lim2)<=lim1+localEps)=lim1;
% 
 not_bound(real(not_bound)==-Inf)=-10;
 not_bound(real(not_bound)==Inf)=10;
