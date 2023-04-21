function bounded_parameter=bound_parameter(parameter,lim1,lim2)
for i=1:length(parameter)
%      bounded_parameter(i)=lim1+(lim2-lim1)/(1+exp(-parameter(i))); 
% bound between a and b+a
a=lim1;
b=lim2-lim1;
bounded_parameter(i)=a+b*sgm(parameter(i));
end
