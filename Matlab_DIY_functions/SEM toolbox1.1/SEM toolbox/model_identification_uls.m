function model_identification_uls(pars)

% obtaining asymmetric, symmetric and filter matrices

[F,A,S]=model_spec(pars);

% checking for necessary condition

num_ob=size(F,2);
num_nr=(num_ob*(num_ob+1))/2;
num_pars=length(pars);

if num_nr<num_pars  
    error('Model is not identified - number of unknown variables exceed number of known variables');    
end 

% estimating the Hessian matrix

[H]=hessiancsd(@model_est_uls,pars);

% estimating the Information matrix

I=-1*H;

% calculating rank of Information matrix

r=rank(I);

% display appropriate messages

if r~=num_pars,
   
    error('Model is not identified - please carry out appropriate diagnostics');
    
else
    
    display('Model is identified!');
    
end

end