function [pars,se_bt,ci_bt,mn,id,gof,res]=sem_fit_bt(model_specs,data) 

% model_specs should be a structure with 4 fields - model_specs.asym_data, 
% model_specs.sym_data, model_specs.obs and model_specs.lat. the
% model_specs.sym_data and model_specs.asym_data fields contain symmetric and asymmetric
% matrices, with 'inf' indicating a free parameter. parameters can also be fixed to a specified 
% value by entering that value in the corresponding matrix element. a fixed value of 0 means no relationship
% between those two variables. The model_specs.obs and model_specs.lat fields simply specify 
% the number of observed and latent variables respectively. The model_specs 
% structure for SEM model 5.11 in file 'SEM models.pdf' is given in SEM toolbox 
% directory. A convenient way of generating the model specification would be 
% to create it using the GUI and save it.

% 'data' should be an MxN matrix, with M=number of observations, N=number of
% observed variables

% the outputs are in following order: vector containing parameter values -
% order is columnwise, first from symmetric matrix and then from
% asymmetric. next outputs are vector containing boot-strapped standard errors, matrix
% containing bootstrapped confidence intervals, flags for multivariate normality
% assumption (1 if violated, 0 if satisfied), model identification (1 if
% not identified, 0 if identified), structure containing fitness measures
% ,including p-value for chi-sq. (2nd row) and 90% CIs for RMSEA (2nd row 
% and 3rd row) and matrix of residuals.


% copying to base workspace

assignin('base','model_specs',model_specs);
assignin('base','data',data);

% checking for multivariate normality

alpha=0.05;
[p1,p1c,p2]=multskewkurt(data,alpha);

% calculates number of free parameters
symel=0;
for idx1=1:size(model_specs.sym_data,1),
    
    for idx2=1:idx1,
        
        if model_specs.sym_data(idx1,idx2)==Inf,
   
        symel=symel+1;
        
        end
        
    end
    
end

par_num=symel+length(find(model_specs.asym_data==Inf))-model_specs.lat;
sample_ratio=(size(data,1)/par_num);

mn=0;

if p2 <=alpha
  mn=1;
  warning('Assumption of multivariate normality is violated (Kurtosis)');    
end

if (sample_ratio >= 10)
    if  p1 <=alpha
        mn=1;
        warning('Assumption of multivariate normality is violated (Skewness)');
    end
else
    if p1c <=alpha
        mn=1;
        warning('Assumption of multivariate normality is violated (Skewness)');   
    end 
end

% estimating parameters using ML

error_flag_id=0;

% calculates number of free parameters

symel=0;
for idx1=1:size(model_specs.sym_data,1),
    
    for idx2=1:idx1,
        
        if model_specs.sym_data(idx1,idx2)==Inf,
   
        symel=symel+1;
        
        end
        
    end
    
end

par_num=symel+length(find(model_specs.asym_data==Inf))-model_specs.lat;
val = ones(1,par_num);

% checking for possibility of identification

num_nr=(model_specs.obs*(model_specs.obs+1))/2;
if num_nr<par_num
    
    errordlg('Model cannot be identified - number of unknowns (parameters) exceeds number of knowns (non-redundant elements)','Invalid model','modal');
    error_flag_id=1;
    
end

% actual parameter estimation

if error_flag_id==0  
     
[pars,minval]=anneal(@model_est_ml,val);

display('Parameters have been estimated');

samp_ratio=size(data,1)/par_num;

if samp_ratio<10   
    warning('Number of samples is likely too few for given number of parameters');
end

% estimating standard errors and confidence intervals

[ci_bt,se_bt]=model_ci_bt(pars,data);

display('Standard errors and CIs have been estimated');

num_pars=length(pars);
  
[H]=hessian(@model_est_ml,pars);   

% estimating the Information matrix
I=-1*H;

% calculating rank of Information matrix
r=rank(I);

% display appropriate message
if r~=num_pars,
    errordlg('Model is not identified','Model ID','modal');
    id=1;
else  
    id=0;
end

if id==0
% Goodness of fit measures
  
gof=goftests_ml(pars,minval); 

display('Goodness-of-fit measures have been estimated');

% calculate residuals

[raw_res,res]=residuals(pars);

end

end

end