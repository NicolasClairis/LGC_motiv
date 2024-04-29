function [gof_measures]=goftests_gls(pars,minval)

input_data=evalin('base','data');
N=size(input_data,1);
vars=size(input_data,2);
dof=(vars*(vars+1))/2-length(pars);
dof_base=(vars*(vars+1))/2-vars;

% calculating chi-square value for appropriate baseline model

[pars_base,minval_base]=anneal(@model_est_baseline_gls,ones(1,vars));
T_base=(N-1)*minval_base;

% calculating chi-square value for model

T=(N-1)*minval;
CHISQ=T;
p=1-chi2cd(T,dof);

display(CHISQ);
display(p);

% calculating NFI

NFI=(T_base-T)/T_base;
display(NFI);

% calculating NNFI

NNFI=(((T_base/dof_base)-(T/dof))/((T_base/dof_base)-1));
display(NNFI);

% calculating AIC
k=length(pars);
AIC=T+2*k;
display(AIC);

%calculating SABIC
SABIC=T+log((N+2)/24)*k;
display(SABIC);

% calculating RMSEA

RMSEA=sqrt(T-dof)/sqrt(dof*(N-1));

if RMSEA<0
    RMSEA=0;
end
       
display(RMSEA);

    % calculate confidence intervals
    
    upp_rmsea=RMSEA;
    p_crit=ncx2cd(T,dof,RMSEA);
    
    while p_crit>0.05
        
       upp_rmsea=upp_rmsea+0.0001; 
       p_crit=ncx2cd(T,dof,(N-1)*dof*upp_rmsea*upp_rmsea);  
        
    end
    
    low_rmsea=RMSEA;
    p_crit=ncx2cd(T,dof,RMSEA);
    
    while p_crit>0.05
      
     low_rmsea=low_rmsea-0.0001;    
     p_crit=ncx2cd(T,dof,(N-1)*dof*low_rmsea*low_rmsea);  
        
    end
    
    if low_rmsea<0
        low_rmsea=0; 
    end

gof_measures=struct('CHISQ',{T,p,NaN},'NFI',{NFI,NaN,NaN},'NNFI',{NNFI,NaN,NaN},'AIC',{AIC,NaN,NaN},'SABIC',{SABIC,NaN,NaN},'RMSEA',{RMSEA,low_rmsea,upp_rmsea});

end