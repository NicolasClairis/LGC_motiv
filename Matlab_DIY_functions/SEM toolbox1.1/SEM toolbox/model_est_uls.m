function [ULS]=model_est_uls(x)

% obtaining asymmetric, symmetric and filter matrices

model=evalin('base','model_specs');

% calculate F matrix

F=zeros(model.obs,model.obs+model.lat);

for i=1:model.obs,
       
    F(i,model.lat+i)=1;
    
end

% generating symmetric matrix

par_num=0;
dimension=size(F,2);
S=zeros(dimension,dimension);

length_sym=size(model.sym,1);

for i=1:length_sym
   
    row=model.sym(i,1);
    col=model.sym(i,2);
    
    par_num=par_num+1;
    S(row,col)=x(1,par_num);
       
    if row~=col 
        S(col,row)=x(1,par_num);
    end
    
end

% fixing variance of latent variables to 1

%for i=1:model.lat,
   
%    S(i,i)=1;
    
%end

% generating asymmetric matrix

dimension=size(F,2);
A=zeros(dimension,dimension);

length_asym=size(model.asym,1);

% getting indices of factor loadings which are set to 1

    cul_indices=[];

for cul_ind=1:model.lat,

    tek1=find(model.asym(:,2)==cul_ind);
    tek2= model.asym(tek1,1)>model.lat;
    tek1=tek1(tek2,1);
    cul_indices=vertcat(cul_indices,tek1(1,1));
    tek1=[];
    tek2=[];

end

for i=1:length_asym
    
   row=model.asym(i,1);
   col=model.asym(i,2); 
   
if length(intersect(i,cul_indices))==1
    
      A(row,col)=1;
    
else
    
      par_num=par_num+1;
      A(row,col)=x(1,par_num);

end
    
end

% obtaining input data

input_data=evalin('base','bt_data');

% checking if number of observations is sufficient

%num_of_obs=size(input_data,1);

%if (num_of_obs/num_pars)<10
   
%    warning('SEM:datapoints','Number of samples might be too few for specified number of free parameters. Either increase the number of samples or decrease the number of free parameters');
    
%end

% estimating model covariance matrix

I=eye(size(A));

model_cov_mat = F * inv(I-A) * S * (inv(I-A))' * F';  

% calculating data covariance matrix

obs_cov_mat=cov(input_data);

% calculating ULS value

p = size(F,1);

ULS = abs(0.5*trace((obs_cov_mat-model_cov_mat)^2));

end