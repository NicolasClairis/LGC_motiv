function [ML]=model_est_baseline_ml(x)

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

for i=model.lat+1:dimension,
    
    par_num=par_num+1;
    S(i,i)=x(1,par_num);
    
end

% generating asymmetric matrix

dimension=size(F,2);
A=zeros(dimension,dimension);

% obtaining input data

input_data=evalin('base','data');

% estimating model covariance matrix

I=eye(size(A));

model_cov_mat = F * inv(I-A) * S * (inv(I-A))' * F';  

% calculating data covariance matrix

obs_cov_mat=cov(input_data);

% calculating maximum likelihood value

p = size(F,1);

ML = log(det(model_cov_mat))-log(det(obs_cov_mat)) + trace(obs_cov_mat*inv(model_cov_mat)) - p;

end