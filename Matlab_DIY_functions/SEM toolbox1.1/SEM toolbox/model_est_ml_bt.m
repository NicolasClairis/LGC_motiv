function [ML]=model_est_ml_bt(x)

% obtaining asymmetric, symmetric and filter matrices

model=evalin('base','model_specs');

% calculate F matrix

F=zeros(model.obs,model.obs+model.lat);

for i=1:model.obs,
       
    F(i,model.lat+i)=1;
    
end

% generating symmetric matrix
sidxs=[];
S=zeros(size(model.sym_data));
S=model.sym_data;
[sidxs(:,1),sidxs(:,2)]=find(model.sym_data==Inf);
sidxsind=find((sidxs(:,2)<=sidxs(:,1))==1);
sidxs=sidxs(sidxsind,:);
snum=length(sidxs);
sidxs=sub2ind(size(S),sidxs(:,1),sidxs(:,2));
S(sidxs)=x(1,1:snum);
for idx1=1:size(S,1),
    for idx2=1:idx1,
        if idx1~=idx2,
        S(idx2,idx1)=S(idx1,idx2);
        end
    end
end

% generating asymmetric matrix

A=zeros(size(model.asym_data));
A=model.asym_data;

for cidx=1:model.lat,

    tek1=find(A(:,cidx)==Inf);
    tek2= find(tek1>model.lat);
    tek1=tek1(tek2,1);
    A(tek1(1,1),cidx)=1;

end

[aidxs]=find(A==Inf);
A(aidxs)=x(1,snum+1:length(x));

% checking if model is identified

% num_ob=size(F,2);
% num_nr=(num_ob*(num_ob+1))/2;
% num_pars=length(x);

% if num_nr<num_pars
    
%    error('Model is unidentified');
    
% end 

% obtaining input data

input_data=evalin('base','bt_data');

% estimating model covariance matrix

I=eye(size(A));

model_cov_mat = F * inv(I-A) * S * (inv(I-A))' * F';  

% calculating data covariance matrix

obs_cov_mat=cov(input_data);

% calculating maximum likelihood value

p = size(F,1);

ML = log(det(model_cov_mat))-log(det(obs_cov_mat)) + trace(obs_cov_mat*inv(model_cov_mat)) - p;

end