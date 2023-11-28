function [raw_res,std_res]=residuals(pars)

% obtaining input data
input_data=evalin('base','data');
obs_cov_mat=cov(input_data);

% calculating new model covariance matrix

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
S(sidxs)=pars(1,1:snum);
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
A(aidxs)=pars(1,snum+1:length(pars));

I=eye(size(A));
model_cov_mat = F * inv(I-A) * S * (inv(I-A))' * F';

raw_res=obs_cov_mat-model_cov_mat;

% calculating mean
[n,m]=size(raw_res);
elements=n*m;
total=0;

for i=1:n,
    
    for j=1:m,
              
        total=total+raw_res(i,j);
        
    end
    
end

avg=(total/elements);

% calculating standard deviation
agg=0;

for i=1:n,
    
    for j=1:m,
   
    agg=agg+(raw_res(i,j)-avg)^2;     
        
    end
    
end

sd=sqrt(agg/(elements-1));

%calculating standardised residuals
for i=1:n,
    
    for j=1:m,
        
    std_res(i,j)=(raw_res(i,j)-avg)/sd;
        
    end
    
end

end