function [conf_intvls,std_errors]=model_ci(pars,minval)

l=length(pars);
conf_intvls=zeros(2,l);
std_errors=zeros(1,l);
new_pars=zeros(1,l);

% obtaining input data

input_data=evalin('base','data');
obs_cov_mat=cov(input_data);

% calculating original chi-square value

T=(size(input_data,1)-1)*minval;

for par_index=1:l,
   
    % initialising parameters
    new_pars(1,:)=pars;
    diff=0;
      
    while diff<3.8416
        
        % incrementing parameter value    
        new_pars(1,par_index)=new_pars(1,par_index)+0.0001;
        
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
    S(sidxs)=new_pars(1,1:snum);
    
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
    A(aidxs)=new_pars(1,snum+1:length(new_pars));
        
    I=eye(size(A));
    model_cov_mat = F * inv(I-A) * S * (inv(I-A))' * F';  
    p = size(F,1);
        
    % calculating new ML value
    new_minval = log(det(model_cov_mat))-log(det(obs_cov_mat)) + trace(obs_cov_mat*inv(model_cov_mat)) - p;
        
    % calculating new chi-square value
    T_new=(size(input_data,1)-1)*new_minval;
        
    % calculating chi-square value
    diff=T_new-T;
        
    end
        
    % calculating confidence intervals
    conf_intvls(1,par_index)=pars(1,par_index)+(new_pars(1,par_index)-pars(1,par_index));
    conf_intvls(2,par_index)=pars(1,par_index)-(new_pars(1,par_index)-pars(1,par_index));
    
    std_errors(1,par_index)=(new_pars(1,par_index)-pars(1,par_index))/1.96;
    
end

end