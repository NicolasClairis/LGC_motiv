function [ci_bt,se_bt]=model_ci_bt(pars,data)

% initialise

replications=1000;
num_subs=size(data,1);

 pars_dist=zeros(replications,length(pars));

% populating dataset of parameter distributions

for i=1:replications,
    
   rand_indices=randi(num_subs,num_subs,1);
   bt_data=data(rand_indices,:);
   
   assignin('base','bt_data',bt_data);
   
   [pars_bt,minval_bt]=anneal(@model_est_ml_bt,pars);
   pars_dist(i,:)=pars_bt;
    
   display(i);
   
end

% estimating standard errors

se_bt=std(pars_dist,1,1);

% estimating confidence intervals

ci_bt=zeros(2,length(pars));

% estimating lower bound

for lb_par=1:length(pars),
   
    dist_dataset=[];
    dist_dataset=pars_dist(:,lb_par);
    lb_thresh=min(dist_dataset)-0.01;
    
    while (length(find(dist_dataset(:,1)>lb_thresh))/replications)>0.975
        
    lb_thresh=lb_thresh+0.001;    
        
    end
    
    lb_thresh=lb_thresh-0.001;
    ci_bt(1,lb_par)=lb_thresh;
    
    display(lb_par);
    
end

% estimating upper bound

for ub_par=1:length(pars),
   
    dist_dataset=[];
    dist_dataset=pars_dist(:,ub_par);
    ub_thresh=max(dist_dataset)+0.01;
    
    while (length(find(dist_dataset(:,1)<ub_thresh))/replications)>0.975
        
    ub_thresh=ub_thresh-0.001;    
        
    end
    
    ub_thresh=ub_thresh+0.001;
    ci_bt(2,ub_par)=ub_thresh;
    
    display(ub_par);
    
end

end