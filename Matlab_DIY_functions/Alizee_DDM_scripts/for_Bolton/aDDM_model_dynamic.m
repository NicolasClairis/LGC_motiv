function [gx]=aDDM_model_dynamic(x,param,uu,xx,fig)
% Dynamic version of the attentional drift diffusion model
% if fig==1, plot a figure of each random walk following imposed fixations.


drift   = param(1); % drift rate (units/sec) 0.1
p.start = param(2); % starting point 
p.s     = param(3); % standard deviation of drift (units/sec)
theta   = param(4); % theta parameter from Krajbich 2010 (weight put on the fixated item compared to the unfixated item. 0.3 in krajbich studies)

drift_borns = 1; % boundaries

u         = uu{1};
DV        = u(1,:)-u(2,:);
fixations = uu{2};

%% Define variables for the random walk:
p.dt        = .016; %step size for simulations (seconds)
nReps = 1;  %number of simulated walks

%% Simulating

for trial=1:size(fixations,2)
    
    %inialize some parameters:
    alive    = true(length(u(1,trial)),nReps);  % index vector for walks that haven't terminated
    y        = p.start*ones(length(u(1,trial)),nReps); % Starting, and to be current position for each walk
    response = zeros(length(u(1,trial)),nReps);  % will be filled with -1 or 1 for incorrect and correct
    RT       = zeros(length(u(1,trial)),nReps);  % will be filled with RT values in seconds
    dy       = zeros(nReps,length(u(1,trial)));
    tStep    = zeros(length(u(1,trial)),1);
    
    for f=1:length(fixations(:,trial))
        % while sum(alive) %loop until all walks are terminated
        tStep = tStep+1;
        
        % SPECIFIC to aDDM
        if     fixations(f,trial) == 0 % no option fixated
            p.u=0;
            
        elseif fixations(f,trial) == 1 % fixation for left
            p.u=(drift*(u(1,trial)-theta*u(2,trial)));
            
        elseif fixations(f,trial) == -1 % fixation for right
            p.u=-(drift*(u(2,trial)-theta*u(1,trial)));
            
        elseif isnan(fixations(f,trial)) % no recording of fixation
            p.u=DRIFT(find(DRIFT(:,trial)~=0,1,'last'),trial);
            if isempty(p.u);
                p.u=0;
            end
        end
        
        %for the 'randn' implementation, use this line:
        dy(alive,1) = p.u'*p.dt +p.s*sqrt(p.dt)*randn(1,sum(alive))';
        
        %increment the 'living' walks
        y(alive) = y(alive)+ dy(alive);
        
        %find the walks that reached the boundary
        aboveA = find(y>=drift_borns & alive);
        if ~isempty(aboveA)
            
            response(aboveA) = 0.5.*ones(length(aboveA),1)+ones(length(aboveA),1).*(p.u(aboveA))'./(drift_borns*2); %correct response
            
            alive(aboveA)    = false; %'kill' the walk
            RT(aboveA)       = tStep(aboveA)*p.dt; %record the RT
        end
        
        %find the walks that reached the '-b' boundary
        belowB = find(y<=-drift_borns & alive);
        if ~isempty(belowB)
            
            response(belowB) = 0.5*ones(length(belowB),1) +(ones(length(belowB),1).*(p.u(belowB))'./(drift_borns*2));  %incorrect response
            
            alive(belowB)    = false; %'kill' the walk
            RT(belowB)       = tStep(belowB)*p.dt; %record the RT
        end
        
        y_tracked(f,trial) = y;
        DRIFT(f,trial)     = p.u;
        
    end
    
    if fig==1
        figure
        
        for fi=1:length(fixations(:,trial))
            if fixations(fi,trial)==1
                h=plot([fi fi], [-drift_borns*fixations(fi,trial) drift_borns*fixations(fi,trial)],'color',[1 0.8 0.8],'linewidth',5);
            elseif fixations(fi,trial)==-1
                h=plot([fi fi], [-drift_borns*fixations(fi,trial) drift_borns*fixations(fi,trial)],'color',[0.8 0.8 1]','linewidth',5);
            end
            hold on
        end
        
        xlim([0 length(fixations(:,trial))])
        hold on
        plot(y_tracked(:,trial),'k','linewidth',2)
        hold on
        plot([0 f],[drift_borns drift_borns],'k')
        hold on
        plot([0 f],[-drift_borns -drift_borns],'k')
        hold on
        plot([0 f],[0 0],'k')
        hold on
        title(['DV=' num2str(DV(trial)) ' , V1=' num2str(u(1,trial)) ' , V2=' num2str(u(2,trial))])
        ylim([-drift_borns-0.01 drift_borns+0.01])
        xlim([0 300])
        
    end
    
    if ~isempty(find(abs(y_tracked(:,trial))>=1,1,'first'))
        gx(2,trial)=find(abs(y_tracked(:,trial))>=1,1,'first')*0.016;
    else
        gx(2,trial)=NaN;
    end
    
    if y_tracked(f,trial)>=0
        gx(1,trial)=1;
    else
        gx(1,trial)=0;
    end
end