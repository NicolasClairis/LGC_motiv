function [ AUC ] = AUC_gaussians_overlap( mu1, sigma1, mu2, sigma2, dispG )
% [ AUC ] = AUC_gaussians_overlap( mu1, sigma1, mu2, sigma2 )
%
% INPUTS
% mu1: mean for first distribution
%
% sigma1: sigma for second distribution
%
% mu2: mean for second distribution
%
% sigma2: sigma for second distribution
%
% dispG: display graph with the two gaussians to give an idea of the parts
% to overlap (by default set to zero if not specified)
%
% OUTPUTS
% AUC: area under the curve for the overlap
%
% See also find_gaussians_overlap.m
%
% Developed by Nicolas Clairis with great help of Quentin Feltgen - march
% 2020

%% if disp not entered, by default at zero
if ~exist('dispG','var') || isempty(dispG)
    dispG = 0;
end

%% control that function makes sense
% if any of the sigma is at zero, then it's not a normal distribution
if sigma1 == 0 || sigma2 == 0
    error('One of the sigma values is at zero, using AUC_gaussians_overlap doesn''t make any sense');
end

%% extract point where overlap
[ X1, X2 ] = find_gaussians_overlap(mu1, sigma1, mu2, sigma2 );

%% check which functions has the higher sigma to order properly
if sigma1 ~= sigma2 || mu1 ~= mu2
    if sigma1 == sigma2
        X = X1; % one single point of overlap since X1 = X2 in that case
        
        %% order distributions according to which one has the higher mean
        pdA = makedist('Normal');
        pdB = makedist('Normal');
        if mu1 < mu2
            pdA.mu      = mu1;
            pdA.sigma   = sigma1;
            pdB.mu      = mu2;
            pdB.sigma   = sigma2;
        elseif mu1 > mu2
            pdA.mu      = mu2;
            pdA.sigma   = sigma2;
            pdB.mu      = mu1;
            pdB.sigma   = sigma1;
        end
        
        %% compute AUC
        left_AUC    = cdf(pdB, X);
        right_AUC   = 1 - cdf(pdA, X);
        AUC = left_AUC + right_AUC;
        
    else
        %% order distributions according to sigma (and not the mean in this case)
        pdA = makedist('Normal');
        pdB = makedist('Normal');
        
        if sigma1 < sigma2
            pdA.mu      = mu1;
            pdA.sigma   = sigma1;
            pdB.mu      = mu2;
            pdB.sigma   = sigma2;
        elseif sigma1 > sigma2
            pdA.mu      = mu2;
            pdA.sigma   = sigma2;
            pdB.mu      = mu1;
            pdB.sigma   = sigma1;
        end
        
        %% compute the AUC
        left_AUC    = cdf(pdA, X1);
        middle_AUC  = cdf(pdB, X2) - cdf(pdB, X1); % since X2 is higher than X1
        right_AUC   = 1 - cdf(pdA, X2);
        AUC = left_AUC + middle_AUC + right_AUC;
    end
    
else
    % a normal distribution is always equal to 1
    % => if the two distributions are exactly the same => AUC = 1
    AUC = 1;
end

%% graph
if dispG == 1
    fig();
    max_val = nanmax(mu1, mu2) + 4*nanmax(sigma1, sigma2);
    x = -max_val:0.001:max_val;
    
    if sigma1 ~= sigma2 || mu1 ~= mu2
        
        % extract Y.values for each normal distribution
        pdfA_normal = pdf(pdA, x);
        pdfB_normal = pdf(pdB, x);
        
        % plot two gaussians
        plot(x, pdfA_normal, 'LineWidth',3);
        hold on;
        plot(x, pdfB_normal, 'LineWidth',3);
        
        if sigma1 == sigma2
            % split into two parts
            left_x = x <= X;
            right_x = x >= X;
            left_curve  = pdfB_normal(left_x);
            right_curve = pdfA_normal(right_x);
            
            % add area of interest
            hL = area(x(left_x), left_curve);
            hR = area(x(right_x), right_curve);
            
            hL.FaceColor = [1 0 0];
            hR.FaceColor = [0.75 0.2 0.2];
            
        else
            % split into three parts
            left_x      = x <= X1;
            middle_x    = x >= X1 & x <= X2;
            right_x     = x >= X2;
            left_curve      = pdfA_normal(left_x);
            middle_curve    = pdfB_normal(middle_x);
            right_curve     = pdfA_normal(right_x);
            
            
            % add area of interest
            hL = area(x(left_x), left_curve);
            hM = area(x(middle_x), middle_curve);
            hR = area(x(right_x), right_curve);
            
            hL.FaceColor = [1 0 0];
            hM.FaceColor = [0.8 0.2 0.2];
            hR.FaceColor = [0.6 0.3 0.3];
        end
        
    else% 1 single distribution
        pd = makedist('Normal');
        pd.mu       = mu1;
        pd.sigma    = sigma1;
        pdf_normal = pdf(pd, x);
        
        % plot area of interest
        area(x, pdf_normal);
    end
end

end % function