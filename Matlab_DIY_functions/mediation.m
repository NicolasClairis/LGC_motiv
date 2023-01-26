function[a,b,c,c_prime, pval] = mediation(X, M, Y, X_nm, M_nm, Y_nm, dispResults)
% [a,b,c,c_prime pval] = mediation(X, M, Y, X_nm, M_nm, Y_nm, dispResults)
% mediation will perform a mediation going from X to Y through M as a
% mediator and will display the corresponding betas (rounded at 3 values
% after the coma) and p.values.
%
% Note: Ideally, vectors should be entered in n*1 format, but the script
% will automatically flip each one of them if they are entered in 1*n
% format.
%
% INPUTS
% X: n*1 vector with variable X
%
% M: n*1 vector with variable M (mediator)
%
% Y: n*1 vector with variable Y
%
% X_nm: name of X for text output (='X' if left empty)
%
% M_nm: name of X for text output (='M' if left empty)
%
% Y_nm: name of X for text output (='Y' if left empty)
%
% dispResults: display results of betas and p.value?
% by default, it will display the results (1) but set to 0 if you don't
% want it to avoid spamming output in command window
%
% OUTPUTS
% a: beta of X=>M path
%
% b: beta of M=>Y path taking X into account
%
% c: beta of X=>Y path
%
% c_prime: beta of X=>Y path taking M into account
%
% pval: structure with p.value for each path
%
% Developped by Nicolas Clairis - 17/08/2022 under Jules Brochard advice

%% by default display the results
if ~exist('dispResults','var') || isempty(dispResults)
    dispResults = 1;
end

%% check names for each member of the mediation
if ~exist('X_nm','var') || isempty(X_nm)
    X_nm = 'X';
end
if ~exist('M_nm','var') || isempty(M_nm)
    M_nm = 'M';
end
if ~exist('Y_nm','var') || isempty(Y_nm)
    Y_nm = 'Y';
end

%% verify that the data are correctly flipped
% check X
if size(X,1) == 1 && size(X,2) > 1
    X = X';
elseif size(X,1) > 1 && size(X,2) == 1
    % nothing to perform in that case
else
    error('problem with X');
end

% check M
if size(M,1) == 1 && size(M,2) > 1
    M = M';
elseif size(M,1) > 1 && size(M,2) == 1
    % nothing to perform in that case
else
    error('problem with M');
end

% check Y
if size(Y,1) == 1 && size(Y,2) > 1
    Y = Y';
elseif size(Y,1) > 1 && size(Y,2) == 1
    % nothing to perform in that case
else
    error('problem with Y');
end

%% perform each path of the mediation

% test correlation between X and M (path a)
[betas_1,~,stats_1] = glmfit(X, M,'normal');
a = betas_1(2);
pval.a = stats_1.p(2);

% test correlation between X, M and Y (path b)
[betas_2,~,stats_2] = glmfit([M, X], Y,'normal');
b       = betas_2(2);
c_prime = betas_2(3);
pval.b          = stats_2.p(2);
pval.c_prime    = stats_2.p(3);

% test also direct path between X and Y (path c)
[betas_3,~,stats_3] = glmfit(X, Y,'normal');
c = betas_3(2);
pval.c = stats_3.p(2);

%% display relevant p.values in the command window to summarize the 
% results of the mediation
if dispResults == 1
    
    %% show results in matlab dialog
    roundingVal = 3;
    
    % path a
    disp([X_nm,' -> ',M_nm,' (path a): ',...
        M_nm,' = ',num2str(round(betas_1(1),roundingVal)),...
        ' + (',num2str(round(a, roundingVal)),')*',X_nm,';']);
    disp(['p(a=',X_nm,') = ',num2str(pval.a)]);
    disp(' ');
    % path b
    disp([M_nm,' -> ',Y_nm,' (path b): ',...
        Y_nm,' = ',num2str(round(betas_2(1), roundingVal)),...
        ' + (',num2str(round(b, roundingVal)),')*',M_nm,...
        ' + (',num2str(round(c_prime,roundingVal)),')*',X_nm,';']);
    disp(['p(b=',M_nm,') = ',num2str(pval.b),';']);
    disp(['p(c''=',X_nm,') = ',num2str(pval.c_prime)]);
    disp(' ');
    % path c
    disp([X_nm,' -> ',Y_nm,' (path c): ',...
        Y_nm,' = ',num2str(round(betas_3(1), roundingVal)),...
        ' + (',num2str(round(c, roundingVal)),')*',X_nm,';']);
    disp(['p(c=',X_nm,') = ',num2str(pval.c)]);
    disp(' ');
    
    %% display results with a figure
    pSize = 40;
    lWidth = 3;
    black = [0 0 0];
    grey = [143 143 143]./255;
    
    % extract relevant data for fit
    X_ascOrder = sort(X);
    % X => Y (direct path c)
    Y_c_fit = glmval(betas_3, X_ascOrder,'identity');
    % X => Y after removing M (direct path c')
    Y_res_without_M = Y - b.*M;
    Y_cPrime_fit = betas_2(1) + c_prime.*X_ascOrder;
    % X => M (path a)
    M_fit = glmval(betas_1, X_ascOrder,'identity');
    % M => Y removing any influence of X (path b)
    M_res_without_X = M - a.*X;
    M_res_ascOrder = sort(M_res_without_X);
    Y_res_without_X = Y - c_prime.*X;
    Y_b_fit = betas_2(1) + b.*M_res_ascOrder;
    
    fig;
    % X => Y direct path without mediation (c)
    subplot(2,2,1);
    scat_hdl = scatter(X, Y);
    scat_hdl.LineWidth = lWidth;
    scat_hdl.MarkerEdgeColor = black;
    hold on;
    fit_hdl = plot(X_ascOrder, Y_c_fit);
    fit_hdl.LineWidth = lWidth;
    fit_hdl.Color = grey;
    xlabel([X_nm,' - direct path (c)']);
    ylabel(Y_nm);
    legend_size(pSize);
    
    % X => Y direct path competition with M (c')
    subplot(2,2,2);
    scat_hdl = scatter(X, Y_res_without_M);
    scat_hdl.LineWidth = lWidth;
    scat_hdl.MarkerEdgeColor = black;
    hold on;
    fit_hdl = plot(X_ascOrder, Y_cPrime_fit);
    fit_hdl.LineWidth = lWidth;
    fit_hdl.Color = grey;
    xlabel([X_nm,' - direct path (c'')']);
    ylabel(Y_nm);
    legend_size(pSize);
    
    % X => M path
    subplot(2,2,3);
    scat_hdl = scatter(X, M);
    scat_hdl.LineWidth = lWidth;
    scat_hdl.MarkerEdgeColor = black;
    hold on;
    fit_hdl = plot(X_ascOrder, M_fit);
    fit_hdl.LineWidth = lWidth;
    fit_hdl.Color = grey;
    xlabel(X_nm);
    ylabel(M_nm);
    legend_size(pSize);
    
    % M => Y path (after removing X => M)
    subplot(2,2,4);
    scat_hdl = scatter(M_res_without_X, Y_res_without_X);
    scat_hdl.LineWidth = lWidth;
    scat_hdl.MarkerEdgeColor = black;
    hold on;
    fit_hdl = plot(M_res_ascOrder, Y_b_fit);
    fit_hdl.LineWidth = lWidth;
    fit_hdl.Color = grey;
    xlabel(M_nm);
    ylabel(Y_nm);
    legend_size(pSize);
end

end % function