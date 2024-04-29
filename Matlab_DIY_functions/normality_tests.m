function[h, pval] = normality_tests(x, dispFig)
% [h, pval] = normality_tests(x, dispFig)
% normality_tests will perform a whole bunch of normality tests to check
% if x can be considered as normal or not. Moreover, the function will
% display a message telling you whether it is the case or not + it will
% display the distribution of x and compare it to the expected normal
% distribution.
%
% INPUTS
% x: vector of data to test
%
% dispFig:
% (0) avoid any display (no graph no message)
% (1) display figure and message regarding data normality
% (2) display only message regarding data normality but no figure
%
% OUTPUTS
% h: structure indicating, for each test, whether null hypothesis (ie that
% data is normal) can be accepted (should be equal to 0 then) or if it 
% failed (equal to 1 in that case).
%
% pval: structure indicating, for each test, the corresponding p.value for
% the test between null hypothesis (ie that data is normal) and alternative
% (data is not normal). If pval.(test) <0.05 then it means that the data is
% not distributed normally with a significance of pval.(test).

%% display data by default
if ~exist('dispFig','var') || ~isempty(dispFig) || ~ismember(dispFig,[0,1,2])
    dispFig = 1;
end

%% Anderson-Darling test
[h.adt,pval.adt] = adtest(x);

%% Jarque-Bera
[h.jbt,pval.jbt] = jbtest(x);

%% Lilliefors normality test
[h.llt,pval.llt] = lillietest(x);

%% One-sample Kolmogorov-Smirnov test
[h.kst,pval.kst] = kstest(x);

%% display
if dispFig > 0
    if dispFig == 1
        %% compare variable x distribution to normal distribution
        figure;
        % show distribution of x
        subplot(1,2,1);
        histogram(x);
        ylabel('Samples');
        title('Variable distribution');
        
        % compare distribution of x to normal distribution
        subplot(1,2,2);
        normplot(x);
    end
    
    %% message indicating whether data passes the normality test or not
    if h.adt == 0 &&...
            h.jbt == 0 &&...
            h.llt == 0 &&...
            h.kst == 0
        disp('The variable passed all the normality tests and can be considered as normal');
    elseif h.adt == 1 &&...
            h.jbt == 1 &&...
            h.llt == 1 &&...
            h.kst == 1
        disp('The variable failed all the normality tests and can be considered as non-normal. You may want to transform it.');
    else
        disp('Diagnostic not clear, the variable passed some normality tests but also failed some others => consider transforming but not mandatory.');
    end

end % display

end % function