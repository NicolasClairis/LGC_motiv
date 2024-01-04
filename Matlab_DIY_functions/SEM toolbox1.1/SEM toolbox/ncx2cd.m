function p = ncx2cd(x,v,delta)
%NCX2CDF Noncentral chi-square cumulative distribution function (cdf).
%   P = NCX2CDF(X,V,DELTA) Returns the noncentral chi-square cdf with V 
%   degrees of freedom and noncentrality parameter, DELTA, at the values 
%   in X.
%
%   The size of P is the common size of the input arguments. A scalar input  
%   functions as a constant matrix of the same size as the other inputs.     
%
%   Some texts refer to this distribution as the generalized Rayleigh,
%   Rayleigh-Rice, or Rice distribution.
%
%   See also NCX2INV, NCX2PDF, NCX2RND, NCX2STAT, CHI2CDF, CDF.

%   Reference:
%      [1]  Evans, Merran, Hastings, Nicholas and Peacock, Brian,
%      "Statistical Distributions, Second Edition", Wiley
%      1993 p. 50-52.

%   Copyright 1993-2007 The MathWorks, Inc. 
%   $Revision: 2.15.4.8 $  $Date: 2007/05/23 19:15:55 $

if nargin <  3, 
    error('stats:ncx2cdf:TooFewInputs','Requires three input arguments.'); 
end

[errorcode x v delta] = distchck(3,x,v,delta);

if errorcode > 0
    error('stats:ncx2cdf:InputSizeMismatch',...
          'Requires non-scalar arguments to match in size.');
end

% Initialize P to zero.
if isa(x,'single') || isa(v,'single') || isa(delta,'single')
   p = zeros(size(x),'single');
   crit = eps('single');
else
   p = zeros(size(x));
   crit = eps;
end
p(isnan(x)) = NaN;
p(x==Inf) = 1;

p(x <= 0) = 0;  
p(delta < 0) = NaN;  % can't have negative non-centrality parameter.
p(v < 0) = NaN;      % can't have negative d.f.
k = (v==0) & (x==0) & (delta>0);
if any(k)            % 0 d.f. requires special treatment
    hdelta = delta(k)/2;
    p(k) = hdelta * exp(-hdelta);
end
k = find((x > 0) & (delta >=0) & isfinite(x));
hdelta = delta(k)/2;  % this is used in the loop, pre-compute.

v = v(k);
x = x(k);

% when non-centrality parameter is very large, the initial values of
% the poisson numbers used in the approximation are very small,
% smaller than epsilon. This would cause premature convergence. To
% avoid that, we start from counter=hdelta, which is the peak of the
% poisson numbers, and go in both directions.
counter = floor(hdelta);

% In principle we are going to sum terms of the form
%     poisspdf(counter,hdelta).*chi2cdf(x,v+2*counter)
% but we will compute these factors from scratch just once,
% and update them each time using a recurrence formula.

% Sum the series from the poisson peak upward
P = poisspd(counter,hdelta);
C = chi2cd(x,v+2*counter);
E = exp((v/2+counter-1).*log(x/2) - x/2 - gammaln(v/2+counter));
p = sumseries(P,C,E,counter,x,v,hdelta,k,p,crit,false);
   
% Now process in the other direction
counter = counter - 1; % set counter back just below the peak;
if ~any(counter>=0)
   return
end
p = sumseries(P,C,E,counter,x,v,hdelta,k,p,crit,true);

% For some points we may never have found anything to sum, because the
% poisson peak was away from the ncx2 peak, so try working up from zero
fromzero = (p(k)==0);
if any(fromzero)
    v = v(fromzero);
    x = x(fromzero);
    hdelta = hdelta(fromzero);
    k(~fromzero) = [];
    counter = zeros(size(x));
    p(k) = poisspd(0,hdelta) .* chi2cd(x,v);
    while(any(counter<hdelta))
        % As long as the series is increasing, roundoff errors are
        % magnified by using recursion so compute each term directly
        counter = counter+1;
        P = poisspd(counter,hdelta);
        C = chi2cd(x,v+2*counter);
        dp = P .* C;
        p(k) = p(k) + dp;
        if all(p(k)>0 & dp<crit*p(k))
            break
        end
    end
end

p = min(p,1);  % fix roundoff

% ------------------------------------------------------
function p = sumseries(P,C,E,counter,x,v,hdelta,k,p,crit,countdown)
  
if countdown
    % If summing downward, get valid indices and compute first values
    j = (counter>=0);
    x = x(j);
    v = v(j);
    hdelta = hdelta(j);
    k = k(j);
    counter = counter(j);
    P = P(j) .* (counter+1) ./ hdelta;
    E = E(j);
    C = C(j) + E;
end

while ~isempty(counter)
   % Compute product to add to the series
   pplus  = P.*C;
   j = isnan(pplus);
   if any(j(:))
      % Do not continue with NaN values
      if all(j(:)), break; end
      j = ~j;
      x = x(j);
      v = v(j);
      hdelta = hdelta(j);
      counter = counter(j);
      k = k(j);
      pplus = pplus(j);
   end
   p(k) = p(k) + pplus; % accumulate p for k indices
   
   % Test to see if we should continue
   j = find(pplus > p(k)*crit);
   if countdown
      j(counter(j)<0) = [];
   end
   if isempty(j) % no more computation needed in this direction
      break;
   end
   
   % Continue with the points for which we continue to make progress
   x = x(j);
   v = v(j);
   hdelta = hdelta(j);
   k = k(j);
   if countdown   % update terms recursively, counting down or up
      counter = counter(j) - 1;
      P = P(j) .* (counter+1) ./ hdelta;
      E = E(j) .* (v/2+counter+1) ./ (x/2);
      C = C(j) + E;
   else
      counter = counter(j) + 1;
      P = P(j) .* hdelta ./ counter;
      E = E(j) .* (x/2) ./ (v/2+counter-1);
      C = C(j) - E;
   end
end
