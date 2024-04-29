function [p1,p1c,p2]=multskewkurt(X,alpha)

if nargin < 2, 
    alpha = 0.05;  %(default)
end 

if nargin < 1, 
    error('Requires at least one input argument.'); 
end 

[n,p] = size(X);

% calculate mean and covariance
mat_mean=mean(X,1);
X_hat=zeros(n,p);

for i=1:n,   
    X_hat(i,:)=X(i,:)-mat_mean;
end

X_cov=cov(X,1);
D=X_hat*inv(X_cov)*X_hat';

b1p = (sum(sum(D.^3)))/n^2;                    % multivariate skewness coefficient
b2p=trace(D.^2)/n;                             % multivariate kurtosis coefficient 

k = ((p+1)*(n+1)*(n+3))/(n*(((n+1)*(p+1))-6)); % small sample correction
v = (p*(p+1)*(p+2))/6;                         % degrees of freedom
g1c = (n*b1p*k)/6;                             % skewness test statistic corrected for small sample:it approximates to a chi-square distribution
g1 = (n*b1p)/6;                                % skewness test statistic:it approximates to a chi-square distribution
p1 = 1 - chi2cd(g1,v);                         % significance value associated to the skewness
p1c = 1 - chi2cd(g1c,v);                       % significance value associated to the skewness corrected for small sample

g2 = (b2p-(p*(p+2)))/(sqrt((8*p*(p+2))/n));    % kurtosis test statistic: it approximates to
                                               % a unit-normal distribution                                                                                        
p2 = 2*(1-normcd(abs(g2)));                    % significance value associated to the kurtosis

end