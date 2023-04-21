function thCMRT = EZ2cmrt(nu, z, a, Ter, s) 

%     this is optimized code for the expression
%     ((exp(4*nu/s^2*a)+exp(2*nu*(x+a)/s^2)-exp(2*nu/s^2*a)-exp(2*nu*x/s^2))*x+(2*exp(2*nu/s^2*a)-2*exp(2*nu*(x+a)/s^2))*a)/nu/(exp(2*nu/s^2*a)-exp(2*nu*x/s^2))/(-1+exp(2*nu/s^2*a))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the routine implements the error response time variance, while most
% people will assume the correct response time variance
z = a - z;
nu = -nu;

if abs(nu) <= 10^(-4) && nu~=0; nu=sign(nu)*10^(-4); 
elseif nu == 0; nu = 10^(-4);
end;

    
expr1 = z*(intermedcalc(z-a,a,nu,s) + intermedcalc(0,z,nu,s)) + 2*a*intermedcalc(z,0,nu,s);
expr2 = nu*intermedcalc(z,a,nu,s)*intermedcalc(-a,0,nu,s);

if abs(expr2) <= 10^(-6) && expr2~=0; expr2=sign(expr2)*10^(-6); 
elseif expr2 == 0; expr2 = 10^(-6);
end;

thCMRT = expr1/expr2 + Ter;