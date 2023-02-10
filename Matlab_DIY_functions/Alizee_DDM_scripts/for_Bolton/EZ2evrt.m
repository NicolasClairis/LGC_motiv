function thEVRT = EZ2evrt(nu, z, a, s) 

%     this is optimized code for the expression
%     (-4*nu*exp(2*nu/s^2*a)*(-1+exp(2*z*nu/s^2))*(exp(4*nu/s^2*a)-exp(2*z*nu/s^2))*a^2-4*exp(2*(z+a)*nu/s^2)*nu*(exp(2*nu/s^2*a)-1)^2*z^2+8*exp(2*(z+a)*nu/s^2)*nu*(exp(2*nu/s^2*a)-1)^2*a*z+2*s^2*exp(2*nu/s^2*a)*(-1+exp(2*z*nu/s^2))*(exp(2*nu/s^2*a)-1)*(-exp(2*nu/s^2*a)+exp(2*z*nu/s^2))*a-s^2*(exp(2*nu/s^2*a)-1)^2*(-exp(4*nu/s^2*a)+exp(4*z*nu/s^2))*z)/(exp(2*nu/s^2*a)-1)^2/nu^3/(exp(2*nu/s^2*a)-exp(2*z*nu/s^2))^2\

s2 = s^2;
if abs(nu) <= 10^(-4) && nu~=0; nu=sign(nu)*10^(-4); 
elseif nu == 0; nu = 10^(-4);
end;

expr1 = -2*a*intermedcalc(0,z,nu,s)*(2*nu*a*intermedcalc(z,2*a,nu,s) + ...
    s2*intermedcalc(0,a,nu,s)*intermedcalc(z,a,nu,s))*exp(2*nu*a/s2);
expr2 = nu^(3)*intermedcalc(0,a,nu,s)^(2)*intermedcalc(z,a,nu,s)^(2);

if abs(expr2) <= 10^(-6) && expr2~=0; expr2=sign(expr2)*10^(-6); 
elseif expr2 == 0; expr2 = 10^(-6);
end;

exprLEFT = expr1/expr2;

expr3 = 4*nu*z*(2*a-z)*exp(2*nu*(z+a)/s2) + z*s2 * intermedcalc(2*z,2*a,nu,s);
expr4 = nu^3*intermedcalc(z,a,nu,s)^(2);


if abs(expr4) <= 10^(-6) && expr4~=0; expr4=sign(expr4)*10^(-6); 
elseif expr4 == 0; expr4 = 10^(-6);
end;

exprRIGHT = expr3/expr4;

thEVRT = exprLEFT + exprRIGHT;
% if thEVRT <= 10^(-6); thEVRT = 10^(-6); end;