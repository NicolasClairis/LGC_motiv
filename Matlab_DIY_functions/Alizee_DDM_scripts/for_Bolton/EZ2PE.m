function thPE = EZ2PE(v, x, a, s)
% v is drift
% x is starting point
% a is boundary
% s is random disturbance (standard deviation of noise)


if (exp(-2 * v * a / (s * s)) - 1) == 0;
     thPE = (exp(-2 * v * a / (s * s)) - exp(-2 * v * x / (s * s) )) / (10^(-6));
else
    thPE = (exp(-2 * v * a / (s * s)) - exp(-2 * v * x / (s * s) )) / (exp(-2 * v * a / (s * s)) - 1);
end