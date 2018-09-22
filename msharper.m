function [y] = msharper(a, mu,Q,rf,weightfun)
g = weightfun(a);
y = (mu*g'-rf)./sqrt(g*Q*g'); % scalar valur for objective function
end