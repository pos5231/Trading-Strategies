function [y] = msharper(a, mu,Q,rf,weightfun)
size(a);
size(mu);
size(Q);
size(weightfun(a));
g = weightfun(a);
y = (mu*g'-rf)./sqrt(g*Q*g');
end