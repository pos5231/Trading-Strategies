function [y] = msharper(a, mu,Q,rf,weightfun)
y = (mu'*weightfun(a)-rf)%./sqrt(weightfun(a)'*Q*weightfun(a));
end