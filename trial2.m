
lb=-inf*ones(1,12);
ub=ones(1,100);
Aeq = ones(1,12);
Beq = 1;

w0=1/1000*rand(12,1);

options = optimoptions(@fmincon,'Display','off');

for i = 3
weightfun=@(a) w2(a,   avrcc, avrco, avroc, avroo, avrrvp,...
    avrtvl, rcc, rco, roc, roo, tvl, rvp, i, j, n);

[a_opt,f] = fmincon(@(a) msharper(a, mu,Q,rf,weightfun),w0,[],[],Aeq,Beq,[],[]);
end