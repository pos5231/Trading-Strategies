clc
lb=[];
ub= [];
Aeq = ones(1,12);
Beq = 1;

a0=1/12*ones(1,12);

options = optimoptions(@fmincon,'Display','off');

for i = 7
    
    weightfun=@(a) w2(a,   avrcc, avrco, avroc, avroo, avrrvp,...
        avrtvl, rcc, rco, roc, roo, tvl, rvp, i, n);
    
    [a_opt,f] = fmincon(@(a) msharper(a, mu,Q,rf,weightfun),a0,[],[],Aeq,Beq,[],[],[],options);
       
    weight_opt = w2(a_opt,   avrcc, avrco, avroc, avroo, avrrvp,...
        avrtvl, rcc, rco, roc, roo, tvl, rvp, i, n);
    
    sum(weight_opt)
    
end