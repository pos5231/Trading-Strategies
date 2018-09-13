function [ceq] = confuneq(a,avrcc, avrco, avroc, avroo, avrrvp,...
    avrtvl, rcc, rco, roc, roo, tvl, rvp, i, n)

for j = 1:n
    y(1,j) = a(1)*(rcc(i-1,j)-avrcc(i-1))./n...
        + a(2)*(roo(i,j)-avroo(i))./n...
        + a(3)*(roc(i-1,j)-avroc(i-1))./n ...
        + a(4)*(rco(i,j)-avrco(i))./n...
        + a(5)*(tvl(i-1,j)/avrtvl(i-1,j))*(rcc(i-1,j)-avrcc(i-1))./n...
        + a(6)*(tvl(i-1,j)/avrtvl(i-1, j))*(roo(i,j)-avroo(i))./n...
        + a(7)*(tvl(i-1,j)/avrtvl(i-1,j))*(roc(i,j)-avroc(i-1))./n...
        + a(8)*(tvl(i-1,j)/avrtvl(i-1,j))*(rco(i,j)-avrco(i))./n...
        + a(9)*(rvp(i-1,j)/avrrvp(i-1,j))*(rcc(i-1,j)-avrcc(i-1))./n...
        + a(10)*(rvp(i-1,j)/avrrvp(i-1,j))*(roo(i,j)-avroo(i))./n...
        + a(11)*(rvp(i-1,j)/avrrvp(i-1,j))*(roc(i-1,j)-avroc(i-1))./n...
        + a(12)*(rvp(i-1,j)/avrrvp(i-1,j))*(rco(i,j)-avrco(i))./n;
end
a
% linear equality constraints
ceq = sum(y) - 1

end
