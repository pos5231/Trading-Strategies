%%
tic
clear
clc
format long g

%% Import Data From File (Skip if second time)

%Uncomment the bottom line if used first time

%masterdata = xlsread('in_sample_data.xlsx');
%save masterdata.mat

warning('off')

%% Importing Data From Matrix (Save time)
% Saving the improted
tic
masterdata = importdata('masterdata.mat');


date = readtable('in_sample_data.xlsx','Sheet','in_sample_data','Range','A:A');

so = masterdata(:,2:6:end);     %Open Price
sh = masterdata(:,3:6:end);     %Daily High
sl = masterdata(:,4:6:end);     %Daily Low
sc = masterdata(:,5:6:end);     %Daily closing price
tvl = masterdata(:,6:6:end);    %Daily Volume
ind = masterdata(:,7:6:end);    %Movement
t = length(so); %Counting the number of trading days
n = size(so,2); %Counting the number of assets

rf = 0.04; % Risk-free rate

%% Part 1
rcc = sc(2:end,:)./sc(1:end-1,:) - 1; %calculating daily returns
rcc = [zeros(1,n);rcc]; %returns on day 1 = 0
avrcc = mean(rcc,2); %equally weighted average daily returns
Q = cov(rcc); %computing the covariance matrix
mu = mean(rcc); %computing the average returns

w_1 = (-1/n)*(rcc(2:end,:)-avrcc(2:end));
w_1 = [zeros(2,n); w_1];
w_1(end,:) =[];

r1 = sum(w_1.*rcc,2);
w1 = sum(abs(w_1));

rp1 = r1./w1; %Daily returns as per this strategy
dailyret_p1 = sum(rp1,2);
cumret_p1 = cumsum(dailyret_p1);

%% Part 2 - Optimizing Over a Set of Parameters

rco = sc(2:end,:)./so(1:end-1,:) - 1; %close to open return
rco = [zeros(1,n);rco];               %open to close return
roc = so(2:end,:)./sc(1:end-1,:) - 1; %open to close return
roc = [zeros(1,n);roc];               %open to close return
roo = so(2:end,:)./sc(1:end-1,:) - 1; %open to close return
roo = [zeros(1,n);roo];               %open to close return

avrco = mean(rco,2);                  %average close to open returns
avroo = mean(roo,2);                  %avergae open to open returns
avroc = mean(roc,2);                  %average open to close returns

e = 1/(4*log(2));

rvp = e*((log(sh(1:end,1:end)))-log(sl(1:end,1:end))).^2; %as mentioned

avrtvl = zeros(t,n);
avrrvp = avrtvl;

for i = 2:t
    for j = 1:n
        if i < 201
            avrtvl(i,j) = mean(tvl(1:i-1,j));
            avrrvp(i,j) = mean(rvp(1:i-1,j));
        else
            avrtvl(i,j) = mean(tvl(i-200:i-1,j));
            avrrvp(i,j) = mean(rvp(i-200:i-1,j));
        end
        
    end
end


Beq = 1;
a0 = rand(1,12);
options = optimoptions(@fmincon,'Display','off');

tic

mu = mean(rcc);
Q = cov(rcc);

for i = 3:t  % trial
    
%     mu = rcc(i-1,:);
%     Q = cov(rcc(i-2:i-1,:));   
    
    c1 = 1/n*sum(rcc(i-1,1:end)-avrcc(i-1));
    c2 = 1/n*sum(roo(i,1:end)-avroo(i));
    c3 = 1/n*sum(roc(i-1,1:end)-avroc(i-1));
    c4 = 1/n*sum(rco(i,1:end)-avrco(i));
    c5 = 1/n*sum((tvl(i-1,1:end)/avrtvl(i-1,1:end)).*(rcc(i-1,1:end)-avrcc(i-1)));
    c6 = 1/n*sum((tvl(i-1,1:end)/avrtvl(i-1,1:end)).*(roo(i-1,1:end)-avroo(i)));
    c7 = 1/n*sum((tvl(i-1,1:end)/avrtvl(i-1,1:end)).*(roc(i-1,1:end)-avroc(i-1)));
    c8 = 1/n*sum((tvl(i-1,1:end)/avrtvl(i-1,1:end)).*(rco(i,1:end)-avrco(i)));
    c9 = 1/n*sum((rvp(i-1,1:end)/avrrvp(i-1,1:end)).*(rcc(i-1,1:end)-avrcc(i-1)));
    c10 = 1/n*sum((rvp(i-1,1:end)/avrrvp(i-1,1:end)).*(roo(i-1,1:end)-avroo(i)));
    c11 = 1/n*sum((rvp(i-1,1:end)/avrrvp(i-1,1:end)).*(roc(i-1,1:end)-avroc(i-1)));
    c12 = 1/n*sum((rvp(i-1,1:end)/avrrvp(i-1,1:end)).*(rco(i,1:end)-avrco(i)));
    
    
    Aeq = [c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12];
    
    a_opt(i,:) = fmincon(@(a) longsharpe(a, mu,Q,  avrcc, avrco, avroc, avroo, avrrvp,...
        avrtvl, rcc, rco, roc, roo, tvl, rvp, i, n),a0,[],[],Aeq,Beq,[],[],[],options);
    
    w_strat2(i,:) = w2(a_opt,   avrcc, avrco, avroc, avroo, avrrvp,...
        avrtvl, rcc, rco, roc, roo, tvl, rvp, i, n);    
end
toc

%%
w_2 = normalize(w_strat2);
r2 = sum(w_2.*roc,2);
w2 = sum(abs(w_2));
rp2 = r2./w2;
dailyret_p2 = sum(rp2,2);
cumret_p2 = cumsum(dailyret_p2);


%% Part 4.1 Maximize Info. Ratio


% lb= -1*ones(1,n);
% ub= ones(1,n);
% Aeq = ones(1,n);
% Beq = 0;
% x0 = 1/n*ones(1,n);
%
% tic
% for i = 3:t
%
%     if i > 30
%         ret4 = rcc(i-30:i,:);
%         mu = mean(ret4);
%         Q = cov(ret4);
%     else
%         ret4 = rcc(i-1:i,:);
%         mu = mean(ret4);
%         Q = cov(ret4);
%     end
%
%     fun = @(x) -(mu*x'-0.15/252)/sqrt(x*(Q*x'));
%     [x] = fmincon(fun,x0,[],[],Aeq,Beq,lb,ub,[],options);
%     w_4(i,:) = x;
%     clear x
%     clear mu
%
% end
% toc

%%
% r4 = sum(w_4.*rcc,2);
% w4 = sum(abs(w_4));
% rp4 = r4./w1;
% dailyret_p4 = sum(rp4,2);
% cumret_p4 = cumsum(dailyret_p4);


%% Functions

function [y] = longsharpe(a, mu,Q,  avrcc, avrco, avroc, avroo, avrrvp,...
    avrtvl, rcc, rco, roc, roo, tvl, rvp, i, n)

% find weight vector for each day

for j = 1:n
        w(1,j) = a(1)*(rcc(i-1,j)-avrcc(i-1))./n...
        + a(2)*(roo(i,j)-avroo(i))./n...
        + a(3)*(roc(i-1,j)-avroc(i-1))./n ...
        + a(4)*(rco(i,j)-avrco(i))./n...
        + a(5)*(tvl(i-1,j)/avrtvl(i-1,j))*(rcc(i-1,j)-avrcc(i-1))./n...
        + a(6)*(tvl(i-1,j)/avrtvl(i-1,j))*(roo(i,j)-avroo(i))./n...
        + a(7)*(tvl(i-1,j)/avrtvl(i-1,j))*(roc(i,j)-avroc(i-1))./n...
        + a(8)*(tvl(i-1,j)/avrtvl(i-1,j))*(rco(i,j)-avrco(i))./n...
        + a(9)*(rvp(i-1,j)/avrrvp(i-1,j))*(rcc(i-1,j)-avrcc(i-1))./n...
        + a(10)*(rvp(i-1,j)/avrrvp(i-1,j))*(roo(i,j)-avroo(i))./n...
        + a(11)*(rvp(i-1,j)/avrrvp(i-1,j))*(roc(i-1,j)-avroc(i-1))./n...
        + a(12)*(rvp(i-1,j)/avrrvp(i-1,j))*(rco(i,j)-avrco(i))./n;
end

%setting up sharpe ratio with risk free rate 4%
y = -(mu*w'-0.04/252)/sqrt(w*(Q*w'));

end