%%
tic
clear
clc
format long g

%% Import Data From File (Skip if second time)

%Uncomment the bottom line if used first time

%masterdata = xlsread('in_sample_data.xlsx');
%save masterdata.mat

%% Importing Data From Matrix (Save time)
% Saving the improted
tic
masterdata = importdata('masterdata.mat');

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

rvp = e*((log(sh(1:end,1:end)))-log(sl(1:end,1:end))).^2;

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


lb=-inf*ones(1,n);
ub=ones(1,n);
Amat=ones(1,n);
Bmat=1;

w0=1/n*ones(n,1);

options = optimoptions(@fmincon,'Display','off');

fun = @(x) -(mu*x-rf)/sqrt(x'*Q*x); %sharpe-ratio function

W = zeros(t,n);

for i = 3:t
        
    for j = 1:n
        lcon = @(a) wconstraint(a, avrcc, avrco, avroc, avroo, avrrvp, avrtvl,rcc,rco,roc,roo, tvl, rvp, i, j, n);
    end % setting up the linear constraint
    
    [x,f] = fmincon(fun,w0,[],[],Amat,Bmat,lb,ub,lcon,options) %solving for max-sharpe ratio
    
    W(i,:) = x';
    
end

r2 = sum(W.*rcc,2);
w_2 = sum(abs(w_1));
rp2 = r2./w_2; %Daily returns as per this strategy
dailyret_p2 = sum(rp2,2);
cumret_p2 = cumsum(dailyret_p2);

%% Part 3 - Imposing Restrictions on Trade Directions

w3 = w_2;

fill3 = zeros(t,j);
filltest_3 = w3.*ind;

for i = 1:t
    for j = 1:n
        if filltest_3(i,j) > 0
            fill3(i,j) = filltest_3(i,j);
        else
            fill3(i,j) = 0;
        end
    end
end

r3 = sum(fill3.*w3.*roc,2);
w_3 = sum(abs(w3));
rp3 = r3./w_3;
dailyret_p3 = sum(rp3,2);
cumret_p3 = cumsum(dailyret_p3);

%% Part 4 - Generalized Portfolio Weights
% 4.1 - Maximum Information Ratio
rspy = 0.08;
clear x f

IR = @(x) (mu*x-rspy)/sqrt(x'*Q*x); %sharpe-ratio function

for i = 3:t
           
    [x,f] = fmincon(IR,w0,[],[],Amat,Bmat,lb,ub,[],options); %solving for max-sharpe ratio    
    W(i,:) = x';
    
end

%% Plot
plot([1:t],cumret_p2)
hold on
plot([1:t],cumret_p3)
toc


