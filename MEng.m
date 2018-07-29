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

w2 = zeros(t,n);

for i=3:t
    for j=1:n
        w2(i,j) = (-1/n)*(rcc(i-1,j)-avrcc(i-1));
    end
end

r2 = sum(w2.*rcc,2);
w1 = abs(w2);
w1 = sum(w1,2);

rp1 = r2./w1; %Daily returns as per this strategy

%% Part 2

rco = sc(2:end,:)./so(1:end-1,:) - 1; 
rco = [zeros(1,n);rco];
roc = so(2:end,:)./sc(1:end-1,:) - 1;
roc = [zeros(1,n);roc];
roo = so(2:end,:)./sc(1:end-1,:) - 1;
roo = [zeros(1,n);roo];
a = ones(1,12);

avrco = mean(rco,2); 
avroo = mean(roo,2);
avroc = mean(roc,2);


e = 1/(4*log(2));

for i = 1:t
    for j = 1:n
        rvp(i,j) = e*((log(sh(i,j)))-log(sl(i,j)))^2;
    end
end

avrtvl = zeros(t,n);
avrrvp = avrtvl;

for i = 3:t
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

options = optimoptions(@fmincon,'Display','off','MaxIterations',3000);


w2  = @(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12)...
            a1*(rcc(i-1,j)-avrcc(i-1))/n + a2*(roo(i,j)-avroo(i))/n...
            + a3*(roc(i-1,j)-avroc(i-1))/n + a4*(rco(i,j)-avrco(i))/n...
            + a5*(tvl(i-1,j)/avrtvl(i-1,j))*(rcc(i-1,j)-avrcc(i-1))/n...
            + a6*(tvl(i-1,j)/avrtvl(i-1,j))*(roo(i,j)-avroo(i))/n...
            + a7*(tvl(i-1,j)/avrtvl(i-1,j))*(roc(i,j)-avroc(i-1))/n...
            + a8*(tvl(i-1,j)/avrtvl(i-1,j))*(rco(i,j)-avrco(i))/n...
            + a9*(rvp(i-1,j)/avrrvp(i-1,j))*(rcc(i-1,j)-avrcc(i-1))/n...
            + a10*(rvp(i-1,j)/avrrvp(i-1,j))*(roo(i,j)-avroo(i))/n...
            + a11*(rvp(i-1,j)/avrrvp(i-1,j))*(roc(i-1,j)-avroc(i-1))/n...
            + a12*(rvp(i-1,j)/avrrvp(i-1,j))*(rco(i,j)-avrco(t))/n;

       
        
for i = 3:7
    for j = 1:n
                        
        f = @(w2) -((mu-rf)/swrt(w2*Q*w2));
        
        c = fmincon(f,w0,[],[],Amat,Bmat,lb,ub,[],options)';
        
    end
    
    w_2(i,:) = c;
end


%%

% for i = 3:t
%     for j = 1:n
%         w2(i,j) = a(1)*(rcc(i-1,j)-avrcc(i-1))/n + a(2)*(roo(i,j)-avroo(i))/n...
%             + a(3)*(roc(i-1,j)-avroc(i-1))/n - a(4)*(rco(i,j)-avrco(i))/n;
%            


%     end
% end
               
%% Part 3

toc