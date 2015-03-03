clear all;
close all;
[Y,Ynames]=xlsread('2008crisisdate.xls');
% data sets must be in format such that every row is a time period

asvar = {'SPX'; 'Bond'; 'MXEF'; 'JPY'; 'GOLD'; 'OIL'; '10YR'};    % variable names
nlag = 1;                   % lags

%format data to returns
PRICES=Y;
T=size(Y,1);
Y=100*(Y(2:T,:)./Y(1:T-1,:)-1);
Y(:,7)=-Y(:,7); %bonds quoted in yields
Y(~Y)=.0001;

%interested in exploring period around the financial crisis
x1=5200;x2=5600;
Y=Y(x1:x2,:);
dates=Ynames(x1+2:x2+2,1);

lt_setvar('data', Y, asvar, nlag);  % set data

lt_setvar('ranseed', 5);       % set ranseed
lt_setvar('intercept', 1);     % set time-varying intercept
lt_setvar('LTb', 1); % set threshold
lt_setvar('LTa', 1); %set threshold