% Computes full predictive density 
% using the BVAR of Giannone, Lenza and Primiceri (2012)
% Based on the default setting, including all priors
% %%%%%%%%%%%%%%%%%%%%
% The model includes the following data 
% RGDP: 4 x logarithm of Real Gross Domestic Product, Quantity Index (2000=100) , SAAR
% PGDP: 4 x logarithm of Gross domestic product Price Index
% Cons: 4 x logarithm of Real Personal Consumption Expenditures, Quantity Index (2000=100) , SAAR
% GPDInv: 4 x logarithm of Real Gross Private Domestic Investment, Quantity Index (2000=100) , SAAR
% Emp. Hours: 4 x logarithm of HOURS OF ALL PERSONS: NONFARM BUSINESS SEC (1982=100,SA)
% Real Comp/Hour: 4 x logarithm of REAL COMPENSATION PER HOUR,EMPLOYEES:NONFARM BUSINESS(82=100,SA)
% FedFunds: INTEREST RATE: FEDERAL FUNDS (EFFECTIVE) (% PER ANNUM,NSA)

clear all
close all
addpath([cd '/subroutines'])  %on a MAC
%addpath([cd '\subroutines']) %on a PC

load DataSW 
% Load data from the dataset of Stock and Watson (2008)
% The variables enter the models in annualized log-levels (i.e. we take logs and multiply by 4), 
% except those already defined in terms of annualized rates, such as
% interest rates, which are taken in levels and are divided by 100.



lags = 5;

% Run the Bayesian VAR
res = bvarGLP(y,lags,'mcmc',1,'MCMCconst',1.6);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Arrange and display some results

hz = res.setpriors.hz; % forecast horizon

FCSTmode = res.postmax.forecast(:,1); %forecasts at the mode (for GDP)

FCSTmode = diag(1./hz)*(res.postmax.forecast(:,1)-y(end,1))*100; %forecasts at the mode (for average GDP growth, annualized)
FCSTmcmc = diag(1./hz)*squeeze(res.mcmc.Dforecast(:,1,:)-y(end,1))*100; %predictive density 


clc
disp('---------------------------------------')
disp(['Forecasting in ',datestr(Time(end),'QQ-YY'),' (GDP annualized average growth)'])
disp(' Qrts ahead    Point   Predictive Density')
disp('            (at Mode)    Mean       Std ')
disp([hz([1 4 8])' FCSTmode([1 4 8],:) median(FCSTmcmc([1 4 8],:),2) std(FCSTmcmc([1 4 8],:),[],2)])


%Forecasting the annualized average quarterly growth of GDP
figure(1)
subplot(2,1,1);
hor = hz(1); % one step ahead
hist(FCSTmcmc(hor,:));
title(['Forecasting in ',datestr(Time(end),'QQ-YY'),': Avg GDP growth ',num2str(hor),' qrts ahead. (Fcst at the mode: ',num2str(FCSTmode(hor,1)),')'])
axis([-15 10 0 inf])

subplot(2,1,2);
hor = hz(4); % one step ahead
hist(FCSTmcmc(hor,:));
title(['Forecasting in ',datestr(Time(end),'QQ-YY'),': Avg GDP growth ',num2str(hor),' qrts ahead. (Fcst at the mode: ',num2str(FCSTmode(hor,1)),')'])
axis([-15 10 0 inf])

