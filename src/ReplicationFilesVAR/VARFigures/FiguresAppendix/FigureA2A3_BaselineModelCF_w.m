clear all
cd ..
cd ..
addpath([cd '/GLPreplicationWeb']);
addpath([cd '/GLPreplicationWeb/subroutines']);

% load the data and assign variable names
LoadDataAssignVariableNames
load medU

% variables included in the VAR
y=[unem nunem corepceinfl gdpinfl 100*(rgdp-log(pop)) 100*(h-log(pop)) aheinfl ls eciWSPIinfl];
series=["unemployment","natural unemployment","core inflation","inflation","GDP","hours","wage inflation (PNSE)","labor share","wage inflation (ECI)","wage inflation (total economy)"]
YLABEL=["percentage points","percentage points","percentage points","percentage points","percent","percent","percentage points","percent","percentage points","percentage points"];

% dimensions and settings
[junk maxind]=max(1-isnan(y));
[junk,n]=size(y);
lags=4;             % # lags
H=20;               % maximum horizon for impulse responses
M=20000;            % # MCMC draws (discard the first M/2)
nshock=1;           % position of the shock (in Cholesky identification)

% first sample
startdate1=1990-lags/4; T01aux=find(Time==startdate1); T01=max([T01aux,max(maxind)]);
if T01~=T01aux; startdate1=Time(T01); end
enddate1=2007.75; T11=find(Time==enddate1);

% second sample
startdate2=1990-lags/4; T02=find(Time==startdate2);
enddate2=2019.5; T12=find(Time==enddate2);

% GLP estimation, sample 1
resGLP1 = bvarGLP([y(T01:T11,:)],lags,'mcmc',1,'MCMCconst',30,'MNpsi',0,'noc',0,'sur',0,'Ndraws',M);
ndraws = size(resGLP1.mcmc.beta,3);
    
% run disturbance smoother, sample 1 
rng(10)
CF1=zeros(n,H+1,ndraws);
for j=1:ndraws
    varZ=zeros(1,n*lags); varZ(1,1)=1;
    varc=0;
    varG=0;
    B=squeeze(resGLP1.mcmc.beta(:,:,j));
    varC=zeros(n*lags,1);
    varT=[B(2:end,:)';[eye(n*(lags-1)) zeros(n*(lags-1),n)]];
    varH=zeros(n*lags,n); varH(1:n,1:n)=chol(squeeze(resGLP1.mcmc.sigma(:,:,j)))';
    
    s00=zeros(n*lags,1);
    P00=zeros(n*lags,n*lags);
    
    DrawStates=DisturbanceSmootherVAR(medU,varc,varZ,varG,varC,varT,varH,s00,P00,H+1,1,n*lags,n,'smoother');
    CF1(:,:,j)=DrawStates(1:n,:);    
end
CF1=permute(CF1,[2,1,3]);
iw1=squeeze(CF1(:,8,:)-CF1(:,6,:)+CF1(:,5,:)+cumsum(CF1(:,4,:)/4));
CF1(:,end+1,:)=4*(iw1-lag(iw1));        % implied response of annualized wage inflation (total economy)
sCF1 = sort(CF1,3);

% GLP estimation, sample 2
resGLP2 = bvarGLP([y(T02:T12,:)],lags,'mcmc',1,'MCMCconst',30,'MNpsi',0,'noc',0,'sur',0,'Ndraws',M);

% run disturbance smoother, sample 2
rng(10)
CF2=zeros(n,H+1,ndraws);
for j=1:ndraws
    varZ=zeros(1,n*lags); varZ(1,1)=1;
    varc=0;
    varG=0;
    B=squeeze(resGLP2.mcmc.beta(:,:,j));
    varC=zeros(n*lags,1);
    varT=[B(2:end,:)';[eye(n*(lags-1)) zeros(n*(lags-1),n)]];
    varH=zeros(n*lags,n); varH(1:n,1:n)=chol(squeeze(resGLP2.mcmc.sigma(:,:,j)))';
    
    s00=zeros(n*lags,1);
    P00=zeros(n*lags,n*lags);
    
    DrawStates=DisturbanceSmootherVAR(medU,varc,varZ,varG,varC,varT,varH,s00,P00,H+1,1,n*lags,n,'smoother');
    CF2(:,:,j)=DrawStates(1:n,:);    
end
CF2=permute(CF2,[2,1,3]);
iw2=squeeze(CF2(:,8,:)-CF2(:,6,:)+CF2(:,5,:)+cumsum(CF2(:,4,:)/4));
CF2(:,end+1,:)=4*(iw2-lag(iw2));        % implied response of annualized wage inflation (total economy)
sCF2 = sort(CF2,3);

% plotg the impulse responses of unemployment and wage inflation
qqq=[.025 .16 .5 .84 .975];     % percentiles of the posterior distribution
indfig=[1 9];                   % variable position of impulse responses to plot

figure('Position', [0, 0, 1050, 300]);
subplot(1,3,1);
quantilePlot([0:H]', squeeze(sCF2(:,1,round(qqq*ndraws)))); hold on; grid on;
quantilePlot([0:H]', squeeze(sCF2(:,1,round(qqq*ndraws))),[.8941, .1020, .1098]);
line([0 H],[0 0],'color','k');
xlabel('horizon')
ylabel(YLABEL(1))
title(series(1));
plot([0:H]', medU,'--','LineWidth',3,'color',[.2157, .4941, .7216]);

subplot(1,3,2);
quantilePlot([0:H]', squeeze(sCF2(:,9,round(qqq*ndraws)))); hold on; grid on;
quantilePlot([0:H]', squeeze(sCF2(:,7,round(qqq*ndraws))),[.8941, .1020, .1098]); 
line([0 H],[0 0],'color','k');
xlabel('horizon')
ylabel(YLABEL(7))
title('wage inflation (PNSE and ECI)');
plots=get(gca, 'Children'); legend(plots([6 3]),{'post 1990 - ECI','post 1990 - PNSE'},'Location','southwest');

subplot(1,3,3);
quantilePlot([0:H]', squeeze(sCF2(:,9,round(qqq*ndraws)))); hold on; grid on;
quantilePlot([0:H]', squeeze(sCF2(:,10,round(qqq*ndraws))),[.8941, .1020, .1098]); 
line([0 H],[0 0],'color','k');
xlabel('horizon')
ylabel(YLABEL(7))
title('wage inflation (total economy and ECI)');
plots=get(gca, 'Children'); legend(plots([6 3]),{'post 1990 - total economy','post 1990 - PNSE'},'Location','southwest');

figure('Position', [0, 0, 700, 300]);
count=0;
for jn = indfig
    count=count+1;
    subplot(1,2,count)
    quantilePlot([0:H]', squeeze(sCF1(:,jn,round(qqq*ndraws)))); hold on; grid on;
    quantilePlot([0:H]', squeeze(sCF2(:,jn,round(qqq*ndraws))),[.8941, .1020, .1098]);
    line([0 H],[0 0],'color','k')
    xlabel('horizon')
    ylabel(YLABEL(jn))
    title(series(jn));
end
aux=subplot(1,2,1); plot([0:H]', medU,'--','LineWidth',3,'color',[.2157, .4941, .7216]); 
plots=get(gca, 'Children'); legend(plots([7 4]),{'1990 - 2007','1990 - 2019'},'Location','southwest');