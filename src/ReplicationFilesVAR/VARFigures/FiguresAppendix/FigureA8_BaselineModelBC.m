clear all
cd ..
cd ..
addpath([cd '/GLPreplicationWeb']);
addpath([cd '/GLPreplicationWeb/subroutines']);

% load the data and assign variable names
LoadDataAssignVariableNames

% variables included in the VAR
y=[unem nunem corepceinfl gdpinfl 100*(rgdp-log(pop)) 100*(h-log(pop)) aheinfl ls];
series=["unemployment","natural unemployment","core inflation","inflation","GDP","hours","wage inflation (PNSE)","labor share","wage inflation (total economy)"]
YLABEL=["percentage points","percentage points","percentage points","percentage points","percent","percent","percentage points","percent","percentage points"];

% dimensions and settings
[junk maxind]=max(1-isnan(y));
[junk,n]=size(y);
lags=4;             % # lags
H=20;               % maximum horizon for impulse responses
M=20000;            % # MCMC draws (discard the first M/2)
f1=2*pi/32;         % lower bound of BC frequency
f2=2*pi/6;          % upper bound of BC frequency

% first sample
startdate1=1960-lags/4; T01aux=find(Time==startdate1); T01=max([T01aux,max(maxind)]);
if T01~=T01aux; startdate1=Time(T01); end
enddate1=1989.75; T11=find(Time==enddate1);

% second sample
startdate2=1990-lags/4; T02=find(Time==startdate2);
enddate2=2019.5; T12=find(Time==enddate2);

% GLP estimation, sample 1
resGLP1 = bvarGLP([y(T01:T11,:)],lags,'mcmc',1,'MCMCconst',5,'MNpsi',0,'noc',0,'sur',0,'Ndraws',M);

% impulse responses, sample 1
ndraws = size(resGLP1.mcmc.beta,3);
Dirf1 = zeros(H+1,size(y,2),ndraws);
Share1=zeros(ndraws,n);
VBC1=zeros(ndraws,n);
for jg = 1:ndraws
    Dirf1(:,:,jg) =  bvarIrfsBC(resGLP1.mcmc.beta(:,:,jg),resGLP1.mcmc.sigma(:,:,jg),H+1,1,f1,f2);
end
iw1=squeeze(Dirf1(:,8,:)-Dirf1(:,6,:)+Dirf1(:,5,:)+cumsum(Dirf1(:,4,:)/4));
Dirf1(:,9,:)=4*(iw1-lag(iw1));      % implied response of annualized wage inflation (total economy)
sIRF1 = sort(Dirf1,3);

% GLP estimation, sample 2
resGLP2 = bvarGLP([y(T02:T12,:)],lags,'mcmc',1,'MCMCconst',5,'MNpsi',0,'noc',0,'sur',0,'Ndraws',M);

% impulse responses, sample 2
Dirf2 = zeros(H+1,size(y,2),ndraws);
Share2=zeros(ndraws,n);
VBC2=zeros(ndraws,n);
for jg = 1:ndraws
    Dirf2(:,:,jg) =  bvarIrfsBC(resGLP2.mcmc.beta(:,:,jg),resGLP2.mcmc.sigma(:,:,jg),H+1,1,f1,f2);
end
iw2=squeeze(Dirf2(:,8,:)-Dirf2(:,6,:)+Dirf2(:,5,:)+cumsum(Dirf2(:,4,:)/4));
Dirf2(:,9,:)=4*(iw2-lag(iw2));      % implied response of annualized wage inflation (total economy)
sIRF2 = sort(Dirf2,3);

% plot the impulse responses of unemployment and inflation
qqq=[.025 .16 .5 .84 .975];     % percentiles of the posterior distribution
indfig=[1 2 3 4 ];              % variable position of responses to plot

figure('Position', [0, 0, 700, 600]);
count=0;
for jn = indfig
    count=count+1;
    subplot(ceil(length(indfig)/2),2,count)
    quantilePlot([0:H]', squeeze(sIRF1(:,jn,round(qqq*ndraws)))); hold on; grid on;
    quantilePlot([0:H]', squeeze(sIRF2(:,jn,round(qqq*ndraws))),[.8941, .1020, .1098]);
    line([0 H],[0 0],'color','k')
    xlabel('horizon')
    ylabel(YLABEL(jn))
    title(series(jn));
end
aux=subplot(2,2,1); limitsy=get(aux,'ylim'); subplot(2,2,2); ylim(limitsy);
subplot(2,2,1); plots=get(gca, 'Children'); legend(plots([6 3]),{'pre 1990','post 1990'},'Location','southwest');