clear all
cd ..
addpath([cd '/GLPreplicationWeb']);
addpath([cd '/GLPreplicationWeb/subroutines']);

% load the data and assign variable names
LoadDataAssignVariableNames

% variables included in the VAR
y=[unem nunem corepceinfl gdpinfl 100*(rgdp-log(pop)) 100*(h-log(pop)) aheinfl ls EBP];
series=["unemployment","natural unemployment","core inflation","inflation","GDP","hours","wage inflation (PNSE)","labor share","excess bond premium","wage inflation (total economy)"]
YLABEL=["percentage points","percentage points","percentage points","percentage points","percent","percent","percentage points","percent","percentage points","percentage points"];

% dimensions and settings
[junk maxind]=max(1-isnan(y));
[junk,n]=size(y);
lags=4;             % # lags
H=20;               % maximum horizon for impulse responses
M=20000;            % # MCMC draws (discard the first M/2)
nshock=9;           % position of the shock (in Cholesky identification)

% first sample
startdate1=1960-lags/4; T01aux=find(Time==startdate1); T01=max([T01aux,max(maxind)]);
if T01~=T01aux; startdate1=Time(T01); end
enddate1=1989.75; T11=find(Time==enddate1);

% second sample
startdate2=1990-lags/4; T02=find(Time==startdate2);
enddate2=2019.5; T12=find(Time==enddate2);

% GLP estimation, sample 1
resGLP1 = bvarGLP([y(T01:T11,:)],lags,'mcmc',1,'MCMCconst',15,'MNpsi',0,'noc',0,'sur',0,'Ndraws',M);

% impulse responses, sample 1
ndraws = size(resGLP1.mcmc.beta,3);
Dirf1 = zeros(H+1,size(y,2),ndraws);
Share1=zeros(ndraws,n);
VBC1=zeros(ndraws,n);
for jg = 1:ndraws
    Dirf1(:,:,jg) =  bvarIrfs(resGLP1.mcmc.beta(:,:,jg),resGLP1.mcmc.sigma(:,:,jg),nshock,H+1);
end
iw1=squeeze(Dirf1(:,8,:)-Dirf1(:,6,:)+Dirf1(:,5,:)+cumsum(Dirf1(:,4,:)/4));
Dirf1(:,end+1,:)=4*(iw1-lag(iw1));      % implied response of annualized wage inflation (total economy)
Dirf1=Dirf1/median(Dirf1(1,nshock,:));  % normalization of the size of the initial impulse
sIRF1 = sort(Dirf1,3);

% GLP estimation, sample 2
resGLP2 = bvarGLP([y(T02:T12,:)],lags,'mcmc',1,'MCMCconst',15,'MNpsi',0,'noc',0,'sur',0,'Ndraws',M);

% impulse responses, sample 2
Dirf2 = zeros(H+1,size(y,2),ndraws);
Share2=zeros(ndraws,n);
VBC2=zeros(ndraws,n);
for jg = 1:ndraws
    Dirf2(:,:,jg) =  bvarIrfs(resGLP2.mcmc.beta(:,:,jg),resGLP2.mcmc.sigma(:,:,jg),nshock,H+1);
end
iw2=squeeze(Dirf2(:,8,:)-Dirf2(:,6,:)+Dirf2(:,5,:)+cumsum(Dirf2(:,4,:)/4));
Dirf2(:,end+1,:)=4*(iw2-lag(iw2));      % implied response of annualized wage inflation (total economy)
Dirf2=Dirf2/median(Dirf2(1,nshock,:));  % normalization of the size of the initial impulse
sIRF2 = sort(Dirf2,3);

% plot the responses of all variables
qqq=[.025 .16 .5 .84 .975];         % percentiles of the posterior distribution
indfig=[1 3 4 7 10 8 5 6 9];        % variable position of responses to plot

figure('Position', [0, 0, 700, 600]);
count=0;
for jn = indfig
    count=count+1;
    subplot(3,3,count)
    quantilePlot([0:H]', squeeze(sIRF1(:,jn,round(qqq*ndraws)))); hold on; grid on;
    quantilePlot([0:H]', squeeze(sIRF2(:,jn,round(qqq*ndraws))),[.8941, .1020, .1098]);
    line([0 H],[0 0],'color','k')
    xlabel('horizon')
    ylabel(YLABEL(jn))
    title(series(jn));
end
aux=subplot(3,3,1); plots=get(gca, 'Children'); legend(plots([6 3]),{'pre 1990','post 1990'},'Location','northeast');