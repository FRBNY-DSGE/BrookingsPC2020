clear all
cd ..
cd ..
addpath([cd '/GLPreplicationWeb']);
addpath([cd '/GLPreplicationWeb/subroutines']);

% load the data and assign variable names
LoadDataAssignVariableNames

% variables included in the VAR
y=[unem nunem corepceinfl gdpinfl 100*(rgdp-log(pop)) 100*(h-log(pop)) aheinfl ls lsnfbs lsnfcs];
series=["unemployment","natural unemployment","core inflation","inflation","GDP","hours","wage inflation (PNSE)","labor share","labor share (NFBS)","labor share (NFCS)"]
YLABEL=["percentage points","percentage points","percentage points","percentage points","percent","percent","percentage points","percent","percent","percent"];

% dimensions and settings
[junk maxind]=max(1-isnan(y));
[junk,n]=size(y);
lags=4;             % # lags
H=20;               % maximum horizon for impulse responses
M=20000;            % # MCMC draws (discard the first M/2)
nshock=1;           % position of the shock (in Cholesky identification)

% first sample
startdate1=1960-lags/4; T01aux=find(Time==startdate1); T01=max([T01aux,max(maxind)]);
if T01~=T01aux; startdate1=Time(T01); end
enddate1=1989.75; T11=find(Time==enddate1);

% second sample
startdate2=1990-lags/4; T02=find(Time==startdate2);
enddate2=2019.5; T12=find(Time==enddate2);

% GLP estimation, sample 1
resGLP1 = bvarGLP([y(T01:T11,:)],lags,'mcmc',1,'MCMCconst',50,'MNpsi',0,'noc',0,'sur',0,'Ndraws',M);
    
% compute median IRF of U, sample 1
ndraws = size(resGLP1.mcmc.beta,3);
Dirf1 = zeros(H+1,size(y,2),ndraws);
for jg = 1:ndraws
    Dirf1(:,:,jg) =  bvarIrfs(resGLP1.mcmc.beta(:,:,jg),resGLP1.mcmc.sigma(:,:,jg),nshock,H+1);
end
medU=median(squeeze(Dirf1(:,1,:))')';

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
sCF1 = sort(CF1,3);

% GLP estimation, sample 2
resGLP2 = bvarGLP([y(T02:T12,:)],lags,'mcmc',1,'MCMCconst',50,'MNpsi',0,'noc',0,'sur',0,'Ndraws',M);

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
sCF2 = sort(CF2,3);

% plot of the impulse responses of unemployment and the labor shares
qqq=[.025 .16 .5 .84 .975];         % percentiles of the posterior distribution
indfig=[1 9 10];                    % variable position of impulse responses to plot

figure('Position', [0, 0, 1050, 300]);
count=0;
for jn = indfig
    count=count+1;
    subplot(1,length(indfig),count)
    quantilePlot([0:H]', squeeze(sCF1(:,jn,round(qqq*ndraws)))); hold on; grid on;
    quantilePlot([0:H]', squeeze(sCF2(:,jn,round(qqq*ndraws))),[.8941, .1020, .1098]);
    line([0 H],[0 0],'color','k')
    xlabel('horizon')
    ylabel(YLABEL(jn))
    title(series(jn));
end
subplot(1,3,1); plot([0:H]', medU,'--','LineWidth',3,'color',[.2157, .4941, .7216]);
plots=get(gca, 'Children'); legend(plots([7 4]),{'pre 1990','post 1990'},'Location','southwest');