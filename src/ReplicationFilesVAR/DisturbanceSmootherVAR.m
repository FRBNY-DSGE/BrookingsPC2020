function [sdraw]=DisturbanceSmoother(y,c,Z,G,C,B,H,s00,P00,T,n,ns,ne,SS);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function performs a draws from the posterior of the disturbances
% and unobservable states of the following state-space model
%
% y(t) = c + Z * s(t) + G * me(t)  ~ N(0,I)
% s(t) = C + T s(t-1) + H * eta(t) ~ N(0,I)
% s(0) ~ N(s00,P00)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Kalman filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V=zeros(n,T);
K=zeros(ns,n,T);
HINV=zeros(n,n,T);
SHAT=zeros(ns,T);
SIG=zeros(ns,ns,T);
shat=s00;
sig=P00;

for t=1:T
    [shat,sig,v,k,hinv]=kfilter_forDS_VAR(y(t,:)',c,Z,G,C,B,H,shat,sig);
    SHAT(:,t)=shat; SIG(:,:,t)=sig;
    V(:,t)=v; K(:,:,t)=k; HINV(:,:,t)=hinv;
end

% disturbance smoother
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epshat=zeros(ne,T);
r=zeros(ns,1);
for t=T:-1:1
    epshat(:,t)=H'*Z'*squeeze(HINV(:,:,t))*V(:,t)+H'*(eye(ns)-squeeze(K(:,:,t))*Z)'*r;
    r=B'*Z'*squeeze(HINV(:,:,t))*V(:,t)+B'*(eye(ns)-squeeze(K(:,:,t))*Z)'*r;
end

if strcmp(SS,'smoother')
    % smoothed states
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sdraw=zeros(ns,T);
    sdraw(:,1)=C + B*s00 + H*epshat(:,1);
    for t=2:T
        sdraw(:,t)=C + B*sdraw(:,t-1) + H*epshat(:,t);
    end
    
elseif strcmp(SS,'simulation')
    % simulating new shocks, states and observables
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    epsplus=randn(ne,T);
    splus=zeros(ns,T);
    yplus=zeros(n,T);
    splus(:,1)=C + B*s00 + H*epsplus(:,1);
    yplus(:,1)=Z*splus(:,1);
    for t=2:T
        splus(:,t)=C + B*splus(:,t-1) + H*epsplus(:,t);
        yplus(:,t)=Z*splus(:,t);
    end
    
    % Kalman filter
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Vplus=zeros(n,T);
    Kplus=zeros(ns,n,T);
    HINVplus=zeros(n,n,T);
    SHATplus=zeros(ns,T);
    SIGplus=zeros(ns,ns,T);
    shat=s00;
    sig=P00;
    
    for t=1:T
        [shat,sig,v,k,hinv]=kfilter_forDS_VAR(yplus(:,t),c,Z,G,C,B,H,shat,sig);
        SHATplus(:,t)=shat; SIGplus(:,:,t)=sig;
        Vplus(:,t)=v; Kplus(:,:,t)=k; HINVplus(:,:,t)=hinv;
    end
    
    % disturbance smoother
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    epshatplus=zeros(ne,T);
    r=zeros(ns,1);
    for t=T:-1:1
        epshatplus(:,t)=H'*Z'*squeeze(HINVplus(:,:,t))*Vplus(:,t)+H'*(eye(ns)-squeeze(Kplus(:,:,t))*Z)'*r;
        r=B'*Z'*squeeze(HINVplus(:,:,t))*Vplus(:,t)+B'*(eye(ns)-squeeze(Kplus(:,:,t))*Z)'*r;
    end
    
    epsdraw=epshat+epsplus-epshatplus;
    
    sdraw=zeros(ns,T);
    sdraw(:,1)=C + B*s00 + H*epsdraw(:,1);
    for t=2:T
        sdraw(:,t)=C + B*sdraw(:,t-1) + H*epsdraw(:,t);
    end
end
