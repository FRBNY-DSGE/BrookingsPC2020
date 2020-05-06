function irf = bvarIrfs(beta,sigma,nshock,hmax);


% computes IRFs using cholesky ordering
% to shock in position nshock
% up to hosizon hmax
% based on beta and sigma


[k,n] = size(beta);

lags = (k-1)/n;

%%% IRFs at the posterior mode
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cholVCM=chol(sigma)';
Y=zeros(lags+hmax,n);
in=lags;
vecshock=zeros(n,1); vecshock(nshock)=1;
for tau=1:hmax
    xT=[reshape(Y([in+tau-1:-1:in+tau-lags],:)',k-1,1)]';
    Y(in+tau,:)=xT*beta(2:end,:)+(tau==1)*(cholVCM*vecshock)';
end

irf = Y(in+1:end,:);