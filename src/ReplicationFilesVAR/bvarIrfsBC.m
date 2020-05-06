function irf = bvarIrfsBC(beta,sigma,hmax,nvar,f1,f2);

% computes IRFs to the shock that maximizes the variance of variable nvar
% in the [f1,f2] frequency band, up to hosizon hmax, given beta and sigma
% This code follows the notes VARspectral

[k,n] = size(beta);
lags = (k-1)/n;

cholVCM=chol(sigma)';

int=abs(f1-f2)/200;
V=0;
for f=f1:int:f2
    eminusif(1,1,:)=exp(-i*f.*[1:lags]);
    sumB=sum(reshape(beta(2:end,:)',n,n,lags).*repmat(eminusif,n,n,1),3);
    invA=(eye(n)-sumB)\cholVCM;
    V=V+real(invA(nvar,:)'*invA(nvar,:))*int/abs(f1-f2); %volatility of variable nvar over frequency band, expressed in terms of the contributions of all the Cholesky-transformed residuals
end
[q junk]=eigs(V,1); 
q=q*sign(q(nvar));


Y=zeros(lags+hmax,n);
in=lags;
vecshock=q;
for tau=1:hmax
    xT=[reshape(Y([in+tau-1:-1:in+tau-lags],:)',k-1,1)]';
    Y(in+tau,:)=xT*beta(2:end,:)+(tau==1)*(cholVCM*vecshock)';
end
irf = Y(in+1:end,:);





% function r=variance_f1f2(q,n,lags,beta,cholVCM,nvar,f1,f2)
% q=real(q); q=q./(q'*q);
% 
% int=abs(f1-f2)/100;
% V=0;
% for f=f1:int:f2
%     eminusif(1,1,:)=exp(-i*f.*[1:lags]);
%     sumB=sum(reshape(beta(2:end,:)',n,n,lags).*repmat(eminusif,n,n,1),3);
%     invA=(eye(n)-sumB)\cholVCM;
%     V=V+invA(nvar,:)'*invA(nvar,:)*int;
% end
% 
% r=-q'*V*q;











% [k,n] = size(beta);
% lags = (k-1)/n;
% 
% cholVCM=chol(sigma)';
% 
% 
% Y=zeros(n,n,lags+hmax);
% in=lags;
% 
% for tau=1:hmax       
%     Y(:,:,in+tau)=beta(2:end,:)'*reshape(permute(Y(:,:,in+tau-1:-1:in+tau-lags),[2 1 3]),n,n*lags)'+(tau==1)*cholVCM;
% end
% 
% irf = Y(:,:,in+1:end);
% 
% int=(max([f1,f2])-min([f1,f2]))/1000;
% 
% OMEGA=0;
% OMEGA1=0;
% 
% eio=nan(1,1,hmax);
% for f=f1:int:f2
%     eio(1,1,:)=exp(-i*f.*[0:hmax-1]);
%     eioplus(1,1,:)=exp(i*f.*[0:hmax-1]);
%     CC=sum(repmat(eio,n,n,1).*irf,3);
%     CCplus=sum(repmat(eioplus,n,n,1).*irf,3);
%     %OMEGA=OMEGA+CCplus(:,nvar)*CC(nvar,:)*int;
%     OMEGA=OMEGA+CC(nvar,:)'*CC(nvar,:)*int;
%     
%     
%     eio1=nan(1,1,lags);
%     eio1(1,1,:)=exp(-i*f.*[1:lags]);
% 
% 
%     BB=sum(reshape(beta(2:end,:)',n,n,lags).*repmat(eio1,n,n,1),3);
% 
% 
%     invA=(eye(n)-BB)\cholVCM;
%     
%     OMEGA1=OMEGA1+ invA(nvar,:)'*invA(nvar,:)*int;
%     
%     
%     
%     
% end
% 
% OMEGA






