function [shatnew,signew,v,k,sigmainv]=kfilter_const(y,c,Z,G,C,T,H,shat,sig)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODEL
% y(t) = c + Z * s(t) + G*me(t)  ~ N(0,I)
% s(t) = C + T s(t-1) + H*eta(t) ~ N(0,I)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=length(y);

omega=T*sig*T'+H*H';
if G(end,1)>=10e10
    sigmainv=[eye(n-1)/(Z(1:end-1,:)*omega*Z(1:end-1,:)') zeros(n-1,1);zeros(1,n)];
else
    sigmainv=eye(n)/(Z*omega*Z'+G*G');
end
k=omega*Z'*sigmainv;
v=y-c-Z*(C+T*shat);
shatnew=C+T*shat+k*v;
signew=omega-k*Z*omega;
