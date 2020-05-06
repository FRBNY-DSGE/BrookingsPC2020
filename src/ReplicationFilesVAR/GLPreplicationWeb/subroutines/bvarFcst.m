function forecast = bvarFcst(y,beta,hz);

% computes the forecasts for y
% at the horizons specifed in hz
% using the coefficients beta

[k,n] = size(beta);

lags = (k-1)/n;

T = size(y,1);

%%% computed the forecats 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Y=[y;zeros(hz(end),n)];
for tau=1:max(hz)
    xT=[1;reshape(Y([T+tau-1:-1:T+tau-lags],:)',k-1,1)]';
    Y(T+tau,:)=xT*beta;
end

forecast = Y(T+hz,:);
