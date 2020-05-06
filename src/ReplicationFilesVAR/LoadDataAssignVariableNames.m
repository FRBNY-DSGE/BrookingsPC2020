% This code loads the data and names the columns of the data matrix
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[data,TXT,RAW]=xlsread('DataVAR');
Time=[1948:.25:2019.5];
data=data(:,2:end); 

unem=data(:,1);                 % unemployment rate
cpiinfl=data(:,2)*4;            % annualized CPI inflation
corecpiinfl=data(:,3)*4;        % annualized core CPI inflation
pceinfl=data(:,4);              % annualized PCE inflation
corepceinfl=data(:,5);          % annualized core PCE inflation
gdpinfl=data(:,6)*4;            % annualized GDP inflation
ls=100*log(data(:,7));          % log labor share (total economy)
lsnfbs=100*log(data(:,8));      % log labor share in the nonfarm business sector
lsnfcs=100*log(data(:,9));      % log labor share in the nonfinancial corporate sector
rcphnfbsinfl=data(:,10)*4;      % 4*growth rate of real compensation per hour in the nonfarm business sector  
aheinfl=data(:,11)*4;           % annualized wage inflation (average hourly earnings of nonsupervisotry and production workers)
nunem=data(:,12);               % netural rate of unemployment (CBO)
coe=log(data(:,13));            % log compensation of employees (total economy)
h=log(data(:,14));              % log hours worked (total economy)
cphnfbsinfl=data(:,15)*4;       % 4*growth rate of compensation per hour in the nonfarm business sector
gdp=log(data(:,16));            % log nominal GDP
eci=log(data(:,17));            % log emplyment cost index
csi=data(:,18);                 % Stock and Watson's CSI index
ltepi=data(:,19);               % annualized longt-term inflation expectations
y2y=data(:,20);                 % 2-year treasury rate
rgdp=log(data(:,21));           % log real GDP
rpotgdp=log(data(:,22));        % log potential real GDP
epr=100*log(data(:,23));        % log employment-population ratio
pcegoodsinfl=data(:,24);        % annualized PCE inflation (goods)    
pceservicesinfl=data(:,25);     % annualized PCE inflation (services)
corecpigoodsinfl=data(:,26);    % annualized CPI inflation (coore goods)
corecpiservicesinfl=data(:,27); % annualized CPI inflation (coore services)
cpifoodinfl=data(:,28);         % annualized CPI inflation (food)
cpishelterinfl=data(:,29);      % annualized CPI inflation (shelter)
hou=data(:,30);                 % hou
oilprice=400*log(data(:,31));   % log oil price
invinfl=4*data(:,32);           % 4*growth rate of investment price index
equipinfl=4*data(:,33);         % 4*growth rate of equipment price index
ffr=data(:,34);                 % federal funds rate
hnfbs=log(data(:,35));          % log hours in the nonfarm business sector
pop=data(:,36);                 % population
MPshocksRR=data(:,37);          % Romer and Romer monetary policy shocks
MPshocksHF=data(:,38);          % "high frequency" monetary policy shocks
GZspread=data(:,39);            % GZ spread
EBP=data(:,40);                 % excess bond premium
eciWSPI_nsa=log(data(:,41));    % log employment cost index (nsa)
eciWSPI=log(data(:,42));        % log employment cost index
MoodysSpread=data(:,43);        % Moody's AAA_BAA spread
stepi=data(:,46);               % 1-year ahead inflation expectations (SPF)
epi1=data(:,47);                % 1-quarter-ahear inflation expectations (SPF)
epi2=data(:,48);                % 2-quarter-ahead inflation expectations, lagged (SPF)
epi3=data(:,49);                % 4-quarter-ahead inflation expectations, lagged (SPF)
epimich=data(:,50);             % 1-year ahead inflation expectations (Michigan survey)


w=coe-h;                        % log compensation per hour (total economy)
winfl=400*(exp(w(2:end))-exp(w(1:end-1)))./exp(w(1:end-1));
winfl=[nan;winfl];              % annualized wage inflation (total economy)
eciWSPIinfl=400*(exp(eciWSPI(2:end))-exp(eciWSPI(1:end-1)))./exp(eciWSPI(1:end-1));
eciWSPIinfl=[nan;eciWSPIinfl];  % annualized wage inflation (emplyment cost index)
rgdpgap=100*(rgdp-rpotgdp);     % real GDP gap