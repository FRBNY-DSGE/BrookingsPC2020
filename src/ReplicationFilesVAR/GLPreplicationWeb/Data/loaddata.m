clear;

[DATA,TXT] = xlsread('DataSW.xls','Quarterly');

TimeQ = datenum(TXT(4:end,1),'mm\dd\yyyy');

YQ = DATA;

LongDescrQ = TXT(2,2:end);

ShortDescrQ = TXT(3,2:end);

%1    'RGDP': 'Real Gross Domestic Product, Quantity Index (2000=100) , SAAR'
%2    'PGDP': 'Gross Domestic Product, Price Index (2000=100) , SAAR'
%3    'Cons': 'Real Personal Consumption Expenditures, Quantity Index (2000=100) , SAAR'
%4    'GPDInv': 'Real Gross Private Domestic Investment, Quantity Index (2000=100) , SAAR'
%5    'Emp. Hours': 'HOURS OF ALL PERSONS: NONFARM BUSINESS SEC (1982=100,SA)'
%6    'Real Comp/Hour':  'REAL COMPENSATION PER HOUR,EMPLOYEES:NONFARM BUSINESS(82=100,SA)'

    

    
    
[DATA,TXT] = xlsread('DataSW.xls','Monthly');

TimeM = datenum(TXT(4:end,1),'mm\dd\yyyy');


YM = DATA;

LongDescrM = TXT(2,2:end);

ShortDescrM = TXT(3,2:end);

% 'FedFunds': 'INTEREST RATE: FEDERAL FUNDS (EFFECTIVE) (% PER ANNUM,NSA)'


if TimeM(2)~=TimeQ(1)
    error('The starting date of Monthly and Quarterly Data does not coincide')
end;


if TimeM(end-1)~=TimeQ(end)
    error('The ending date of Monthly and Quarterly Data does not coincide')
end;



YMQ = filter([1 1 1]/3,1,YM); YMQ = YMQ(3:3:end,:); % makes quarterly by taking the average 

Y = [YQ YMQ];
Time = TimeQ;
LongDescr  = [LongDescrQ  LongDescrM];
ShortDescr = [ShortDescrQ ShortDescrM];
Dates = datevec(Time);

y = [log(YQ)*400 YMQ]/100; 


%1    'RGDP': 'Real Gross Domestic Product, Quantity Index (2000=100) , SAAR'
%2    'PGDP': 'Gross Domestic Product, Price Index (2000=100) , SAAR'
%3    'Cons': 'Real Personal Consumption Expenditures, Quantity Index (2000=100) , SAAR'
%4    'GPDInv': 'Real Gross Private Domestic Investment, Quantity Index (2000=100) , SAAR'
%5    'Emp. Hours': 'HOURS OF ALL PERSONS: NONFARM BUSINESS SEC (1982=100,SA)'
%6    'Real Comp/Hour':  'REAL COMPENSATION PER HOUR,EMPLOYEES:NONFARM BUSINESS(82=100,SA)'
%7    'FedFunds': 'INTEREST RATE: FEDERAL FUNDS (EFFECTIVE) (% PER ANNUM,NSA)'

    


save DataSW y Time LongDescr ShortDescr Dates
