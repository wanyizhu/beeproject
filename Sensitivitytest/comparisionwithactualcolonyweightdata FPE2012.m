Predictedmaxhiveweight=[Y1(225),Y1(233),Y1(241),Y1(248),Y1(255),Y1(259),Y1(269),Y1(279)];foundationweight = 50.2 * 453.6 /1000;
Simmaxcolonyweight=Predictedmaxhiveweight+foundationweight*ones(1,8);


Predictedmeanhiveweight=[Y1(225),Y1(233),Y1(241),Y1(248),Y1(255),Y1(259),Y1(269),Y1(279)];foundationweight = 50.2 * 453.6 /1000;
Simmeancolonyweight=Predictedmeanhiveweight+foundationweight*ones(1,8);

Predictedminhiveweight=[Y1(225),Y1(233),Y1(241),Y1(248),Y1(255),Y1(259),Y1(269),Y1(279)];foundationweight = 50.2 * 453.6 /1000;
Simmincolonyweight=Predictedminhiveweight+foundationweight*ones(1,8);


Actualweight=[45536.4,44352,45588.614,43732.483,42968.621,39974.861,38795.501,41107.046]/1000;
figure 
startDate=datenum('07-16-2012');
endDate=datenum('09-10-2012');

xData=linspace(startDate,endDate,8);
plot(xData,realmeanweight,':ks')
set(gca,'XTick',xData)
datetick('x','mm/dd','keepticks')
hold on 
plot(xData,realmaxweight,':^k')
plot(xData,Simmeancolonyweight,'-or')
plot(xData,realminweight,':vk')


legend('Realmeanweight','Realmaxweight','SimmeanWeight','Realminweight')
xlabel('Date of 2012')
ylabel('Colony Weight (kg)')
title('Dynamics of a healthy honey bee colony weight')

%%%%%%%%%%%5Dimethoate treatment%%%%%%%%%%%%%
startDate=datenum('08-02-2012');
endDate=datenum('09-14-2012');

Predictedmeanhiveweight=[Y1(241),Y1(248),Y1(255),Y1(259),Y1(269),Y1(279),Y1(283)];foundationweight = 50.2 * 453.6 /1000;
Simmeancolonyweight=Predictedmeanhiveweight+foundationweight*ones(1,7);

xData=linspace(startDate,endDate,7);
plot(xData,dimethoatemaxcolonyweight,':^k')
set(gca,'XTick',xData)
datetick('x','mm/dd','keepticks')
hold on 
plot(xData,dimethoatemeancolonyweight,':ks')
plot(xData,Simmeancolonyweight,'-or')
plot(xData,dimethoatemincolonyweight,':vk')


legend('Dimmaxweight','Dimmeanweight','Simmeanweight','DimminWeight')
xlabel('Date of 2012')
ylabel('Colony Weight (kg)')
title('Dynamics of a dimethoate-treated honey bee colony weight')




%%%%%%%%%%%%%%%%%%%%%Fungicide Wax Treatment%%%%%%%%%%%%%%

startDate=datenum('08-02-2012');
endDate=datenum('09-14-2012');

Predictedmeanhiveweight=[Y1(241),Y1(248),Y1(255),Y1(259),Y1(269),Y1(279),Y1(283)];foundationweight = 50.2 * 453.6 /1000;
Simmeancolonyweight=Predictedmeanhiveweight+foundationweight*ones(1,7);

xData=linspace(startDate,endDate,7);
plot(xData,FungicideWaxMaxcolonyweight,':ks')
set(gca,'XTick',xData)
datetick('x','mm/dd','keepticks')
hold on 
plot(xData,FungicideWaxMeancolonyweight,':^k')
plot(xData,Simmeancolonyweight,'-or')
plot(xData,FungicideWaxMincolonyweight,':vk')


legend('Fungmaxweight','Fungmeanweight','Simmeanweight','FungminWeight')
xlabel('Date of 2012')
ylabel('Colony Weight (kg)')
title('Dynamics of a Fungicide-wax treated honey bee colony weight')




  figure;
Controlpollen=A;
plot(Control(1:360));hold on
plot(Pollendimethoate(1:360)) ; 
t=[0:30:360];months=['Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec';]; 
set(gca,'xtick',t)
set(gca,'xticklabel',months)
xlabel('Date')
ylabel('Pollen Weight')