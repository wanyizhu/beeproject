% plot(BARatio(1,60:240),'r:')
% hold on 
% plot(BARatio(1,780:960),'g*')
% hold on 
% plot(BARatio(1,420:600),'ks')
% hold on 
% plot(REM(60:240,11),'mx');
% hold on 
% plot(REM(420:600,11),'c+');


Time=[1:1080];
plot(Time(60:360),BARatio(1,60:360),'r:',Time(420:720),BARatio(1,420:720),'ks',Time(780:1080),BARatio(1,780:1080),'g*');

figure 
plot(Time(90:360),BNRatio(1,90:360),'r:',Time(450:720),BNRatio(1,450:720),'ks',Time(810:1080),BNRatio(1,810:1080),'g*');

figure
plot(Time(60:360),BNRatio(1,60:360),'r:',Time(420:720),BNRatio(1,420:720),'ks',Time(780:1080),BNRatio(1,780:1080),'g*');
