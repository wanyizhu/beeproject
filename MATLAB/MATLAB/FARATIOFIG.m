%  X1:  vector of x data
%  YMATRIX1:  matrix of y data
%  X2:  vector of x data
%  YMATRIX2:  matrix of y data
%  X3:  vector of x data
%  YMATRIX3:  matrix of y data
FARatio=FARatio(1:1080);
DRFAE=(FAREM(13,:)-FARatio)./FARatio;
DRFAL=(FARLM(14,:)-FARatio)./FARatio;
DRFAP=(FARPM(15,:)-FARatio)./FARatio;
DRFAN=(FARNM(16,:)-FARatio)./FARatio;
DRFAH=(FARHM(16,:)-FARatio)./FARatio;
DRFAF=(FARFM(13,:)-FARatio)./FARatio;
%DRFAN=(BNRBR(5,:)-BNRatio)./BNRatio;

X1=[1:420];
YMatrix1=[DRFAE(1:420);DRFAL(1:420);DRFAP(1:420);DRFAN(1:420);DRFAH(1:420);DRFAF(1:420)];
X2=[440:720];
YMatrix2=[DRFAE(440:720);DRFAL(440:720);DRFAP(440:720);DRFAN(440:720);DRFAH(440:720);DRFAF(440:720)];
X3=[800:1080];
YMatrix3=[DRFAE(800:1080);DRFAL(800:1080);DRFAP(800:1080);DRFAN(800:1080);DRFAH(800:1080);DRFAF(800:1080)];
createfigureBNRatio(X1, YMatrix1, X2, YMatrix2, X3, YMatrix3)

%BNRMatrix=[DREB(1:1080);DRLB(1:1080);DRPB(1:1080);DRNB(1:1080);DRHB(1:1080);DRFB(1:1080)];
% % subplot(1,2,1)
% % plot(BNRatio(150:240),'rs')% August-November
% % % hold on 
% % % plot(BNREM(150:240),'g+')% August-November
% % % subplot(1,2,1)
% % hold on
% % plot(BNRatio(510:600),'m+')
% % plot(BNRatio(870:960),'b:')
% % 
% % % subplot(1,2,2)
% % % DREM=BNREM(11,:)-BNRatio;
% % % plot(BNREM(150:240));
% % % hold on 
% % % DRE=BNREM(10,:)-BNRatio;
% % % plot(DRE(150:240));
% % % hold on 
% % % DREb=BNREM(9,:)-BNRatio;
% % % plot(DREb(150:240));
% % % hold on 
% % % DREa=BNREM(12,:)-BNRatio;
% % % plot(DREa(150:240),'rs');
% % 
%  %subplot(1,2,1)
% % DRLA=(BNRLM(13,:)-BNRatio)./BNRatio;
% % plot(DRLA(1:360),'ks');
% % hold on 
%  %DRLB=(BNRLM(14,:)-BNRatio)./BNRatio;
% % plot(DRLB(1:400),'r*');
% % hold on 
% % % DRLA=(BNRLM(13,:)-BNRatio)./BNRatio;
% % % plot(DRLA(1:400),'g+');
% % % % DREa=(BNREM(11,:)-BNRatio)./BNRatio;
% % % % plot(DREa(1:360),'gs');
% % % hold on 
% 
% 
% DREb=(BNREM(13,:)-BNRatio)./BNRatio;
% plot(Time(1:420),DREb(1:420),'m+',Time(440:720),DREb(440:720),'m*',Time(800:1080),DREb(800:1080),'m*');
% %legend('Egg')
% hold on 
% DRLB=(BNRLM(14,:)-BNRatio)./BNRatio;
% plot(Time(1:420),DRLB(1:420),'r+', Time(440:720),DRLB(440:720),'r*',Time(800:1080),DRLB(800:1080),'r*');
% hold on 
% DRPB=(BNRPM(15,:)-BNRatio)./BNRatio;
% plot(Time(1:420),DRPB(1:420),'k+',Time(440:720),DRPB(440:720),'k*',Time(800:1080),DRPB(800:1080),'k*');
% hold on 
% DRNB=(BNRNM(16,:)-BNRatio)./BNRatio;
% plot(Time(1:420),DRNB(1:420),'g+',Time(440:720),DRNB(440:720),'g*',Time(800:1080),DRNB(800:1080),'g*');
% hold on 
% DRHB=(BNRHM(16,:)-BNRatio)./BNRatio;
% plot(Time(1:420),DRHB(1:420),'b+',Time(440:720),DRHB(440:720),'b*',Time(800:1080),DRHB(800:1080),'b*');
% hold on 
% DRFB=(BNRFM(13,:)-BNRatio)./BNRatio;
% plot(Time(1:420),DRFB(1:420),'c+',Time(440:720),DRFB(440:720),'c*',Time(800:1080),DRFB(800:1080),'c*');
% 
% %legend('Egg','Larva','Pupa','Nurse','House','Forager')
% % DRFBa=(BNRFM(12,:)-BNRatio)./BNRatio;
% % plot(Time(1:420),DRFBa(1:420),'y+',Time(440:720),DRFBa(440:720),'y*',Time(800:1080),DRFBa(800:1080),'y*');
% 
% 
% %DRHBa=(BNRHM(15,:)-BNRatio)./BNRatio;
% % plot(Time(1:420),DRHBa(1:420),'y+',Time(440:720),DRHBa(440:720),'y*',Time(800:1080),DRHBa(800:1080),'y*');
% % 
% %plot(Time(1:420),DRHBa(1:420),'y+',Time(440:720),DRHBa(440:720),'y*');
% 
% % hold on 
% % DRNa=(BNRNM(16,:)-BNRatio)./BNRatio;
% % plot(DRNa(440:720),'c+');
% %  subplot(1,2,2)
% % DRLA=(BNRLM(14,:)-BNRatio)./BNRatio;
% % plot(DRLA(510:600),'ks');
% % hold on 
% % DRLB=(BNRLM(15,:)-BNRatio)./BNRatio;
% % plot(DRLB(510:600),'k*');
% % hold on 
% % DREa=(BNREM(12,:)-BNRatio)./BNRatio;
% % plot(DREa(510:600),'rs');
% % hold on 
% % DREb=(BNREM(11,:)-BNRatio)./BNRatio;
% % plot(DREb(510:600),'r^');
% 
% 
% %figure
% % 
% % % plot(BNRatio(90:240),BNRatio(450:600),BNRatio(810:960))
% % %BARatio=(pop(4,1:360*numyears)+ pop(5,1:360*numyears)+pop(6,1:360*numyears))./(pop(1,1:360*numyears)+pop(2,1:360*numyears)+pop(3,1:360*numyears)); 
% % 
% % % figure
% % subplot(1,3,1)
% % plot(BNRatio(150:240),'rs')
% % hold on 
% % % plot(BNREM(11,150:240),'k+')
% % plot(BNREM(12,150:240),'gx')
% % % plot(BNREM(13,150:240),'m*')
% % % plot(BNRLM(11,150:240),'k+')
% % % plot(BNRLM(12,150:240),'gx')
% % plot(BNRLM(14,150:240),'m*')
% % 
% % 
% % subplot(1,3,2)
% % plot(BNRatio(510:600),'rs')
% % hold on 
% % % plot(BNREM(11,510:600),'k+')
% %  plot(BNREM(12,510:600),'gx')
% % % plot(BNREM(13,510:600),'m*')
% % % plot(BNRLM(11,510:600),'k+')
% % % plot(BNRLM(12,510:600),'gx')
% % plot(BNRLM(14,510:600),'m*')
% % 
% % subplot(1,3,3)
% % plot(BNRatio(870:960),'rs')
% % hold on 
% % % plot(BNREM(11,870:960),'k+')
% %  plot(BNREM(12,870:960),'gx')
% % % plot(BNREM(13,870:960),'m*')
% % % plot(BNRLM(11,870:960),'k+')
% % % plot(BNRLM(12,870:960),'gx')
% % plot(BNRLM(14,870:960),'m*')
% % 
% % % % figure
% % % % % subplot(1,3,1)
% % % % plot(ELPARatio(90:360),'rs')
% % % % hold on 
% % % % % plot(RatioEM(12,90:360),'g+')
% % % % % subplot(1,3,2)
% % % % plot(ELPARatio(450:720),'rs')
% % % % hold on 
% % % % plot(RatioEM(12,450:720),'g+')
% % % % subplot(1,3,3)
% % % %plot(ELPARatio(810:1080),'rs')
% % % % hold on 
% % % % plot(RatioEM(12,810:1080),'g+')
% % % 
% % % % figure 
% % % % subplot(1,3,1)
% % % % plot(BARatio(90:360),'rs')
% % % % hold on 
% % % % plot(RatioEM(13,90:360),'g+')
% % % % subplot(1,3,2)
% % % % plot(BARatio(450:720),'rs')
% % % % hold on 
% % % % plot(RatioEM(13,450:720),'g+')
% % % subplot(1,3,3)
% % % plot(BARatio(810:1080),'rs')
% % % hold on 
% % % plot(RatioEM(13,810:1080),'g+')
% % % 
% % % 
% % % 
% % % figure 
% % % subplot(1,3,1)
% % % plot(BARatio(90:360),'rs')
% % % hold on 
% % % plot(RatioEM(14,90:360),'g+')
% % % subplot(1,3,2)
% % % plot(BARatio(450:720),'rs')
% % % hold on 
% % % plot(RatioEM(14,450:720),'g+')
% % % subplot(1,3,3)
% % % plot(BARatio(810:1080),'rs')
% % % hold on 
% % % plot(RatioEM(14,810:1080),'g+')
% % % 
% % % 
% % % 
% % % figure 
% % % subplot(1,3,1)
% % % plot(BARatio(90:360),'rs')
% % % hold on 
% % % plot(RatioEM(11,90:360),'g+')
% % % subplot(1,3,2)
% % % plot(BARatio(450:720),'rs')
% % % hold on 
% % % plot(RatioEM(11,450:720),'g+')
% % % subplot(1,3,3)
% % % plot(BARatio(810:1080),'rs')
% % % hold on 
% % % plot(RatioEM(11,810:1080),'g+')
% % % 
% % % 
% % % 
% % % 
% % % 
