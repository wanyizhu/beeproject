% % 
% % 
% figure 
% plot(BNRRange(1,100:360),'b:')
% 
% hold on 
% plot(BNRRange(2,100:360),'b-')
% 
% plot(BNRRange(1,460:720),'c:')
% 
% plot(BNRRange(2,460:720),'c-')
% 
% plot(BNRRange(1,820:1080),'r:')
% 
% hold on 
% plot(BNRRange(2,820:1080),'r-')
% 
% % figure 
% % title('Forager vs. Adult Ratio at the field season')
% % plot(FARange(1,80:240),'b:')
% % hold on 
% % plot(FARange(2,80:240),'b-')
% % plot(FARange(1,440:600),'g:')
% % plot(FARange(2,440:600),'g-')
% % plot(FARange(1,800:960),'r:')
% % plot(FARange(2,800:960),'r-')
% 
% legend('Max-First season','Min-First season','Max-Second season','Min-Second season','Max-Third season','Min-Third season')
% 
% 
FAYMATRIX=[FARange(1,80:240);FARange(2,80:240);FARange(1,440:600);FARange(2,440:600);FARange(1,800:960);FARange(2,800:960)]'; 
FNCREATFIGURE(FAYMATRIX)

%FARange(1,180:240)+BNRange(1:180:240)

BNMATRIX=[BNRRange(1,80:240);BNRRange(2,80:240);BNRRange(1,440:600);BNRRange(2,440:600);BNRRange(1,800:960);BNRRange(2,800:960)]';

BFAF1=FARange(1,80:240)+BNRRange(1,80:240);
BFAF2=FARange(2,80:240)+BNRRange(2,80:240);
BFAF=[BFAF1;BFAF2];
BFAS1=FARange(1,440:600)+BNRRange(1,440:600);
BFAS2=FARange(2,440:600)+BNRRange(2,440:600);
BFAS=[BFAS1;BFAS2];
BFAT1=FARange(1,800:960)+BNRRange(1,800:960);
BFAT2=FARange(2,800:960)+BNRRange(2,800:960);
BFAT=[BFAT1;BFAT2];



