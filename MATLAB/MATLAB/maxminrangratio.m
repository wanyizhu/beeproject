
figure 
% subplot(1,3,1)
plot(BNRnop1(1,60:340)','r-')
hold on 
plot(BNRnop1(2,60:340)','r:')
plot(BNRnop2(1,:)','g-')
hold on 
plot(BNRnop2(2,:)','g:')
hold on 
plot(BNRnop3(1,:)','m-')
hold on 
plot(BNRnop3(2,:)','m:')

% subplot(1,3,3)
% hold on 
% % subplot(1,3,2)
% plot(BNRnop2')
% % subplot(1,3,3)
% plot(BNRnop3')


plot(BNRRange(1,1:420))
hold on 
plot(BNRRange(1,440:780))
plot(BNRRange(1,800:1080))