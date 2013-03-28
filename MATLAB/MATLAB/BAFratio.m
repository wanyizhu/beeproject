figure1 = figure('Color',[1 1 1]);
%i=3;YMATRIX1=[maxBAF(180:240) minBAF(180:240) maxBAF(540:600) minBAF(540:600) maxBAF(900:960) minBAF(900:960)] ; 
% 
% maxvBAF=zeros(1080,1);minvBAF=zeros(1080,1);
% for i=1:1080
% maxvBAF(i)=max(TyvBAF(:,i));
% minvBAF(i)=min(TyvBAF(:,i));
% end
% 
% %YMATRIX2=[minBAF(180:240) minBAF(540:600) minBAF(900:960)] ; 
% YMATRIX1=[maxvBAF(540:600) minvBAF(540:600) maxvBAF(900:960) minvBAF(900:960)] ;
%YMATRIX1=[maxuBAF(540:600) minuBAF(540:600) maxuBAF(900:960) minuBAF(900:960)] ; 
 YMATRIX1=[maxuvBAF(540:600) minuvBAF(540:600) maxuvBAF(900:960) minuvBAF(900:960)] ;
axes1 = axes('Parent',figure1,...
    'XTickLabel',{'Sep','Oct','Nov','Dec'},...
    'XTick',[0 30 60 90]);
box('on');
hold('all');

% Create multiple lines using matrix input to plot
%subplot(1,2,1)
plot1 = plot(YMATRIX1,'Parent',axes1,'LineWidth',2.0);
% set(plot1(1),'LineStyle',':','Color',[0 0 1],...
%     'DisplayName','Max First season');
% set(plot1(2),'LineStyle',':','Color',[0 1 1],...
%     'DisplayName','Min First season');
set(plot1(1),'LineStyle','-','Color',[0 0 1],'DisplayName','Max Second season');
set(plot1(2),'LineStyle',':','Color',[0 1 1],'DisplayName','Min Second season');

set(plot1(3),'LineStyle','-.','Color',[1 0 0],...
    'DisplayName','Max Third season');
set(plot1(4),'LineStyle',':','Color',[1 0 1],...
    'DisplayName','Min Third season');
xlabel('Time of the season (in Months)');

% Create ylabel
ylabel('Brood & Forager vs. In-hive Bees Ratio');

% Create title
title('Range of the brood and Forager-Hivebees Ratio in a three-year healthy colony');
% %subplot(1,2,2);
% plot2 = plot(YMATRIX2,'Parent',axes1,'LineWidth',2.0);
% 
% set(plot2(1),'LineStyle',':','Color',[0 0 1],...
%     'DisplayName','Min First season');

figure
subplot(1,2,1)
scatter3(X_uvBAF(:,1),X_uvBAF(:,2),TyuvBAF(:,360))
hold on 
scatter3(X_uvBAF(:,1),X_uvBAF(:,2),TyuvBAF(:,720))
scatter3(X_uvBAF(:,1),X_uvBAF(:,2),TyuvBAF(:,1080))

subplot(1,2,2)
scatter3(X_uvBAF(:,1),X_uvBAF(:,2),TyuvBAF(:,240))
hold on 
scatter3(X_uvBAF(:,1),X_uvBAF(:,2),TyuvBAF(:,600))
scatter3(X_uvBAF(:,1),X_uvBAF(:,2),TyuvBAF(:,960))
