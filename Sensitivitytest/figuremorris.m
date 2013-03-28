figure1 = figure('Color',[1 1 1]);
i=50;
%YMatrix1=[Ty(i,80:240);Ty(i,440:600);Ty(i,800:960)]'; 
YMatrix2=[minR2;minR3;maxR3]'; 
axes1 = axes('Parent',figure1,...
    'XTickLabel',{'June','July','Aug','Sep','Oct','Nov','Dec'},...
    'XTick',[0 30 60 90 120 150 180]);
box('on');
hold('all');

% Create multiple lines using matrix input to plot
plot1 = plot(YMatrix2,'Parent',axes1,'LineWidth',2.0);
set(plot1(1),'LineStyle',':','Color',[0 0 1],...
    'DisplayName','First season');
set(plot1(2),'Color',[0 0 1],'DisplayName','Second season');
set(plot1(3),'LineStyle',':','Color',[0 1 0],...
    'DisplayName','Third season');


%YMatrix2=[maxR;minR;maxR2]'; 