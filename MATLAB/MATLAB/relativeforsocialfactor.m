
BRRS=[RSSF1(:,1), RSSF2(:,1),RSSF3(:,1),RSSW1(:,1),RSSW2(:,1),RSSW3(:,1)];
PRRS=[RSSF1(:,2), RSSF2(:,2),RSSF3(:,2),RSSW1(:,2),RSSW2(:,2),RSSW3(:,2)];
RTRS=[RSSF1(:,3), RSSF2(:,3),RSSF3(:,3),RSSW1(:,3),RSSW2(:,3),RSSW3(:,3)];
RDRS=[RSSF1(:,4), RSSF2(:,4),RSSF3(:,4),RSSW1(:,4),RSSW2(:,4),RSSW3(:,4)];
FRRS=[RSSF1(:,5), RSSF2(:,5),RSSF3(:,5),RSSW1(:,5),RSSW2(:,5),RSSW3(:,5)];


figure
%subplot(2,1,1)
linestyles = cellstr(char('-',':','-.','--','-',':','-.','--','-',':','-',':','-.','--','-',':','-.','--','-',':','-.'));
    
MarkerEdgeColors=Lines(6);  % n is the number of different items you have
Markers=['o','x','+','*','s','d','v','^','<','>','p','h','.','+','*','o','x','^','<','h','.','>','p','s','d','v','o','x','+','*','s','d','v','^','<','>','p','h','.'];% [...]

hold on;
%plot1 = plot(Survival,RSF1(1:20,1:6));

for i=1:1:6
plot(RTM.*8,BRRS(1:20,i),[linestyles{i} Markers(i)],'Color',MarkerEdgeColors(i,:));
end 
legend('First Nov','First End of Season','Second Nov','Second End of Season','Third Nov','Third End of Season')
Title('Relative Sensitivity of adult bees number to Brood Rearing Efficiency')
xlabel('Brood Rearing Efficiency');
ylabel('Relative Sensitivity'); 


%%%%%

% subplot(3,1,2)
linestyles = cellstr(char('-',':','-.','--','-',':','-.','--','-',':','-',':','-.','--','-',':','-.','--','-',':','-.'));
MarkerEdgeColors=Lines(6);  % n is the number of different items you have
Markers=['o','x','+','*','s','d','v','^','<','>','p','h','.','+','*','o','x','^','<','h','.','>','p','s','d','v','o','x','+','*','s','d','v','^','<','>','p','h','.'];% [...]
figure 
hold on;

%plot1 = plot(Survival,RSF1(1:20,1:6));
for i=1:1:6
plot(RTM,PRRS(1:20,i),[linestyles{i} Markers(i)],'Color',MarkerEdgeColors(i,:));
end 
legend('First Nov','First End of Season','Second Nov','Second End of Season','Third Nov','Third End of Season')
Title('Relative Sensitivity of adult bees number to Precocious behavior')
xlabel('Precocious Probability');
ylabel('Relative Sensitivity'); 



%%%%%%
%subplot(3,1,3)
linestyles = cellstr(char('-',':','-.','--','-',':','-.','--','-',':','-',':','-.','--','-',':','-.','--','-',':','-.'));
    
MarkerEdgeColors=Lines(6);  % n is the number of different items you have
Markers=['o','x','+','*','s','d','v','^','<','>','p','h','.','+','*','o','x','^','<','h','.','>','p','s','d','v','o','x','+','*','s','d','v','^','<','>','p','h','.'];% [...]
    
figure
hold on;

for i=1:1:6
plot(RTM,RTRS(1:20,i),[linestyles{i} Markers(i)],'Color',MarkerEdgeColors(i,:));
end 
legend('First Nov','First End of Season','Second Nov','Second End of Season','Third Nov','Third End of Season')
Title('Relative Sensitivity of adult bees number to Reversed Behavior')
xlabel('Reversed probability');
ylabel('Relative Sensitivity'); 



%%%%%

%subplot(3,1,1)
linestyles = cellstr(char('-',':','-.','--','-',':','-.','--','-',':','-',':','-.','--','-',':','-.','--','-',':','-.'));
    
MarkerEdgeColors=Lines(6);  % n is the number of different items you have
Markers=['o','x','+','*','s','d','v','^','<','>','p','h','.','+','*','o','x','^','<','h','.','>','p','s','d','v','o','x','+','*','s','d','v','^','<','>','p','h','.'];% [...]
    
figure
hold on;

for i=1:1:6
plot(RTM,RDRS(1:20,i),[linestyles{i} Markers(i)],'Color',MarkerEdgeColors(i,:));
end 
legend('First Nov','First End of Season','Second Nov','Second End of Season','Third Nov','Third End of Season')
Title('Relative Sensitivity of adult bees number to Retarded development')
xlabel('Retarded Development probability');
ylabel('Relative Sensitivity'); 


%%%%%
%subplot(3,1,2)
linestyles = cellstr(char('-',':','-.','--','-',':','-.','--','-',':','-',':','-.','--','-',':','-.','--','-',':','-.'));
    
MarkerEdgeColors=Lines(6);  % n is the number of different items you have
Markers=['o','x','+','*','s','d','v','^','<','>','p','h','.','+','*','o','x','^','<','h','.','>','p','s','d','v','o','x','+','*','s','d','v','^','<','>','p','h','.'];% [...]
    
figure
hold on;

for i=1:1:6
plot(RTM,FRRS(1:20,i),[linestyles{i} Markers(i)],'Color',MarkerEdgeColors(i,:));
end 
legend('First Nov','First End of Season','Second Nov','Second End of Season','Third Nov','Third End of Season')
Title('Relative Sensitivity of adult bees number to Foraging Efficiency')
xlabel('Foraging Efficiency');
ylabel('Relative Sensitivity'); 

%%%%%
%subplot(3,1,3)
% linestyles = cellstr(char('-',':','-.','--','-',':','-.','--','-',':','-',':','-.','--','-',':','-.','--','-',':','-.'));
%     
% MarkerEdgeColors=Lines(6);  % n is the number of different items you have
% Markers=['o','x','+','*','s','d','v','^','<','>','p','h','.','+','*','o','x','^','<','h','.','>','p','s','d','v','o','x','+','*','s','d','v','^','<','>','p','h','.'];% [...]
%     
% figure
% hold on;
% 
% for i=1:1:6
% plot(Survival,ForagerRS(1:20,i),[linestyles{i} Markers(i)],'Color',MarkerEdgeColors(i,:));
% end 
% legend('First Nov','First End of Season','Second Nov','Second End of Season','Third Nov','Third End of Season')
% Title('Relative Sensitivity of adult bees number to Forager survivorship')
% xlabel('Survival Rate');
% ylabel('Relative Sensitivity'); 
