% set up dates 

years= [ones(1,12)*2010,ones(1,12)*2011,ones(1,12)*2012];

months=[1:12,1:12,1:12];
days= ones(1,36);

for k= 1:36
    dates(k)=datenum(years(k),months(k),days(k));
end 

%plot data 
plot(100+20*rand(1,36),'o')

%Set ticks 
labels=datestr(dates,12);
set(gca,'XTick',dates);
set(gca,'XTickLabel',labels);

%Label Axes and set Title 
xlabel('Date')
ylabel('Colony Weight')

%months=['Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec';];
months=['Jan';'Mar';'May';'Jul';'Sep';'Nov';];

% axis annotation to months 
plot(Y1);
t=[0:30:360*3];
set(gca,'xtick',t)
set(gca,'xticklabel',months)

% axis annotation to months 
plot(Y1);
t=[0:60:360*3];
set(gca,'xtick',t)
set(gca,'xticklabel',months)

