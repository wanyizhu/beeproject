

RBRSensitivityW1=((BRAdultBeesW1-NOPAdultBeesW1)./NOPAdultBeesW1)./((BRM*8-5)/5);
RBRSensitivityW2=((BRAdultBeesW2-NOPAdultBeesW2)./NOPAdultBeesW2)./((BRM*8-5)/5);
RBRSensitivityW3=((BRAdultBeesW3-NOPAdultBeesW3)./NOPAdultBeesW3)./((BRM*8-5)/5);
RBRSensitivityF1=((BRAdultBeesF1-NOPAdultBeesF1)./NOPAdultBeesF1)./((BRM*8-5)/5);
RBRSensitivityF2=((BRAdultBeesF2-NOPAdultBeesF2)./NOPAdultBeesF2)./((BRM*8-5)/5);
RBRSensitivityF3=((BRAdultBeesF3-NOPAdultBeesF3)./NOPAdultBeesF3)./((BRM*8-5)/5);



RPRSensitivityW1=((PRAdultBeesW1-NOPAdultBeesW1)./NOPAdultBeesW1)./((UM-0.00000001)/0.00000001);
RPRSensitivityW2=((PRAdultBeesW2-NOPAdultBeesW2)./NOPAdultBeesW2)./((UM-0.00000001)/0.00000001);
RPRSensitivityW3=((PRAdultBeesW3-NOPAdultBeesW3)./NOPAdultBeesW3)./((UM-0.00000001)/0.00000001);
RPRSensitivityF1=((PRAdultBeesF1-NOPAdultBeesF1)./NOPAdultBeesF1)./((UM-0.00000001)/0.00000001);
RPRSensitivityF2=((PRAdultBeesF2-NOPAdultBeesF2)./NOPAdultBeesF2)./((UM-0.00000001)/0.00000001);
RPRSensitivityF3=((PRAdultBeesF3-NOPAdultBeesF3)./NOPAdultBeesF3)./((UM-0.00000001)/0.00000001);


RRTSensitivityW1=((RTAdultBeesW1-NOPAdultBeesW1)./NOPAdultBeesW1)./((RTM-0.00000001)/0.00000001);
RRTSensitivityW2=((RTAdultBeesW2-NOPAdultBeesW2)./NOPAdultBeesW2)./((RTM-0.00000001)/0.00000001);
RRTSensitivityW3=((RTAdultBeesW3-NOPAdultBeesW3)./NOPAdultBeesW3)./((RTM-0.00000001)/0.00000001);
RRTSensitivityF1=((RTAdultBeesF1-NOPAdultBeesF1)./NOPAdultBeesF1)./((RTM-0.00000001)/0.00000001);
RRTSensitivityF2=((RTAdultBeesF2-NOPAdultBeesF2)./NOPAdultBeesF2)./((RTM-0.00000001)/0.00000001);
RRTSensitivityF3=((RTAdultBeesF3-NOPAdultBeesF3)./NOPAdultBeesF3)./((RTM-0.00000001)/0.00000001);


RRDSensitivityW1=((RDAdultBeesW1-NOPAdultBeesW1)./NOPAdultBeesW1)./((RD-0.00000001)/0.00000001);
RRDSensitivityW2=((RDAdultBeesW2-NOPAdultBeesW2)./NOPAdultBeesW2)./((RD-0.00000001)/0.00000001);
RRDSensitivityW3=((RDAdultBeesW3-NOPAdultBeesW3)./NOPAdultBeesW3)./((RD-0.00000001)/0.00000001);
RRDSensitivityF1=((RDAdultBeesF1-NOPAdultBeesF1)./NOPAdultBeesF1)./((RD-0.00000001)/0.00000001);
RRDSensitivityF2=((RDAdultBeesF2-NOPAdultBeesF2)./NOPAdultBeesF2)./((RD-0.00000001)/0.00000001);
RRDSensitivityF3=((RDAdultBeesF3-NOPAdultBeesF3)./NOPAdultBeesF3)./((RD-0.00000001)/0.00000001);



RFRSensitivityW1=((FRAdultBeesW1-NOPAdultBeesW1)./NOPAdultBeesW1)./((FRM-0.1)/0.1);
RFRSensitivityW2=((FRAdultBeesW2-NOPAdultBeesW2)./NOPAdultBeesW2)./((FRM-0.1)/0.1);
RFRSensitivityW3=((FRAdultBeesW3-NOPAdultBeesW3)./NOPAdultBeesW3)./((FRM-0.1)/0.1);
RFRSensitivityF1=((FRAdultBeesF1-NOPAdultBeesF1)./NOPAdultBeesF1)./((FRM-0.1)/0.1);
RFRSensitivityF2=((FRAdultBeesF2-NOPAdultBeesF2)./NOPAdultBeesF2)./((FRM-0.1)/0.1);
RFRSensitivityF3=((FRAdultBeesF3-NOPAdultBeesF3)./NOPAdultBeesF3)./((FRM-0.1)/0.1);


RSSF1=[RBRSensitivityF1; RPRSensitivityF1; RRTSensitivityF1; RRDSensitivityF1; RFRSensitivityF1]'; 
RSSF2=[RBRSensitivityF2; RPRSensitivityF2; RRTSensitivityF2; RRDSensitivityF2; RFRSensitivityF2]';
RSSF3=[RBRSensitivityF3; RPRSensitivityF3; RRTSensitivityF3; RRDSensitivityF3; RFRSensitivityF3]'; 
RSSW1=[RBRSensitivityW1; RPRSensitivityW1; RRTSensitivityW1; RRDSensitivityW1; RFRSensitivityW1]'; 
RSSW2=[RBRSensitivityW2; RPRSensitivityW2; RRTSensitivityW2; RRDSensitivityW2; RFRSensitivityW2]';
RSSW3=[RBRSensitivityW3; RPRSensitivityW3; RRTSensitivityW3; RRDSensitivityW3; RFRSensitivityW3]'; 


%LarvaSensitivity=[LarvaSensitivityF1 LarvaSensitivityW1 LarvaSensitivityF2 LarvaSensitivityW2 LarvaSensitivityF3 LarvaSensitivityW3];


% subplot(3,2,1)
% plot(REggSensitivityW3, 'DisplayName', 'REggSensitivityW3', 'YDataSource', 'REggSensitivityW3'); hold all; plot(RLarvaSensitivityW3, 'DisplayName', 'RLarvaSensitivityW3', 'YDataSource', 'RLarvaSensitivityW3'); plot(RNurseSensitivityW3, 'DisplayName', 'RNurseSensitivityW3', 'YDataSource', 'RNurseSensitivityW3'); plot(RPupaSensitivityW3, 'DisplayName', 'RPupaSensitivityW3', 'YDataSource', 'RPupaSensitivityW3'); plot(RHouseSensitivityW3, 'DisplayName', 'RHouseSensitivityW3', 'YDataSource', 'RHouseSensitivityW3'); plot(RForagerSensitivityW3, 'DisplayName', 'RForagerSensitivityW3', 'YDataSource', 'RForagerSensitivityW3');  hold off; figure(gcf)
% subplot(3,2)
% plot(REggSensitivityW2, 'DisplayName', 'REggSensitivityW2', 'YDataSource', 'REggSensitivityW3'); hold all; plot(RLarvaSensitivityW3, 'DisplayName', 'RLarvaSensitivityW3', 'YDataSource', 'RLarvaSensitivityW3'); plot(RNurseSensitivityW3, 'DisplayName', 'RNurseSensitivityW3', 'YDataSource', 'RNurseSensitivityW3'); plot(RPupaSensitivityW3, 'DisplayName', 'RPupaSensitivityW3', 'YDataSource', 'RPupaSensitivityW3'); plot(RHouseSensitivityW3, 'DisplayName', 'RHouseSensitivityW3', 'YDataSource', 'RHouseSensitivityW3'); plot(RForagerSensitivityW3, 'DisplayName', 'RForagerSensitivityW3', 'YDataSource', 'RForagerSensitivityW3');  hold off; figure(gcf)



