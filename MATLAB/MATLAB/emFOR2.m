
agemax = 60; % +1 because of matlab indexing
global mt1 mt2 mt3 mt4 mt5 mt6 BR u v rt foragingsuccess Q nt ws; 
 %mt1=0.96; % 0.86--survivorship for egg stage 
 %%% Matrix Type%%%%%
 %BRM=[0 0.2 0.4 0.6 0.8 1 2 4 5 6 ];
 u=0.0;%precocious prob
 %RTM=zeros(1,10);
 EM=zeros(1,10);
 v=0;% reversed prob. between foragers and house bees;
% NTM=zeros(1,10);
 BR=5;
 Q=3;
 FRM=zeros(1,10);

% theta = rt*ones(agemax-1,1); % theta = probabilities of development retardation
rt=0;
 Honeystore1=zeros(1,10);
 AdultBees1=zeros(1,10);
  Honeystore2=zeros(1,10);
 AdultBees2=zeros(1,10);
% for i=1:1:10
mt2= (0.85)^(1/7); %---survivorship for larval stage 
% mt2= 0.87;
%mt3= (0.86)^(1/7); % 86%--survivorship for pupal stage 
mt3= (0.86)^(1/7);
          mt4= (0.85)^(1/16); % 0.99-85%--survivorship for nurse bee stage 
           %mt4=0.84; 
              mt5=(0.85)^(1/6); % 0.96-88.6%--survivorship for house bee stage 
              mt6=(0.78)^(1/12); % 78.5%--survivorship for forager bee stage
              % mt6=(0.9)^(1/12);
foragingsuccess=0.09;
nt=8;
ws=0.8;
% 
for i=1:1:10
    mt1=((i-1)/10)^(1/3);
    EM(i)=(i-1)/10;
% 
%     WM(i)=ws;
% %     nt=i;
% %     NTM(i)=nt; 
% % % %      mt5= (i/10); %---survivorship for larval stage 
% % % %      HouseM(i)=mt5;
% %      BR=BRM(i);
% foragingsuccess=(i-1)/50;
% FRM(i)=foragingsuccess; 
%rt=(i-1)/10;
% v=(i-1)/100; 
%mt2=((i-1)/10)^(1/7);
% mt1=((i-1)/10)^(1/3);
% EM(i)=(i-1)/10;
%    % FRM(i)=foragingsuccess;
%    % mt6=((i-1)/10)^(1/12);
% %    mt4=((i-1)/10)^(1/16);
% %    NM(i)= (i-1)/10; 
%    
   %RTM(i)=rt;
   %FM(i)=(i-1)/10; 
   

   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
G = zeros(6,agemax);
tx=240;
G(1,1:3)=1; G(2,4:11)=1; G(3,12:26)=1; G(4,27:42)=1;G(5,43:48)=1;G(6,49:agemax)=1;

P0 = 1000;

V0 = 300000 - P0;

H0=0;

R0=0;

N = zeros(agemax,1);

N(1:3)=100;

N(4:11)=100;

N(12:26)=300;

N(27:42)=200;

N(43:48)=200;

N(49:agemax)=100;

X = [ V0; P0; H0;R0; N ];

res=zeros(6,tx);

V=zeros(1,tx);

P=zeros(1,tx);

H=zeros(1,tx);

R=zeros(1,tx);


numyears =2;
summerdays = 240;
yeardays = 360;

pop=zeros(6,yeardays*numyears);
Vpop=zeros(1,yeardays*numyears);
Ppop=zeros(1,yeardays*numyears);
Hpop=zeros(1,yeardays*numyears);
Rpop=zeros(1,yeardays*numyears);

for T = 0:(numyears-1)

          for t=(yeardays*T+1):(yeardays*T+summerdays)

		     X = bees(X,t);

		      res(1:6,t-yeardays*T)=G*X(5:end);
 
		     V(1,t-yeardays*T)= X(1);

		     P(1,t-yeardays*T)= X(2);
        
             H(1,t-yeardays*T)= X(3);

		     R(1,t-yeardays*T)= X(4);
 
          end
     
%     pollen=P(1,yeardays*T+summerdays);
%     
%     Honey= H(1,yeardays*T+summerdays);
    
	pop(:,(yeardays*T+1):(yeardays*T+summerdays))=res;
    Vpop(:,(yeardays*T+1):(yeardays*T+summerdays))=V;
    Ppop(:,(yeardays*T+1):(yeardays*T+summerdays))=P;
    Hpop(:,(yeardays*T+1):(yeardays*T+summerdays))=H;
    Rpop(:,(yeardays*T+1):(yeardays*T+summerdays))=R;
	
%     if pop(4,(yeardays*T+1):(yeardays*T+summerdays))+pop(5,(yeardays*T+1):(yeardays*T+summerdays))>=1000
    
    % First Season Winter Dynamics 

	agemaxwinter=150; 

	W = zeros(4,agemaxwinter);

	W(1,1:3)=1; W(2,4:11)=1; W(3,12:26)=1; W(4,27:agemaxwinter)=1;

	N = zeros(agemaxwinter,1);

	N(1:3)=res(1,summerdays)/3;

	N(4:11)=res(2,summerdays)/8;

	N(12:26)=res(3,summerdays)/15;

	N(27:agemaxwinter)=(res(4,summerdays)+res(5,summerdays)+res(6,summerdays))/100;

	P0 = P(1,summerdays);

    V0 = V(1,summerdays);

    H0 = H(1,summerdays);

    R0 = R(1,summerdays);
    

    Y = [ V0; P0; H0;R0; N ];

	clear res V P H R;

	res=zeros(6,yeardays-summerdays);
    
    V=zeros(1,yeardays-summerdays);

    P=zeros(1,yeardays-summerdays);

    H=zeros(1,yeardays-summerdays);

    R=zeros(1,yeardays-summerdays);
    	

	for t=(yeardays*T+summerdays+1):(yeardays*(T+1))

		Y = winterbeesR(Y,t);

		res(1:4,(t-(yeardays*T+summerdays)))=W*Y(5:end);
      
        V(1,(t-(yeardays*T+summerdays)))= Y(1);

        P(1,(t-(yeardays*T+summerdays)))= Y(2);
    
        H(1,(t-(yeardays*T+summerdays)))= Y(3);

        R(1,(t-(yeardays*T+summerdays)))= Y(4);
        
     
    end 
    
	pop(:, (yeardays*T+summerdays+1):(yeardays*(T+1))) =res;
    
    Vpop(1,(yeardays*T+summerdays+1):(yeardays*(T+1)))= V;
    
    Ppop(1,(yeardays*T+summerdays+1):(yeardays*(T+1)))= P;
    
    Hpop (1,(yeardays*T+summerdays+1):(yeardays*(T+1)))= H;
    
    Rpop (1,(yeardays*T+summerdays+1):(yeardays*(T+1)))= R;
    
        
    
%     
%      
           
	%Second Season Summer Dynamics 

	N = zeros(agemax,1);

	N(1:3)=pop(1,yeardays*(T+1))/3;

	N(4:11)=pop(2,yeardays*(T+1))/8;

	N(12:26)=pop(3,yeardays*(T+1))/15;

	N(27:42)= pop(4,yeardays*(T+1))/34;

	N(43:48)= pop(4,yeardays*(T+1))/34 ;

	N(49:agemax)=pop(4,yeardays*(T+1))/34;

	P0 = Ppop(1,yeardays*(T+1));

	V0 = Vpop(1,yeardays*(T+1));

	R0= Rpop(1,yeardays*(T+1));
    
    H0= Hpop(1,yeardays*(T+1)); 

	X = [ V0; P0; H0; R0; N];

	res=zeros(6,summerdays);

	R=zeros(1,summerdays);

	V=zeros(1,summerdays);

	P=zeros(1,summerdays);

    H= zeros(1,summerdays); 
%     else 
%          fprintf('colony collapse');
%         break 
%        
%     end 
end
Honeystore1(i)=Hpop(241);
AdultBees1(i)=pop(4,241);
Honeystore2(i)=Hpop(601);
AdultBees2(i)=pop(4,601);
 end 


HoneyE1=Honeystore1;
HoneyE2=Honeystore2;
AdultE1=AdultBees1;
AdultE2=AdultBees2;

 %FM=ForagerM.^6;
 %VP=VM.^16; 
% plot3(PupalM,Honeystore,AdultBees)
% figure(9);
% stem3 (VM, Honeystore2, AdultBees2, 'DisplayName', 'BRM, Honeystore, AdultBees'); figure(gcf);
% xlabel('Precicous Behavior Probability');
% % xlabel('Nectar foraging time');
% % xlabel('Brood Rearing Efficiency');
% ylabel('Honeystore at November');
% zlabel('Adult bees Number at November');

% % HM=HouseM.^16;
% % subplot(1,2,1);
% %
% % hold on 

% % plot(NM,HoneystoreL,'--bs','LineWidth',2);
% % legend('Honeystores chaged with NurseBee Survival','HoneyStores chaged with Larval Survival')
% % subplot(1,2,2);
% % plot(NM,AdultBees,'-ro','LineWidth',2);
% % hold on 
% % plot(NM,AdultBeesL,'-bo','LineWidth',2);
% % legend('AdultBees chaged with pupal Survival','AdultBees chaged with NurseBee Survival')
% 
% 
% 
% % 

% plot(EM,AdultBeesFM1,'-yo','LineWidth',2);
% hold on 
% plot(EM,AdultBeesHM1,'-mo','LineWidth',2);
% plot(EM,AdultBeesNM1,'-ro','LineWidth',2);
% plot(EM,AdultBeesPM1,'-bo','LineWidth',2);
% plot(EM,AdultBeesLM1,'-go','LineWidth',2);
%  plot(HM,AdultBeesEM1,'-co','LineWidth',2);
% legend('Change with ForagerBee','Change with HouseBee','Change with NurseBee','Change with Pupa','Change with Larva','Change with Egg')
% xlabel('Survival Rate');
% ylabel('AdultBees Number at First November')
% 

plot(WM,AdultBees2,'--rs','LineWidth',2);
hold on 
plot(WM,AdultBees1,'--gs','LineWidth',2);


