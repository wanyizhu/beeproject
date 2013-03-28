%nrm = [1/3,1/8,1/15,1/16,1/6,1/12]'; % normalize frequencies for each
%stage

% First Season Summer Dynamics 

agemax = 60; % +1 because of matlab indexing
global mt1 mt2 mt3 mt4 mt5 mt6 BR u v rt foragingsuccess Q nt; 
 global hsurfX hsurfY hsurf;
 %BRM=[0 0.2 0.4 0.6 0.8 1 2 4 5 6 ];
 u=0.0;%precocious prob
%  %RTM=zeros(1,10);
%  % EM=zeros(1,10);
v=0;% reversed prob. between foragers and house bees;
% NTM=zeros(1,10);
 BR=5;
 Q=3;

% theta = rt*ones(agemax-1,1); % theta = probabilities of development retardation
rt=0.;
%  Honeystore1=zeros(10,10);
%  AdultBees1=zeros(10,10);
%   Honeystore2=zeros(10,10);
%  AdultBees2=zeros(10,10);
 
 HoneystoreF1=zeros(1,20);% First November-F1
 AdultBeesF1=zeros(1,20);

HoneystoreW1=zeros(1,20);% First End of One season-Next Feb-W1
AdultBeesW1=zeros(1,20);

HoneystoreF2=zeros(1,20);% Second November-F2
AdultBeesF2=zeros(1,20);

HoneystoreW2=zeros(1,20);% Second End of One season-Second Feb-W2
AdultBeesW2=zeros(1,20);


HoneystoreF3=zeros(1,20);% Third November-F3
AdultBeesF3=zeros(1,20);

HoneystoreW3=zeros(1,20);% Third End of One season-Second Feb-W2
AdultBeesW3=zeros(1,20);
%  for i=1:1:10
mt1=(0.96)^(1/3); % 0.86--survivorship for egg stage 
mt2= (0.85)^(1/7); %---survivorship for larval stage 
mt3= (0.86)^(1/14);
mt4= (0.85)^(1/16); % 0.99-85%--survivorship for nurse bee stage 
mt5=(0.85)^(1/6); % 0.96-88.6%--survivorship for house bee stage 
mt6=(0.78)^(1/12); % 78.5%--survivorship for forager bee stage

foragingsuccess=0.1;
nt=8;


% for i=1:1:10
%    for j=1:1:10
%        mt1=((i-1)/10)^(1/3);
%        EM(i)=(i-1)/10;
%        mt2=((j-1)/10)^(1/7);
%        LM(j)=(j-1)/10; 
%         u=(i-1)/10;
%         UM(i)=u;
%    foragingsuccess=(j-1)/10;
%    FGM(j)=foragingsuccess; 
%       rt=(j-1)/10;
%       RTM(j)=rt;
       
%        u=(i-1)/10;
%        UM(i)=u;
%        v=(j-1)/10;
%        VM(j)=v;
%  mt1=((i-1)/10)^(1/3);
%  EM(i)=mt1;
% %     mt2=((j-1)/10)^(1/7);
% %     LM(j)=mt2;
%     
% %    mt4=((j-1)/10)^(1/16);
% %    NM(j)=mt4;
%       mt5=((j-1)/10)^(1/6);
%       HM(j)= mt5;
      
%       mt5=((i-1)/10)^(1/6);
%       HM(i)= mt5;
%       mt6=((j-1)/10)^(1/12);
%       FM(j)= mt6;
      
%     ws=(i-1)/10;
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


numyears =3;
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
% 	
%     if pop(4,(yeardays*T+1):(yeardays*T+summerdays))+pop(5,(yeardays*T+1):(yeardays*T+summerdays))<1000
%         
%           
%          fprintf('colony collapse');
%         break 
%     else 
%  
    
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
%     end    
    end 
% Honeystore1=Hpop(241);
% AdultBees1=pop(4,241);
% Honeystore2=Hpop(481);
% AdultBees2=pop(4,481);
% 
% HoneystoreF1(i,j)=Hpop(241);% First November-F1
% AdultBeesF1(i,j)=pop(4,241);
% 
% HoneystoreW1(i,j)=Hpop(360);% First End of One season-Next Feb-W1
% AdultBeesW1(i,j)=pop(4,360);
% 
% HoneystoreF2(i,j)=Hpop(601);% Second November-F2
% AdultBeesF2(i,j)=pop(4,601);
% 
% HoneystoreW2(i,j)=Hpop(720);% Second End of One season-Second Feb-W2
% AdultBeesW2(i,j)=pop(4,720);
% 
% 
% HoneystoreF3(i,j)=Hpop(961);% Third November-F3
% AdultBeesF3(i,j)=pop(4,961);
% 
% HoneystoreW3(i,j)=Hpop(1080);% Third End of One season-Second Feb-W2
% AdultBeesW3(i,j)=pop(4,1080);
% 
%    end 
%    end

% save UVAdultBees1.dat AdultBees1 -ascii
% save UVAdultBees2.dat AdultBees2 -ascii
% save UVHoneystore1.dat Honeystore1 -ascii
% save UVHoneystore2.dat Honeystore2 -ascii
% load UVAdultBees1.dat
% load UVAdultBees2.dat
% load UVHoneystore1.dat
% load UVHoneystore2.dat
% 
% save URTAdultBees1.dat AdultBees1 -ascii
% save URTAdultBees2.dat AdultBees2 -ascii
% save URTHoneystore1.dat Honeystore1 -ascii
% save URTHoneystore2.dat Honeystore2 -ascii
% load URTAdultBees1.dat
% load URTAdultBees2.dat
% load URTHoneystore1.dat
% load URTHoneystore2.dat

% save UFGAdultBees1.dat AdultBees1 -ascii
% save UFGAdultBees2.dat AdultBees2 -ascii
% save UFGHoneystore1.dat Honeystore1 -ascii
% save UFGHoneystore2.dat Honeystore2 -ascii
% load UFGAdultBees1.dat
% load UFGAdultBees2.dat
% load UFGHoneystore1.dat
% load UFGHoneystore2.dat

% save ELAdultBeesF1.dat AdultBeesF1 -ascii
% save ELAdultBeesF2.dat AdultBeesF2 -ascii
% save ELAdultBeesF3.dat AdultBeesF3 -ascii
% 
% save ELAdultBeesW1.dat AdultBeesW1 -ascii
% save ELAdultBeesW2.dat AdultBeesW2 -ascii
% save ELAdultBeesW3.dat AdultBeesW3 -ascii
% 
% 
% save ELHoneystoreW1.dat HoneystoreW1 -ascii
% save ELHoneystoreW2.dat HoneystoreW2 -ascii
% save ELHoneystoreW3.dat HoneystoreW3 -ascii
% 
% save ELHoneystoreF1.dat HoneystoreW1 -ascii
% save ELHoneystoreF2.dat HoneystoreW2 -ascii
% save ELHoneystoreF3.dat HoneystoreW3 -ascii
% 
% load ELAdultBeesF1.dat
% load ELAdultBeesF2.dat
% load ELAdultBeesF3.dat
% load ELAdultBeesW1.dat
% load ELAdultBeesW2.dat
% load ELAdultBeesW3.dat
% 
% 
% load ELHoneystoreW1.dat
% load ELHoneystoreW2.dat
% load ELHoneystoreW3.dat
% load ELHoneystoreF1.dat
% load ELHoneystoreF2.dat
% load ELHoneystoreF3.dat


% save EHAdultBees1.dat AdultBees1 -ascii
% save EHAdultBees2.dat AdultBees2 -ascii
% save EHHoneystore1.dat Honeystore1 -ascii
% save EHHoneystore2.dat Honeystore2 -ascii
% load EHAdultBees1.dat
% load EHAdultBees2.dat
% load EHHoneystore1.dat
% load EHHoneystore2.dat


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
% 
% plot(EM,AdultBeesFM2,'--ys','LineWidth',2);
% hold on 
% plot(EM,AdultBeesHM2,'--ms','LineWidth',2);
% plot(EM,AdultBeesNM2,'--rs','LineWidth',2);
% plot(EM,AdultBeesPM2,'--bs','LineWidth',2);
% plot(EM,AdultBeesLM2,'--gs','LineWidth',2);
%  plot(HM,AdultBeesEM2,'--cs','LineWidth',2);
% legend('Change with ForagerBee','Change with HouseBee','Change with NurseBee','Change with Pupa','Change with Larva','Change with Egg')
% xlabel('Survival Rate');
% ylabel('AdultBees Number at Second November')






% xlabel('Survival Rate');
% ylabel('HoneyStore at November')
% plot(HM,HoneystoreF,'--ys','LineWidth',2);
% hold on 
% plot(HM,HoneystoreH,'--ms','LineWidth',2);
% plot(HM,HoneystoreN,'--rs','LineWidth',2);
% plot(HM,HoneystoreP,'--bs','LineWidth',2);
% plot(HM,HoneystoreL,'--gs','LineWidth',2);
% legend('Change with ForagerBee','Change with HouseBee','Change with NurseBee','Change with Pupa','Change with Larva')
% % 
% 
% figure(16);
% plot(HM,AdultBees1,'--bo','LineWidth',2);
% hold on 
% plot(HM,AdultBees2,'--rs','LineWidth',2);
% xlabel('HouseBee Survival')
% legend('First November','Second November')
% ylabel('AdultBees')
% 
% figure(3);
% plot(EM,AdultBees1,'--bo','LineWidth',2);
% hold on 
% plot(EM,AdultBees2,'--rs','LineWidth',2);
% xlabel('Egg Survival')
% legend('First November','Second November')
% ylabel('AdultBees Number')
% % % 
% figure(5);
% plot(NTM,Honeystore2,'--rs','LineWidth',2);
% hold on 
% plot(NTM,AdultBees2,'--bo','LineWidth',2);
% legend('Honeystore','Adultbees')
% % % % % % % % 
figure(1);
YMatrix1=pop';
A=Ppop.*0.23/1000;
B=Hpop*0.5/1000;
YMatrix2= [A;B]';
Y1=Rpop;
createfigure1(YMatrix1, YMatrix2, Y1); 
% 
% stem3 (NTM, Honeystore2, AdultBees2, 'DisplayName', 'QP, Honeystore, AdultBees'); figure(gcf);
% % xlabel('precicous behavior probability');
% xlabel('Nectar foraging time');
% ylabel('Honeystore at November');
% zlabel('Adult bees Number at November');
% 
% 
% plot(BRM,Honeystore1,'--bo','LineWidth',2);
% hold on 
% plot(BRM,Honeystore2,'--rs','LineWidth',2);
% % xlabel('Nectar Foraging Time')
% % legend('First November','Second November')
% % ylabel('Honeystore')
% 
% plot(BRM,AdultBees1,'--bo','LineWidth',2);
% hold on 
% plot(BRM,AdultBees2,'--rs','LineWidth',2);
% legend('First November','Second November')
% xlabel('Brood Rearing Efficiency')
% ylabel('AdultBees Number');
% 
% plot(UM,Honeystore1,'--bo','LineWidth',2);
% hold on 
% plot(UM,Honeystore2,'--rs','LineWidth',2);
% xlabel('Nectar Foraging Time')
% legend('First November','Second November')
% ylabel('Honeystore')
% 
% 
% plot(UM,AdultBees1,'--bo','LineWidth',2);
% hold on 
% plot(UM,AdultBees2,'--rs','LineWidth',2);
% xlabel('Nectar Foraging Time')
% legend('First November','Second November')
% ylabel('AdultBees Number')
% 
% figure(13);
% plot(RTM,Honeystore1,'--bo','LineWidth',2);
% hold on 
% plot(RTM,Honeystore2,'--rs','LineWidth',2);
% xlabel('Retarded Behavior Probability')
% legend('First November','Second November')
% ylabel('Honeystore')
% 
% figure(1);
% plot(FRM,AdultBees1,'--bo','LineWidth',2);
% hold on 
% plot(FRM,AdultBees2,'--rs','LineWidth',2);
% xlabel('Foraging Efficiency')
% legend('First November','Second November')
% ylabel('AdultBees')
% 
% figure(2);
% plot(FRM,Honeystore1,'--bo','LineWidth',2);
% hold on 
% plot(FRM,Honeystore2,'--rs','LineWidth',2);
% xlabel('Foraging Efficiency')
% legend('First November','Second November')
% ylabel('Honeystore')


