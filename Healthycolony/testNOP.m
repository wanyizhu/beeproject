% First Season Summer Dynamics 

agemax = 60; % +1 because of matlab indexing


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
G = zeros(6,agemax);

tx=240; % number of summer days of the year

G(1,1:3)=1; G(2,4:11)=1; G(3,12:26)=1; G(4,27:42)=1;G(5,43:48)=1;G(6,49:agemax)=1;

P0 = 200;
%P0 = 1000; %initial cells of pollen

V0 = 300000 - P0; %intial vacant cells, total number cells is 140000
%subtract to leave room for eggs and pollen

H0=0; %initial  honey

R0=0; %initial ...? egg laying?

N = zeros(agemax,1);

N(1:3)=100; %  initial number of eggs/3 days
%SHOULD BE ZERO, BUT THIS CAUSES CODE TO CRASH WITH ERROR "DEAD HIVE
% Undefined function or variable "storedhoney".
% 
% Error in bees (line 177)
% Ht1= min(V0,Ht-honeyeaten+storedhoney); % The net honey storage
% at the end of the day.
% 
% Error in testNOP (line 110)
% 		     X = bees(X,t);  % call to bees.m function, which
%              outputs new state of hive"

N(4:11)=200; % initial number of larva = 1600/8 days

N(12:26)=160; % initial number of pupa = 2400/15 days

N(27:42)=187; % initial number of nurse bees = 3000/16 days

N(43:48)=500; % initial number of house bees= 3000/ 6 days

N(49:agemax)=250; % initial number of forager bees = 3000 / 12 days

X = [ V0; P0; H0;R0; N ]; % This hold the initial bee populion that goes into bees.m

res=zeros(6,tx); % res will hold bee population by stage for each day of summer

V=zeros(1,tx); %vector will hold # vacant cells for each day of summer

P=zeros(1,tx);

H=zeros(1,tx);

R=zeros(1,tx);

numyears = 2;
summerdays = 240;
yeardays = 360;

%these super long vectors hold the vacant cells, pollen, honey, and egg
%filled cells for every day of the years in our time series
pop=zeros(6,yeardays*numyears);
Vpop=zeros(1,yeardays*numyears);
Ppop=zeros(1,yeardays*numyears);
Hpop=zeros(1,yeardays*numyears);
Rpop=zeros(1,yeardays*numyears);

%each year starts with a field season, goes through one winter, and then
%one more field season
for T = 0:(numyears-1) %T tells us what year we are in 0,1, 2...

          for t=(yeardays*T+1):(yeardays*T+summerdays) %sets the date, goes through all field season days
   
		     X = bees(X,t);  % call to bees.m function, which outputs new state of hive
              
             %G is 6 x agemax, and X = [V,P,H,R,N]
		     res(1:6,t-yeardays*T)=G*X(5:end); 
 
		     V(1,t-yeardays*T)= X(1);

		     P(1,t-yeardays*T) = X(2);
        
             H(1,t-yeardays*T)= X(3);
             % disp([t,X(3)]);

		     R(1,t-yeardays*T)= X(4);
 
          end %END OF LOOP THROUGH THIS SUMMER
     
	pop(:,(yeardays*T+1):(yeardays*T+summerdays))=res;
    Vpop(:,(yeardays*T+1):(yeardays*T+summerdays))=V;
    Ppop(:,(yeardays*T+1):(yeardays*T+summerdays))=P;
    Hpop(:,(yeardays*T+1):(yeardays*T+summerdays))=H;
    Rpop(:,(yeardays*T+1):(yeardays*T+summerdays))=R;
    
	
    
    % First Season Winter Dynamics 

	agemaxwinter=150; %max life of winter bee

	W = zeros(4,agemaxwinter);

	W(1,1:3)=1; W(2,4:11)=1; W(3,12:26)=1; W(4,27:agemaxwinter)=1;

	N = zeros(agemaxwinter,1);

	N(1:3) = res(1,summerdays)/3;

	N(4:11) = res(2,summerdays)/8;

	N(12:26) = res(3,summerdays)/15;

    N(27:agemax) = (res(4,summerdays)+res(5,summerdays)+res(6,summerdays))/34;
    % this doesn't make sense- it's like artificially aging many of them.  
	%N(27:agemaxwinter)=(res(4,summerdays)+res(5,summerdays)+res(6,summerdays))/124; 

	P0 = P(1,summerdays);

    V0 = V(1,summerdays); 

    H0 = H(1,summerdays);

    R0 = R(1,summerdays);
    

    Y = [ V0; P0; H0; R0; N ];

	clear res V P H R;

	res=zeros(6,yeardays-summerdays);
    
    V=zeros(1,yeardays-summerdays);

    P=zeros(1,yeardays-summerdays);

    H=zeros(1,yeardays-summerdays);

    R=zeros(1,yeardays-summerdays);
    	

    for t = (yeardays*T+summerdays+1):(yeardays*(T+1))
        
        Y = winterbeesR(Y,t);
        
        wintpop = W*Y(5:end); %this is a 4xmaxagewinter matrix %this pieces causes second season crash- not just bad pic
        res(1:3,(t-(yeardays*T+summerdays))) = wintpop(1:3);
        res(5,(t-(yeardays*T+summerdays))) = wintpop(4);
        
        V(1,(t-(yeardays*T+summerdays))) = Y(1);
        if Y(1)== 0
            disp('ran out of space, on day:')
            disp(t)
            break
        end
        P(1,(t-(yeardays*T+summerdays))) = Y(2);
        if Y(2) == 0
            disp('Hive starved, no pollen, on day:')
            disp(t)
            break
        end
        H(1,(t-(yeardays*T+summerdays))) = Y(3);
        
        if Y(3) == 0
            disp('Hive staved, no honey, on day:')
            disp(t)
            break
        end
        R(1,(t-(yeardays*T+summerdays))) = Y(4);
    end %END OF LOOP THROUGH WINTER
    
	pop(:, (yeardays*T+summerdays+1):(yeardays*(T+1))) =res;
    
    Vpop(1,(yeardays*T+summerdays+1):(yeardays*(T+1)))= V;
    
    Ppop(1,(yeardays*T+summerdays+1):(yeardays*(T+1)))= P;
    
    Hpop (1,(yeardays*T+summerdays+1):(yeardays*(T+1)))= H;
    
    Rpop (1,(yeardays*T+summerdays+1):(yeardays*(T+1)))= R;
    
        
	%Second Season Summer Dynamics 

	N = zeros(agemax,1);

	N(1:3)=pop(1,yeardays*(T+1))/3;

	N(4:11)=pop(2,yeardays*(T+1))/8;

	N(12:26)=pop(3,yeardays*(T+1))/15;

	N(27:42)= pop(5,yeardays*(T+1))/34;

	N(43:48)= pop(5,yeardays*(T+1))/34 ;

	N(49:agemax)=pop(5,yeardays*(T+1))/34;

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

end %END OF LOOP THROUGH MULTIPLE YEARS

%for each day, this gives the ratio of eggs+larvae/nurse+house bees
% BARatio=(pop(1,1:360*numyears)+pop(2,1:360*numyears))./(pop(4,1:360*numyears)+pop(5,1:360*numyears)); 
%for each day, this gives the ratio of foragers/nurse+house bees
% FARatio=pop(6,1:360*numyears)./(pop(4,1:360*numyears)+pop(5,1:360*numyears));

YMatrix1=pop';
A=Ppop; %pollen storage throughout all seaseons
B=Hpop;  %honey storage throught all seasons
% A=Ppop.*0.23/1000;
% B=Hpop*0.5/1000;
YMatrix2= [A;B]';
 Y3=Rpop;
%Y3=pop(3)*0.1552/1000+pop(4)*0.2189/1000+pop(5)*0.2189/1000+A+B;
createfigure1(YMatrix1, YMatrix2, Y3); 
 % figure;

% plot(Y3);
% foundationweight = 50.2 * 453.6 /1000;
% 
% Y1=(pop(2)+pop(3))*0.1552/1000+pop(4)*0.2189/1000+pop(5)*0.2189/1000+Ppop.*0.23/1000+Hpop*0.5/1000;
% plot(Y1(1:360));
% t=[0:30:360];months=['Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';'Aug';'Sep';'Oct';'Nov';'Dec';];
% set(gca,'xtick',t)
% set(gca,'xticklabel',months)
% xlabel('Date')
% ylabel('Colony Weight')
% 
% BNy=(BARatio+FARatio)';

