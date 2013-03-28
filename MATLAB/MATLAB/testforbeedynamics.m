%nrm = [1/3,1/8,1/15,1/16,1/6,1/12]'; % normalize frequencies for each
%stage

% First Season Summer Dynamics 


global st1 st2 st3 mt4 st5 st6; % st1,2,3,5,6: survival rate for each stage of bees; mt4-the mortality rate of nurse bee stage; 
global  BR; % The brood rearing efficiency: the constant ratio of egg to nurse bees as signal for queen to decide the egg laying rate. 
global u v rt; % The probability of precocious foraging (u), reversed behavior of foragers to in-hive behaviors(v), the delayed development of adult bees at each age(rt); 
global a1 a2 a4 a5 h2 h4 h5 h6; %%%cosumption rate of honey and pollen  for each stage of bees 
global hsurfX hsurfY hsurf; % the interpolated surface of NectarODE and honeycollection 

 
 %% Parameter Value 
%  Fraction of a cellful honey consumed by a bee in one day, a cellful of honey weighs~~0.5g
h1=0;
h2=0.012; 
h3=0;
h4=0.042;
h5=0.036;
h6=0.09; 

% Fraction of a cellful pollen consumed by a bee in one day,a cellful of
% pollen weighs~~0.23g
a1 = .02;  
a2 = .02; % fraction of a cell's pollen consumed by a larva in one day
a3 = .0;
a4 = .063; % fraction of a cell's pollen consumed by a nurse in one day
a5 = .003; % fraction of a cell's pollen removed by a house bee in one day
a6=.00;


u=0.0000000;%precocious prob
v=0.0000000;% reversed prob. between foragers and house bees;
BR=4;
rt=0.0000000;

%% The time-independent base survival of bees at different stage
st1= 0.97; % 0.97--survivorship for egg stage-normal condition 
st2= 0.99; % 0.99 survivorship for larval stage-normal condition 
st3= 0.999;% 0.999 survivorship for pupal stage-normal condition 
mt4= 0.005; % 0.005 mortality for nurse stage-normal condition 
st5=0.995; % 0.995  --survivorship for house bee stage-normal condition 
st6=0.965; % 0.965--survivorship for forager bee stage-normal condition 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
agemax=60; % indexing in matlab starts at 1, so add an extra day
G = zeros(6,agemax);
tx=240;
G(1,1:3)=1; G(2,4:11)=1; G(3,12:26)=1; G(4,27:42)=1;G(5,43:48)=1;G(6,49:agemax)=1;

P0 = 1000;

V0 = 100000 - P0;

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


numyears =5;
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
%NOP-normal operating point 


BNRatio=(pop(1,1:360*numyears)+pop(2,1:360*numyears))./(pop(4,1:360*numyears)+pop(5,1:360*numyears)); 
FARatio=pop(6,1:360*numyears)./(pop(4,1:360*numyears)+pop(5,1:360*numyears));


figure(1);
YMatrix1=[pop(1,:)+pop(2,:)+pop(3,:);pop(4,:)+pop(5,:)+pop(6,:)]';
A=Ppop;
B=Hpop;
YMatrix2= [A;B]';
Y1=Rpop';
createfigure1(YMatrix1, YMatrix2,Y1); 
