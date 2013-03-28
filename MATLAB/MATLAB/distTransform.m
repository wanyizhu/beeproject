function Y = distTransform(X)
global k ; % the bounder of x value 
%Inverse Sampling: Transform distributions from uniform to general
%   x is a matrix storing CDF p.
% p(p<0)=NaN;
% p(p>1)=NaN;
a=0.5;b=1;
r=size(X,1);
%Y = a*ones(r,k) + (b-a)*X;
Y = (b-a)*X + a*ones(r,k);
%Y = (b-a)*X;
return

%X=65*10; k=10;