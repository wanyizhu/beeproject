function Ty=function2evaluate(X)
global k;
%fhandle = @testNOP;
% Pamx k+1,K matrix 
Ty=zeros(1,k+1);
for j=1:k+1 
    surX=X(j,:);
%st1=Pamx(i,1);st2=Pamx(i,2);st3=Pamx(i,3);mt4=1-Pamx(i,4);st5=Pamx(i,5);st6=Pamx(i,6);
%Pamx_i=Pamx(i,:);

%Pamx=[0.8 0.8 0.8 0.05 0.8 0.8];
Ty(j)=feval(@testNOP, surX)';
    %Y(i,:) = feval(@testNOP, surX);
end 
