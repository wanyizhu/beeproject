
k=6;
tX=ones(k+1,k);
Ty2=zeros(1,k+1);
for i=1:k+1 
    surX=tX(i,:);
Ty2(i)=feval(@testNOP, surX)';
end 
% for i=1:k+1 
%     surX=X(i,:);
% %st1=Pamx(i,1);st2=Pamx(i,2);st3=Pamx(i,3);mt4=1-Pamx(i,4);st5=Pamx(i,5);st6=Pamx(i,6);
% %Pamx_i=Pamx(i,:);
% 
% %Pamx=[0.8 0.8 0.8 0.05 0.8 0.8];
%     ttY(i,:) = feval(@testNOP, surX);
% end 
% Ty2=ttY';