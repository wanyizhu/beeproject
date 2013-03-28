N=100;k=2;
X_uvBAF=zeros(N,k);
for i=1:k 
    X_uvBAF(:,i)=0.2*(randperm(N)'-rand(N,1))/N; 
end % X is N*k input factor matrix 
%Ty=zeros(N,1080);
      for j=1:N
        % X= xmat((i-1)*(k+1)+1:i*(k+1),:);
         %Pamx_j=X(i)
        surX=X_uvBAF(j,:);
        TyuvBAF(j,:)=feval(@testNOP, surX);
    %Y(i,:) = feval(@testNOP, surX);
      end
   

maxuvBAF=zeros(1080,1);minuvBAF=zeros(1080,1);
for i=1:1080
maxuvBAF(i)=max(TyuvBAF(:,i));
minuvBAF(i)=min(TyuvBAF(:,i));
end
 YMATRIX1=[maxuvBAF(540:600) minuvBAF(540:600) maxuvBAF(900:960) minuvBAF(900:960)] ; 
     
%figure(1)
%scatter(X(:,1),Ty)
% hold on 
% plot(X(:,2),Ty)
% plot(X(:,3),Ty)
% plot(X(:,4),Ty)
% plot(X(:,5),Ty)
% plot(X(:,6),Ty)
% plot(X(:,7),Ty)
% plot(X(:,8),Ty)
% 
% plot(X(:,9),Ty)
% plot(X(:,10),Ty)
% 
% figure1 = figure('Color',[1 1 1]);
% i=3;
% YMATRIX=[Ty(i,80:240);Ty(i,440:600);Ty(i,800:960)]'; 
% axes1 = axes('Parent',figure1,...
%     'XTickLabel',{'June','July','Aug','Sep','Oct','Nov','Dec'},...
%     'XTick',[0 30 60 90 120 150 180]);
% box('on');
% hold('all');
% 
% % Create multiple lines using matrix input to plot
% plot1 = plot(YMATRIX,'Parent',axes1,'LineWidth',2.0);
% set(plot1(1),'LineStyle',':','Color',[0 0 1],...
%     'DisplayName','First season');
% set(plot1(2),'Color',[0 0 1],'DisplayName','Second season');
% set(plot1(3),'LineStyle',':','Color',[0 1 0],...
%     'DisplayName','Third season');
