N=20;k=10;
X=zeros(N,k);
for i=1:k 
    X(:,i)=(randperm(N)'-rand(N,1))/N; 
end % X is N*k input factor matrix 
Ty=zeros(N,1);
      for j=1:N
        % X= xmat((i-1)*(k+1)+1:i*(k+1),:);
         %Pamx_j=X(i)
        surX=X(j,:);
        Ty(j)=feval(@testNOP, surX);
    %Y(i,:) = feval(@testNOP, surX);
      end
