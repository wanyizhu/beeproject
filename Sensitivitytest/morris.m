function [mu sigma] = morris(k,p,r)
%[MU SIGMA] = MORRIS(K,P,R)
% Comutes the mean and standard deviation of
% the elementary effects.
% Method of Morris.
% Max D. Morris. 1991.
% ‘Factorial Sampling Plans for Preliminary
% Computational Experiments.’
% Technometrics 33:161-174
% Input:
% K : no. of factors
% P : no. of levels
% R : no. of samples
% Output:
% MU : mean of elementary effects
% SIGMA: standard deviation of e.f.
%
% No. of model runs is R*(K+1)
%
% Get the delta value
delta = p/2/(p-1);
% Initial sample matrix B, (k+1,k),
% strictly lower triangular matrix of 1’s.
B = tril(ones(k+1,k),-1);
% Matrix J, (k+1,k) matrix of 1’s.
Jk = ones(k+1,k);
% Get the set of x values, [0,1-delta].
xset = 0:1/(p-1):1-delta;
% Initialize the elementary effects to zero.
eei = zeros(r,k);
indmat = []; xmat = [];
% Generate input data matrix.
for i=1:r
% Get the matrix D*, k-dimensional
% diagonal matrix in which each
% nonzero element is either +1 or
% -1 with equal probability.
Dstar = diag(sign(rand(k,1)*2-1));
% Get matrix P*, (k,k) random
% permutation matrix in which each
% column contains one element equal
% to 1 and all others equal to 0 and
% no two columns have 1s in the same
% position.
I = eye(k);
Pstar = I(:,randperm(k));
% Get vector x*, k-dimensional ’base value’
% of X. Elements are randomly drawn from
% the allowed set.
xstar = xset(ceil(rand(k,1)*floor(p/2)));
% Get matrix B*, (k+1,k) orientation matrix.
% One elementary effect per input.
Bstar = (Jk(:,1)*xstar+(delta/2)*((2*B-Jk)*Dstar+Jk))*Pstar;
% i’th component has been decreased by delta.
[idec,jdec] = find(Dstar*Pstar<0);
% i’th component has been increased by delta.
[iinc,jinc] = find(Dstar*Pstar>0);
% Build up an index matrix.
indmat = [indmat;jdec idec idec+1;jinc iinc+1 iinc];
% Store the x-values.
xmat = [xmat;Bstar];
end
% Transform nonuniform indata from u(0,1)
% to p quantiles.
%ind = find(<nonuniform input factors>);
ind = find(xmat<0 & xmat>1);
xmat(:,ind) = 1/2/p+((p-1)/2/p-1/2/p)* xmat(:,ind);
% % Transform distributions from standard
% % uniform to general.
%xmat = distTransform(xmat);
a=0.98;b=1;
xmat = a*ones(r*(k+1),k) + (b-a)*xmat;
% Perform the r*(k+1) evaluations.
Ty=zeros(r,k+1);
for i=1:r
      for j=1:k+1 
         X= xmat((i-1)*(k+1)+1:i*(k+1),:);
         %Pamx_j=X(i)
        surX=X(j,:);
        Ty(i,j)=feval(@testNOP, surX)';
    %Y(i,:) = feval(@testNOP, surX);
      end
end
        % Y(i,:) = function2evaluate(X);

% Compute the r elementary effects.
for i=1:r
ind = indmat((i-1)*k+1:i*k,:);% 6*3 matrix, 6 is K
eei(i,ind(:,1)) = (Ty(i,ind(:,2))- ...
Ty(i,ind(:,3)))/delta;
end
% Take the absolute value of the r
% elementary effects (Simlab does this..).
eei = abs(eei);
% Get the mean and standard deviation
% of the elementary effects.
mu = mean(eei);
sigma = std(eei);
return 