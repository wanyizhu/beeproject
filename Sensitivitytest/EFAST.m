function [Si,Sti] = EFAST(k,wantedN)
%[SI,STI] = EFAST(K,WANTEDN)
% First order and total effect indices for a given
% model computed with ‘Extended Fourier Amplitude
% Sensitivity Test’ (EFAST).
% Andrea Saltelli, Stefano Tarantola and Karen Chan.
% 1999
% ‘A quantitative model-independent method for global
% sensitivity analysis of model output.’
% Technometrics 41:39-56
% Input:
% K : no. of input factors
% WANTEDN : wanted no. of sample points
% Output:
% SI[] : first order sensitivity indices
% STI[] : total effect sensitivity indices
% Other used variables/constants:
% OM[] : vector of ‘k’ frequencies
% OMI : frequency for the group of interest
% OMCI[] : set of freq. used for the compl. group
% X[] : parameter combination rank matrix
% AC[],BC[]: fourier coefficients
% FI[] : random phase shift
% V : total output variance (for each curve)
% VI : partial var. of par. i (for each curve)
% VCI : part. var. of the compl. set of par...
% AV : total variance in the time domain
% AVI : partial variance of par. i
% AVCI : part. var. of the compl. set of par.
% Y[] : model output
% N : no. of runs on each curve
NR = 1; %: no. of search curves
MI = 4; %: maximum number of fourier
global N;
% coefficients that may be
% retained in calculating the
% partial variances without
% interferences between the
% assigned frequencies
%
% Computation of the frequency for the group
% of interest ‘OMi’ and the no. of sample points ‘N’.
OMi = floor((wantedN/NR-1)/(2*MI)/k);
N = 2*MI*OMi+1;
if(N*NR < 65)
fprintf(['Error: sample size must be >= ' ...
'65 per factor.\n']);
return;
end
% Algorithm for selecting the set of frequencies.
% OMci(i), i=1:k-1, contains the set of frequencies
% to be used by the complementary group.
OMci = SETFREQ(N-1,OMi/2/MI);
% Loop over the ‘k’ input factors.
for i=1:k
% Initialize AV,AVi,AVci to zero.
AV = 0;
AVi = 0;
AVci = 0;
% Loop over the ‘NR’ search curves.
for L=1:NR
% Setting the vector of frequencies ‘OM’
% for the ‘k’ factors.
cj = 1;
for j=1:k
if(j==i)
% For the factor of interest.
OM(i) = OMi;
else
% For the complementary group.
OM(j) = OMci(cj);
cj = cj+1;
end
end
% Setting the relation between the scalar
% variable ‘S’ and the coordinates
% {X(1),X(2),...X(k)} of each sample point.
FI = rand(1,k)*2*pi; % random phase shift
S_VEC = pi*(2*(1:N)-N-1)/N;
OM_VEC = OM(1:k);
FI_MAT = FI(ones(N,1),1:k)';
ANGLE = OM_VEC'*S_VEC+FI_MAT;
X = 0.5+asin(sin(ANGLE'))/pi;
% Transform distributions from standard
% uniform to general.
%X = distTransform(X);
a=0.98;b=1;
r=size(X,1);
%Y = a*ones(r,k) + (b-a)*X;
X = (b-a)*X + a*ones(r,k);
% Do the ‘N’ model evaluations.
      for j=1:N
         %Pamx_j=X(i)
        surX=X(j,:);
        Y(j)=feval(@testNOP, surX)';
    %Y(i,:) = feval(@testNOP, surX);
      end
%Y = function2evaluate(X);
% Subtract the average value.
Y = Y-mean(Y);
% Fourier coeff. at [1:OMi/2].
NQ = (N-1)/2;
N0 = NQ+1;
COMPL = 0;
Y_VECP = Y(N0+(1:NQ))+Y(N0-(1:NQ));
Y_VECM = Y(N0+(1:NQ))-Y(N0-(1:NQ));
for j=1:OMi/2
ANGLE = j*2*(1:NQ)*pi/N;
C_VEC = cos(ANGLE);
S_VEC = sin(ANGLE);
AC(j) = (Y(N0)+Y_VECP*C_VEC')/N;
BC(j) = Y_VECM*S_VEC'/N;
COMPL = COMPL+AC(j)^2+BC(j)^2;
end
% Computation of V_{(ci)}.

Vci = 2*COMPL;
AVci = AVci+Vci;
% Fourier coeff. at [P*OMi, for P=1:MI].
COMPL = 0;
Y_VECP = Y(N0+(1:NQ))+Y(N0-(1:NQ));
Y_VECM = Y(N0+(1:NQ))-Y(N0-(1:NQ));
for j=OMi:OMi:OMi*MI
ANGLE = j*2*(1:NQ)*pi/N;
C_VEC = cos(ANGLE');
S_VEC = sin(ANGLE');
AC(j) = (Y(N0)+Y_VECP*C_VEC)/N;
BC(j) = Y_VECM*S_VEC/N;
COMPL = COMPL+AC(j)^2+BC(j)^2;
end
% Computation of V_i.
Vi = 2*COMPL;
AVi = AVi+Vi;
% Computation of the total variance
% in the time domain.
AV = AV+Y*Y'/N;
end % L=1:NR
% Computation of sensitivity indicies.
AV = AV/NR;
AVi = AVi/NR;
AVci = AVci/NR;
Si(i) = AVi/AV;
Sti(i) = 1-AVci/AV;
end % i=1:k
% Algorithm for selection of a frequency
% set for the complementary group. Done
% recursively as described in:
% Appendix of ‘Sensitivity Analysis’
% [Saltelli et al.]


