NR = 1; %: no. of search curves
MI = 4; %: maximum number of fourier
k=6;wantedN=700;
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