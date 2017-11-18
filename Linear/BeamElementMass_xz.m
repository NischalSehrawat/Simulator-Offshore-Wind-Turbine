function y = BeamElementMass_xz(p,A,L,I,alfa,nu)
%BeamElementMass This function returns the element
% mass matrix for a Timoshenko beam
% element with density p (rho),
% Area A, and length L, moment of inertia I, shear coefficient alfa, poison
% ratio nu
% The size of the element mass
% matrix is 4 x 4.

r=sqrt(I/A);

fi=24*alfa*(1+nu)*(r/L)^2;

a=13/35+7*fi/10+1/3*fi^2+6/5*(r/L)^2;

b=9/70+3*fi/10+fi^2/6-6/5*(r/L)^2;

c=(11/210+11*fi/120+fi^2/24+(1/10-fi/2)*(r/L)^2)*L;

d=(13/420+3*fi/40+fi^2/24-(1/10-fi/2)*(r/L)^2)*L;

e=(1/105+fi/60+fi^2/120+(2/15+fi/6+fi^2/3)*(r/L)^2)*L^2;

f=(1/140+fi/60+fi^2/120+(1/30+fi/6-fi^2/6)*(r/L)^2)*L^2;

y=(p*A*L/(1+fi)^2)*[a -c b d;...
                   -c e -d -f;...
                    b -d a c;...
                    d -f c e];
end

