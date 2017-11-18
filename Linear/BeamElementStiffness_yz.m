function y = BeamElementStiffness_yz(E,I,L,P,A,alfa,nu)
%BeamElementStiffness  function returns the element 
% stiffness matrix for a Timoshenko beam 
% element with modulus of elasticity E, 
% moment of inertia I, length L, axial force P [N], cross sectional area A, shear coefficient alfa, poison
% ratio nu
% The size of the 
% element stiffness matrix is 4 x 4.
r=sqrt(I/A);
fi=24*alfa*(1+nu)*(r/L)^2;
% fi=0;
y1 = E*I/(L*L*L*(1+fi))*[(12) (6*L) (-12) (6*L);...
                         (6*L) ((4+fi)*L*L) (-6*L) ((2-fi)*L*L) ;...
                         (-12) (-6*L) (12) (-6*L);...
                         (6*L) ((2-fi)*L*L) (-6*L) ((4+fi)*L*L)];
% Effect of axial force P on beam is to reduce the stiffness
y2=-P/(L*(1+fi)^2)*[(6/5+2*fi+fi^2) (L/10) (-6/5-2*fi-fi^2) (L/10);...
                    (L/10) (2*L^2/15+L^2*fi/6+L^2*fi^2/12) (-L/10) (-L^2/30-L^2*fi/6-L^2*fi^2/12);...
                    (-6/5-2*fi-fi^2) (-L/10) (6/5+2*fi+fi^2) (-L/10);...
                      (L/10) (-L^2/30-L^2*fi/6-L^2*fi^2/12) (-L/10) (2*L^2/15+L^2*fi/6+L^2*fi^2/12)];
y=y1+y2;
end
