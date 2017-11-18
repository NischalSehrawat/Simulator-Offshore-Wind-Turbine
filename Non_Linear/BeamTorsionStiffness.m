function K = BeamTorsionStiffness(J,G,L)

% BeamTorsionStiffness calculates the torsion stiffness matrix associated with the beam element 

K=(J*G/L)*[1 -1;-1 1];

end

