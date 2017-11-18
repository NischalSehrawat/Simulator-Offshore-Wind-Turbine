function [M] = BeamTorsionMass(p,L,J)
% BEAMTORSIONMASS calculates the mass matrices associates with the rotating
% inertia of a beam element
M=p*L*J*[1/3 1/6;1/6 1/3];
end

