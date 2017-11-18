function [M_system,K_system] = system_assemble(M_xy,K_xy,M_xz,K_xz,M_tor,K_tor)
%system_assemble assebles all the mass and martices together to generate a
% mass and a stiffness matrix of the whole system i.i. XZ plane, YZ plane
% and rotations about z axis
p1=size(M_xy,1);
p2=size(M_xz,1);
p3=size(M_tor,1);
p=p1+p2+p3;
M_system=zeros(p,p);
M_system(1:p1,1:p1)=M_xy;
M_system(p1+1:p1+p2,p1+1:p1+p2)=M_xz;
M_system(p1+p2+1:p,p1+p2+1:p)=M_tor;

K_system=zeros(p,p);
K_system(1:p1,1:p1)=K_xy;
K_system(p1+1:p1+p2,p1+1:p1+p2)=K_xz;
K_system(p1+p2+1:p,p1+p2+1:p)=K_tor;


