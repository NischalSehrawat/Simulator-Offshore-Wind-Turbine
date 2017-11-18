function [M_system,K_system] = system_assemble(M_xz,K_xz,M_yz,K_yz,M_tor,K_tor)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
p1=size(M_xz,1);
p2=size(M_yz,1);
p3=size(M_tor,1);
p=p1+p2+p3;
M_system=zeros(p,p);
M_system(1:p1,1:p1)=M_xz;
M_system(p1+1:p1+p2,p1+1:p1+p2)=M_yz;
M_system(p1+p2+1:p,p1+p2+1:p)=M_tor;

K_system=zeros(p,p);
K_system(1:p1,1:p1)=K_xz;
K_system(p1+1:p1+p2,p1+1:p1+p2)=K_yz;
K_system(p1+p2+1:p,p1+p2+1:p)=K_tor;


