function [M,K] = Assemble_tor(n1,n2,n3,n4,n5,M1,K1,M2,K2,M3,K3,M4,K4,M5,K5)
 % Assemble_tor assembles the torsional mass and stiffness matrices
                                % from top to bottom
                                % It takes the following inputs
                                % n1=number of element in tower section
                                % n2=number of element in transition piece above water
                                % n3=number of element in transition piece below water
                                % n4=number of element in UNembedded pile
                                % n5=number of element in embedded pile
                                % M1= tower torion mass matrix
                                % K1=tower torsion stiffness matrix
                                % M2=transition piece above water torsion mass matrix
                                % K2=transition piece above water torsion stiffness matrix
                                % M3= transition piece below water torsion mass matrix
                                % K3=transition piece below water torsion stiffness matrix
                                % M4=unembed pile torsion mass matrix
                                % K4=unembed pile torsion stiffness matrix
                                % M5=embed pile torsion mass matrix
                                % K5=embed pile torsion stiffness matrix
                                
%% Begin assembling
N=n1+n2+n3+n4+n5+1; % total number of dof in the final model for torsional matrices
% Initialising assembly matrices
M11=zeros(N,N);M22=zeros(N,N);M33=zeros(N,N);M44=zeros(N,N);M55=zeros(N,N);
K11=zeros(N,N);K22=zeros(N,N);K33=zeros(N,N);K44=zeros(N,N);K55=zeros(N,N);

% Assembling 1st torsional matrix in a global matrix
n11=size(M1,1);
M11(1:n11,1:n11)=M1;
K11(1:n11,1:n11)=K1;

% Assembling 2nd matrix in a global matrix
n22=size(M2,1);
M22(n11:n11+n22-1,n11:n11+n22-1)=M2;
K22(n11:n11+n22-1,n11:n11+n22-1)=K2;
% Assembling 3rd matrix in a global matrix
n33=size(M3,1);
M33(n11+n22-1:n11+n22+n33-2,n11+n22-1:n11+n22+n33-2)=M3;
K33(n11+n22-1:n11+n22+n33-2,n11+n22-1:n11+n22+n33-2)=K3;
% Assembling 4th matrix in a global matrix
n44=size(M4,1);
M44(n11+n22+n33-2:n11+n22+n33+n44-3,n11+n22+n33-2:n11+n22+n33+n44-3)=M4;
K44(n11+n22+n33-2:n11+n22+n33+n44-3,n11+n22+n33-2:n11+n22+n33+n44-3)=K4;

% Assembling 5th matrix in a global matrix
n55=size(M5,1);
M55(n11+n22+n33+n44-3:n11+n22+n33+n44+n55-4,n11+n22+n33+n44-3:n11+n22+n33+n44+n55-4)=M5;
K55(n11+n22+n33+n44-3:n11+n22+n33+n44+n55-4,n11+n22+n33+n44-3:n11+n22+n33+n44+n55-4)=K5;
% Assembling final matrices
M=M11+M22+M33+M44+M55;
K=K11+K22+K33+K44+K55;
end

