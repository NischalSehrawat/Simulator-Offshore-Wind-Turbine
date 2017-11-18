function [M,K,N1,N2,N_t,N_tp_aw,N_tp_uw,N_pile_uw,N_wind] = Assemble_M_K( n1,n2,n3,n4,n5,M1,K1,M2,K2,M3,K3,M4,K4,M5,K5 )
                                % Assemble_M_K assembles the total system
                                % mass and stiffness matrices from top to
                                % bottom
                                % It takes the following inputs
                                % n1=number of element in tower section
                                % n2=number of element in transition piece above water
                                % n3=number of element in transition piece below water
                                % n4=number of element in UNembedded pile
                                % n5=number of element in embedded pile
                                % M1= tower mass matrix
                                % K1=tower stiffness matrix
                                % M2=transition piece above water mas matrix
                                % K2=transition piece above water stiffness matrix
                                % M3= transition piece below water mas matrix
                                % K3=transition piece below water stiffness matrix
                                % M4=unembed pile mas matrix
                                % K4=unembed pile stiffness matrix
                                % M5=embed pile mas matrix
                                % K5=embed pile stiffness matrix
 %% Begin assembling
N=2*(n1+n2+n3+n4+n5+1); % total number of dof in the final model
% Initialising assembly matrices
M11=zeros(N,N);M22=zeros(N,N);M33=zeros(N,N);M44=zeros(N,N);M55=zeros(N,N);
K11=zeros(N,N);K22=zeros(N,N);K33=zeros(N,N);K44=zeros(N,N);K55=zeros(N,N);

% Assembling 1st matrix in a global matrix
n11=size(M1,1);
M11(1:n11,1:n11)=M1;
K11(1:n11,1:n11)=K1;
% Assembling 2nd matrix in a global matrix
n22=size(M2,1);
M22(n11-1:n11+n22-2,n11-1:n11+n22-2)=M2;
K22(n11-1:n11+n22-2,n11-1:n11+n22-2)=K2;
% Assembling 3rd matrix in a global matrix
n33=size(M3,1);
M33(n11+n22-3:n11+n22+n33-4,n11+n22-3:n11+n22+n33-4)=M3;
K33(n11+n22-3:n11+n22+n33-4,n11+n22-3:n11+n22+n33-4)=K3;
% Assembling 4th matrix in a global matrix
n44=size(M4,1);
M44(n11+n22+n33-5:n11+n22+n33+n44-6,n11+n22+n33-5:n11+n22+n33+n44-6)=M4;
K44(n11+n22+n33-5:n11+n22+n33+n44-6,n11+n22+n33-5:n11+n22+n33+n44-6)=K4;
% Assembling 5th matrix in a global matrix
n55=size(M5,1);
M55(n11+n22+n33+n44-7:n11+n22+n33+n44+n55-8,n11+n22+n33+n44-7:n11+n22+n33+n44+n55-8)=M5;
K55(n11+n22+n33+n44-7:n11+n22+n33+n44+n55-8,n11+n22+n33+n44-7:n11+n22+n33+n44+n55-8)=K5;

% Assembling final matrices
M=M11+M22+M33+M44+M55;
K=K11+K22+K33+K44+K55;

N1=[n11+n22-3,n11+n22+n33+n44-6]; % index of the starting and end nodes where hydrodynamic forces will act

N2=[n11+n22+n33+n44-7,n11+n22+n33+n44+n55-8]; % index of the starting and end nodes where spring forces will act

N_t=[1,n11]; % start and end positions of the tower nodes in main system matrix
N_tp_aw=[n11-1,n11+n22-2]; % start and end positions of the above water tp nodes in main system matrix
N_tp_uw=[n11+n22-3,n11+n22+n33-4]; % start and end positions of the under water tp nodes in main system matrix
N_pile_uw=[n11+n22+n33-5,n11+n22+n33+n44-6]; % start and end positions of the under water pile nodes in main system matrix

N_wind=[1,n11+n22-2]; % start and end positions of nodes in air
end

