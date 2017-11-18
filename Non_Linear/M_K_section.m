function [M_xz,K_xz,M_yz,K_yz,M_tor,K_tor] = M_K_section(N,L,p,A,M_axial,m,E,I,alfa,nu,J,G)
                                           % M_K_section  gives mass & stiffness matrices of any section in xy and xy plane including the
                                           % effects of axial compression
                                           % it takes the following inputs
                                           %N= number of elements; L=Length of the section;p=density of steel;
                                           %A= area of section; M_axial=axial mass [Ton] on top of that section;
                                           % m= mass per unit length; E= young's modulus;I=moment of inertia
%% Begin of function
n_nodes=N+1; % number of nodes are equal to number of elements +1
% A beam element has 2 DOF per node i.e. translation (v) and rotation (theta) 
% Therefore total number of degrees of freedom 

n_dof=2*n_nodes;  
l=L/N;% length of one element [m]

%================================================Mass matrix for XZ plane=========================

M_section_xz=zeros(n_dof,n_dof); % initializing section mass matrix for xz plane
m_el_xz=BeamElementMass_xz(p,A,l,I,alfa,nu); % Generating element mass matrix for xz plane

%=================================================Mass matrix for YZ plane=========================

M_section_yz=zeros(n_dof,n_dof); % initializing section mass matrix for yz plane
m_el_yz=BeamElementMass_yz(p,A,l,I,alfa,nu); % Generating element mass matrix for yz plane
i=1;j=2; % initialising counters for matrix assembly

% Mass Matrix Assembly

while((i<n_nodes)&&(j<=n_nodes))
         M_section_xz=BeamAssemble(M_section_xz,m_el_xz,i,j);
         M_section_yz=BeamAssemble(M_section_yz,m_el_yz,i,j);
        i=i+1;j=j+1;
end

% calculating axial compressive force in each element

P_el=zeros(N,1);
for k=1:1:N
P_el(k)=(M_axial +m*l/2+(k-1)*m*l)*9.8*1000; %  axial force in  element number 'k' [N]
end

% Stiffness matrices & assembly

%===========================================Stiffness matrix for XZ plane=========================

K_section_xz=zeros(n_dof,n_dof); % initializing system mass matrix for xz plane

%===========================================Stiffness matrix for YZ plane=========================

K_section_yz=zeros(n_dof,n_dof); % initializing system mass matrix for yz plane

i=1;j=2; % initialising counters for matrix assembly

while((i<n_nodes)&&(j<=n_nodes))
            k_el_xz=BeamElementStiffness_xz(E,I,l,P_el(i),A,alfa,nu);
            K_section_xz=BeamAssemble(K_section_xz,k_el_xz,i,j);
            
            k_el_yz=BeamElementStiffness_yz(E,I,l,P_el(i),A,alfa,nu);
            K_section_yz=BeamAssemble(K_section_yz,k_el_yz,i,j);
        i=i+1;j=j+1;
end
%=============================================== Torsion matrices=================================

% Torsion stiffness matrix
K_torsion=zeros(n_nodes,n_nodes);
k_tor=BeamTorsionStiffness(J,G,l);

% Torsion stiffness Matrix Assembly
i=1;j=2; % initialising counters for matrix assembly

while((i<n_nodes)&&(j<=n_nodes))
     
            K_torsion=BeamTorsionAssemble(K_torsion,k_tor,i,j);
        
            i=i+1;j=j+1;
end

% Torsion Mass Matrix Assembly
% Torsion mass matrix

M_torsion=zeros(n_nodes,n_nodes);

m_tor=BeamTorsionMass(p,l,J);

i=1;j=2; % initialising counters for matrix assembly

while((i<n_nodes)&&(j<=n_nodes))
     
            M_torsion=BeamTorsionAssemble(M_torsion,m_tor,i,j);
        
            i=i+1;j=j+1;
end

%================================================XZ plane mass and stiffness matrices================

M_xz=M_section_xz;
K_xz=K_section_xz;

%========================================YZ plane mass and stiffness matrices=========================

M_yz=M_section_yz;
K_yz=K_section_yz;

%=======================================Torsion mass and stiffness matrices===========================
M_tor=M_torsion;
K_tor=K_torsion;
end
