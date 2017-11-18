function [M_ad] = added_mass(Ca,D,M,Num_el,L)

% added_mass adds the added mass of water to the appropriate place in the system mass mastrix

n_nodes=Num_el+1;
l=L/Num_el;
% Added mas on 1st node
m_1=1025*Ca*pi*D^2/4*l/2;
% Added mass on last node 
m_last=1025*Ca*pi*D^2/4*l/2;

% Added mass on intermediate nodes
m_int(1:n_nodes-2,1)=1025*Ca*pi*D^2/4*l;

m_added=[m_1;m_int;m_last];

j=1;
for i=1:1:n_nodes
    M(j,j)=M(j,j)+m_added(i,1);
    j=j+2;
end
M_ad=M;
end
