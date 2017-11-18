function [vel,acc] = vel_acc(zeta_A,omega,phase,wavenum,t,depth,N_tp_uw,L_tp_uw,N_pile_uw,L_pile_uw)

% assumed wave profile is n=acos(wt-kx)

k=wavenum;
x=0;

n_nodes_tp_uw=N_tp_uw+1; % number of nodes = number of elements +1 for underwater transition piece

n_nodes_pile_uw=N_pile_uw+1; % number of nodes = number of elements +1 for un-embedded pile section

l_tp_uw=L_tp_uw/N_tp_uw; % length of one element for underwater transition piece

l_pile_uw=L_pile_uw/N_pile_uw; % length of one element for un-embedded pile section

% Profile for underwater transition piece

for i=1:1:n_nodes_tp_uw
  for j=1:length(t)
  
     vel_tp_uw(i,j)=sum(omega.*zeta_A.*cos(omega.*t(j)+phase-k.*x).*cosh(k.*(depth-(i-1).*l_tp_uw))./sinh(k.*depth));
     
     acc_tp_uw(i,j)=-sum(omega.^2.*zeta_A.*sin(omega.*t(j)+phase-k.*x).*cosh(k.*(depth-(i-1).*l_tp_uw))./sinh(k.*depth));
  
 end
end

% Profile for un-embedded pile section 

for i=1:1:n_nodes_pile_uw
  for j=1:length(t)
  
     vel_pile_uw(i,j)=sum(omega.*zeta_A.*cos(omega.*t(j)+phase-k.*x).*cosh(k.*(depth-L_tp_uw-(i-1).*l_pile_uw))./sinh(k.*depth));
     
     acc_pile_uw(i,j)=-sum(omega.^2.*zeta_A.*sin(omega.*t(j)+phase-k.*x).*cosh(k.*(depth-L_tp_uw-(i-1).*l_pile_uw))./sinh(k.*depth));
  
 end
end


vel=[vel_tp_uw;vel_pile_uw(2:end,:)];

acc=[acc_tp_uw;acc_pile_uw(2:end,:)];

end

