function [Force_MacCamy] = Mac(zeta_A,omega,wavenum,D_out,depth,t,phase,N_tp_uw,L_tp_uw,N_pile_uw,L_pile_uw)

k=wavenum;

n_nodes_tp_uw=N_tp_uw+1; % number of nodes = number of elements +1 for underwater transition piece

n_nodes_pile_uw=N_pile_uw+1; % number of nodes = number of elements +1 for un-embedded pile section

N_tot=N_tp_uw+N_pile_uw+1; % total number of nodes in the pile+transition piece system

l_tp_uw=L_tp_uw/N_tp_uw; % length of one element for underwater transition piece

l_pile_uw=L_pile_uw/N_pile_uw; % length of one element for un-embedded pile section

%% Bessel function evaluations
x=k.*D_out*0.5;
dx=0.000001;
J1dx=besselj(1,x+dx);
Y1dx=bessely(1,x+dx);
J1mdx=besselj(1,x-dx);
Y1mdx=bessely(1,x-dx);
J1der=(J1dx-J1mdx)/(2*dx);
Y1der=(Y1dx-Y1mdx)/(2*dx);
A=1./sqrt(J1der.^2+Y1der.^2);
alfa=-atan(Y1der./J1der);

%%  Underwater transition piece

% force on 1st node of underwater transition piece 

for i=1:length(t)

    f_1_tp(i)=0.5*l_tp_uw*4*1025*9.81.*sum((zeta_A.*A.*cos(omega.*t(i)-alfa+phase))./(k));
end

% force on intermediate nodes of underwater transition piece 

for i=2:1:(n_nodes_tp_uw-1)
 for j=1:length(t)
  
     f_int_tp(i-1,j)=l_tp_uw*4*1025*9.81.*sum((zeta_A.*A.*(cosh(k.*(depth-(i-1).*l_tp_uw))./cosh(k.*depth)).*cos(omega.*t(j)-alfa+phase))./(k));

 end
end

% force on last node of underwater transition piece

for i=1:length(t)

    f_last_tp(i)=0.5*l_tp_uw*4*1025*9.81.*sum((zeta_A.*A.*(cosh(k.*(depth-L_tp_uw))./cosh(k.*depth)).*cos(omega.*t(i)-alfa+phase))./(k));
end

F_tp=[f_1_tp;f_int_tp;f_last_tp];

%%  Unembedded pile

% force on 1st node of unembedded pile 

for i=1:length(t)

    f_1_pile(i)=0.5*l_pile_uw*4*1025*9.81.*sum((zeta_A.*A.*(cosh(k.*(depth-L_tp_uw))./cosh(k.*depth)).*cos(omega.*t(i)-alfa+phase))./(k));
end

% force on intermediate nodes of unembedded pile 

for i=2:1:(n_nodes_pile_uw-1)
 for j=1:length(t)
  
     f_int_pile(i-1,j)=l_pile_uw*4*1025*9.81.*sum((zeta_A.*A.*(cosh(k.*(depth-L_tp_uw-(i-1).*l_tp_uw))./cosh(k.*depth)).*cos(omega.*t(j)-alfa+phase))./(k));

 end
end

% force on last node of unembedded pile 

for i=1:length(t)

    f_last_pile(i)=0.5*l_pile_uw*4*1025*9.81.*sum((zeta_A.*A.*(cosh(k.*(depth-L_tp_uw-L_pile_uw))./cosh(k.*depth)).*cos(omega.*t(i)-alfa+phase))./(k));
end

F_pile=[f_1_pile;f_int_pile;f_last_pile];

% Assembling all the forces 

f_tp=zeros(N_tot,length(t));

f_tp(1:n_nodes_tp_uw,:)=F_tp;

f_pile=zeros(N_tot,length(t));

f_pile(n_nodes_tp_uw:end,:)=F_pile;


Force_MacCamy=f_tp+f_pile;


end
