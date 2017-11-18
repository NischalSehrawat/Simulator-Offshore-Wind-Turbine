function [A_hub,U_nodes,t,a222] = kaimal_sim(I_ref,L_u,U_mean,fs,nfreq,rand_seed,N_t,L_t,N_tp,L_tp)
                %  KAIMAL_SIM 1st generates time series of wind at hub and then gives the
                %  profile on exposed nodes and time vector
                %  It takes the following inputs
                % I_ref = turbulence class coefficient {0.12 0.14 0.16}
                % L_u= length parameter [m]
                % U_mean= mean wind velocity at the hub [m/s]
                % fs = sampling frequency or number of samples taken per second [Hz]
                % nfreq= number of frequencies in the spectrum, decide the
                % maximum simulation time
                % rand_seed = set seed number so that each time same set of
                % random numbers are generated
                % N_t= number of elements in the tower section
                % L_t= length of tower section [m]
                % N_tp= number of elements in the transition piece section 
                % L_tp= length of transition piece [m]

f_max=fs/3; % fs>2f_max Nyquist rule therefore instead of 2, 3 is taken

df=f_max/nfreq; % Frequency increments [Hz]

f_vec=0:df:f_max; % defining frequency vector

t_sim_max=1/(df);% (= nfreq/f_max) maximum simulation time depends on number of frequencies in the spectrum

sigma=I_ref*(0.75*U_mean+5.6);

S_f=((4.*L_u.*sigma.^2./U_mean)./(1+6.*L_u.*f_vec./U_mean).^(5/3)); % calculating frequency spectrum

qq=1;

while f_vec(qq)<0.02
        
    S_f(qq)=0;
    
    qq=qq+1;

end

a222(1,:)=f_vec;
a222(2,:)=S_f;

plot(f_vec,S_f)
hold on

xlabel('\textbf{Frequency}  \textbf{[Hz]}','Interpreter','latex','FontSize',12)

plot([0.27 0.27],[0 15],'b')

amp=sqrt(2.*S_f.*df); % calculating amplitudes of the frequencies

% Generating time series of wind velocities 

rng(rand_seed,'v5uniform'); % set random generator

phase=2*pi*rand(1,nfreq+1);% create random phase angles

t=0:0.05:t_sim_max; % time vector

for i=1:length(t)

    u(i)=sum(amp.*sin(2*pi.*f_vec.*t(i)+phase)); % Turbulent component of wind velocities at hub
    
    u_dot(i)=sum(amp.*2*pi.*f_vec.*cos(2*pi.*f_vec.*t(i)+phase)); % wind accelerations at the hub
end

U_hub=U_mean+u;

A_hub=u_dot;% Time series of wind accelerations at the hub

n_nodes_t=N_t+1; % number of nodes = number of elements +1 for tower section

n_nodes_tp=N_tp+1; % number of nodes = number of elements +1 for transition piece section

L_tot=L_t+L_tp; % total length of the section exposed to wind

l_t=L_t/N_t; % length of one element for tower section

l_tp=L_tp/N_tp; % length of one element for transition piece section

% velocity profile for tower section 

for i=1:1:n_nodes_t
    U_nodes_t(i,:)=U_hub(1,:)*((L_tot-(i-1)*l_t)/L_tot)^0.14;
end

% velocity profile for transition piece section 
for i=1:1:n_nodes_tp
        U_nodes_tp(i,:)=U_hub(1,:)*((L_tot-L_t-(i-1)*l_tp)/L_tot)^0.14;
end
 U_nodes=[U_nodes_t;U_nodes_tp(2:end,:)];

end

