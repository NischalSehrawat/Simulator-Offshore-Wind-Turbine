tic

clc;clearvars;close all;

HS=[2;4;6;8];UU_wind=[12;15;20;24];beta_waves=[0;30;60;90];theta_pitch=[3.83;10.45;17.43;22.35];ind_fac=[0.14;0.07;0.04;0.03];

for ii=3:1:3
%% Assemble the OWT model from top to bottom using Timoshenko beam elements

% The OWT has 3 components i.e. tower , transition piece and the pile.

% However the entire model will be divided into 5 parts accounting for

% submerged & un-submerged transition piece,embedded and un-embedded pile

%% Discretization definition

N_TWR=15; % number of elements for discretizing tower

N_TP_AW=3; % number of elements for discretizing above water transition piece

N_TP_UW=3; % number of elements for discretizing underwater transition piece

N_pile_UNembed=3; % number of elements for discretizing UN-EMBEDED pile

N_pile_embed=10; % number of elements for discretizing EMBEDED pile

%% hydrodynamic parameters

Ca=1;Cm=1+Ca;Cdw=0.6;

spectrum_type=3; % PM spectrum;

hs=HS(ii); % Significant wave height [m]

T0=5*sqrt(hs) ; % peak period [s]

D_out_pile=5.75; % pile outer dia [m]

depth=25; %  water depth at the site [m]

nfreq=60; % number of frequencies in the spectrum decide the maximum simulation time

freq_cutoff=5; % max frequency [rad/s]

rand_seed=2;

theta_waves=beta_waves(ii); % incoming wave direction angle in degrees

% Getting wave surface elevation parameters from spectrum 

[zeta_A,omega,phase,wavenum,t_sim_max,a111]=create_waves(spectrum_type,hs,T0,depth,nfreq,freq_cutoff,rand_seed);

% plot wave spectrum

% figure
% 
% plot(omega./(2*pi),zeta_A)

t1=0:0.05:t_sim_max-2;

fprintf('\n Maximum time that you can simulate in seconds %0.2f \n\n',t_sim_max);

fprintf(' Time being simulated in seconds %0.2f \n',t_sim_max-2);


Surf_el = surf_el(zeta_A,omega,phase,wavenum,t1)
%% Acc,vel & MacCamy force

L_TP_UW=13.7;% Length of transition piece underwater [m]

L_pile_UNembed=11.3; % UNembedded Pile length [m]

[vel,acc] = vel_acc(zeta_A,omega,phase,wavenum,t1,depth,N_TP_UW,L_TP_UW,N_pile_UNembed,L_pile_UNembed);

% plot(t1,vel(1,:),t1,acc(1,:));

[Force_MacCamy] = Mac(zeta_A,omega,wavenum,D_out_pile,depth,t1,phase,N_TP_UW,L_TP_UW,N_pile_UNembed,L_pile_UNembed);

% figure
%  
% plot(t1,sum(Force_MacCamy));

%% Aerodynamic loading and parameters

I_ref=0.14; % turbulence class

U_mean=UU_wind(ii); % mean wind velocity at the hub [m/s]

fs=5; % sampling frequency [Hz]

a_ind=ind_fac(ii);

red_fac=1-a_ind;

L_blade = 61.5; % blade length [m]

L_u=340.2; % length parameter [m]

L_TWR_total=68; % total length of the Tower [m]

L_TP_AW=5;% Length of transition piece above water [m]

Cd_a=0.6; % drag coefficient for aerodynamic loading

theta_wind=90; % Angle of approach of the wind [degrees]

rpm=12.1; % Revolving speed of the turbine [revolutions per minute '+' for clockwise and - for anticlockwise rotation]

pitch=theta_pitch(ii); % Pitch angle of the blades [degrees]

overhang=5; % overhang of the rotor from the tower [m]

load Blade_properties.mat; % (17x4) vector containing blade element length [m], twist angle [degrees], chord length [m] & element radius [m] in its 1st, 2nd, 3rd and 4th column respectively

[A_hub,U_wind,t_wind,a222] = kaimal_sim(I_ref,L_u,U_mean,fs,7*nfreq,rand_seed,N_TWR,L_TWR_total,N_TP_AW,L_TP_AW); % generating wind velocity time series at nodes subjected to wind loading

%% FEM model

% Parameter definition for Tower

p=7850;% density [kg/m^3]

nu=0.33; % poisson ratio

E=210e9; % Young's modulus  [N/m^2]

G=80e9;  % Shear modulus  [N/m^2]

Mass_TWR_total=348; % Total mass of Tower in [Ton]

% L_TWR_total=68; % total length of the Tower [m] already defined in
% aerodynamic loading

m_TWR=Mass_TWR_total/L_TWR_total; % mass per unit length of the tower in [Ton/m]

M_axial_TWR=(240+110); % Mass of rotor+nacelle [Ton]

D_out_TWR=4.8; % outer dia in [m]

t_TWR=0.043629; % thickness in [m]

D_in_TWR=D_out_TWR-2*t_TWR; % inner dia [m]

alfa_tower=shear_coefficient(D_out_TWR,D_in_TWR,nu); % shear coefficient for hollow cylinders

A_TWR=area(D_out_TWR,D_in_TWR); % Cross sectional area  of tower [m^2]

[I_TWR,J_TWR]=moment_of_inertia(D_out_TWR,D_in_TWR); % area moment of inertia of tower [m^4]

[M_TWR_xz,K_TWR_xz,M_TWR_yz,K_TWR_yz,M_TWR_tor,K_TWR_tor]=M_K_section(N_TWR,L_TWR_total,...
    p,A_TWR,M_axial_TWR,m_TWR,E,I_TWR,alfa_tower,nu,J_TWR,G);

%% Parameter definition for Above Water Transition Piece

D_out_TP=5.75; % outer dia in [m]

t_TP=0.056; % thickness in [m]

D_in_TP=D_out_TP-2*t_TP; % inner dia [m]

alfa_tp=shear_coefficient(D_out_TP,D_in_TP,nu);

Mass_TP_total=147; % Total mass of  transition piece [Ton]

L_TP_total=18.7; % total length of the transition piece [m]

m_TP=Mass_TP_total/L_TP_total; % mass per unit length of the transition piece [Ton/m]

M_axial_TP_AW=(240+110+348);% Mass of rotor+nacelle+tower [Ton]

A_TP=area(D_out_TP,D_in_TP); % Cross sectional area  of transition piece [m^2]

[I_TP,J_TP]=moment_of_inertia(D_out_TP,D_in_TP); % area moment of inertia of tranition piece [m^4]

% L_TP_AW=5;% Length of transition piece above water [m] already defined in
% aerodynamic loading

[M_TP_AW_xz,K_TP_AW_xz,M_TP_AW_yz,K_TP_AW_yz,M_TP_AW_tor,K_TP_AW_tor]=M_K_section(N_TP_AW,L_TP_AW,p,...
    A_TP,M_axial_TP_AW,m_TP,E,I_TP,alfa_tp,nu,J_TP,G);

%% Parameter definition for Underwater Transition Piece

M_axial_TP_UW=(240+110+348+m_TP*5);% Mass of rotor+nacelle+tower+transition piece above water [Ton]

% L_TP_UW=13.7;% Length of transition piece underwater [m] Already define
% above

[M_TP_UW_xz1,K_TP_UW_xz,M_TP_UW_yz1,K_TP_UW_yz,M_TP_UW_tor,K_TP_UW_tor]=M_K_section(N_TP_UW,L_TP_UW,p,...
    A_TP,M_axial_TP_UW,m_TP,E,I_TP,alfa_tp,nu,J_TP,G);

%% Compensating for added mass in underwater transition piece xy plane matrices

M_TP_UW_xz=added_mass(Ca,D_out_TP,M_TP_UW_xz1,N_TP_UW,L_TP_UW);

%% Compensating for added mass in underwater transition piece in xz plane matrices

M_TP_UW_yz=added_mass(Ca,D_out_TP,M_TP_UW_yz1,N_TP_UW,L_TP_UW);

%% Parameter definition for UN-embedded Pile

D_out_pile=5.75; % outer dia in [m]

t_pile=0.1104; % thickness in [m]

D_in_pile=D_out_pile-2*t_pile; % inner dia [m]

alfa_pile= shear_coefficient(D_out_pile,D_in_pile,nu);

A_pile=area(D_out_pile,D_in_pile); % Cross sectional area of the pile [m^2]

[I_pile,J_pile]=moment_of_inertia(D_out_pile,D_in_pile); % area moment of inertia of transition piece [m^4]

Mass_pile_total=542;% Total mass of  pile in [Ton]

L_pile_total=(24+11.3);% total length of the pile [m]

m_pile=Mass_pile_total/L_pile_total; % mass per unit length of the pile in [Ton/m]

M_axial_pile_UNembed=(240+110+348+147); % Mass of rotor+nacelle+tower+transition piece [Ton]

% L_pile_UNembed=11.3; % UNembedded Pile length [m] Already defined above

[M_PILE_UNEMBED_xz1,K_PILE_UNEMBED_xz,M_PILE_UNEMBED_yz1,K_PILE_UNEMBED_yz,M_PILE_UNEMBED_tor,K_PILE_UNEMBED_tor]=M_K_section(N_pile_UNembed,L_pile_UNembed,p,...
    A_pile,M_axial_pile_UNembed,m_pile,E,I_pile,alfa_pile,nu,J_pile,G);

%% Compensating for added mass in unembedded pile xy plane matrices

M_PILE_UNEMBED_xz=added_mass(Ca,D_out_pile,M_PILE_UNEMBED_xz1,N_pile_UNembed,L_pile_UNembed);

%% Compensating for added mass in unembedded pile in xz plane matrices

M_PILE_UNEMBED_yz=added_mass(Ca,D_out_pile,M_PILE_UNEMBED_yz1,N_pile_UNembed,L_pile_UNembed);

%% Parameter definition for Embedded Pile

M_axial_pile_embed=(240+110+348+147+m_pile*11.3); % Mass of rotor+nacelle+tower+transition piece+UNembedded pile [Ton]

L_pile_embed=24;% embedded Pile length [m]

k_soil=10*10^6; % Pa/m, taken from the DNV curves assuming angle of friction for sand is 30 degrees

% k_soil=0;

[M_PILE_EMBED_xz,K_PILE_EMBED_1_xz,M_PILE_EMBED_yz,K_PILE_EMBED_1_yz,M_PILE_EMBED_tor,K_PILE_EMBED_tor]=M_K_section(N_pile_embed,L_pile_embed,p,...
    A_pile,M_axial_pile_embed,m_pile,E,I_pile,alfa_pile,nu,J_pile,G);

%% Soil springs

% Adding springs to XZ plane stiffness matrices

K_PILE_EMBED_xz = linear_spring(K_PILE_EMBED_1_xz,N_pile_embed,L_pile_embed,k_soil);

% Adding springs to YZ plane stiffness matrices

K_PILE_EMBED_yz = linear_spring(K_PILE_EMBED_1_yz,N_pile_embed,L_pile_embed,k_soil);

%% Assembled System Matrices

% XZ PLANE matrix assembly

[M_system_xz,K_system_xz,N_xz,N_t,N_tp_aw,N_tp_uw,N_pile_uw]=Assemble_M_K(N_TWR,N_TP_AW,N_TP_UW,N_pile_UNembed,N_pile_embed,...
    M_TWR_xz,K_TWR_xz,M_TP_AW_xz,K_TP_AW_xz,M_TP_UW_xz,K_TP_UW_xz,M_PILE_UNEMBED_xz,K_PILE_UNEMBED_xz,...
    M_PILE_EMBED_xz,K_PILE_EMBED_xz);

% YZ PLANE matrix assembly

[M_system_yz,K_system_yz,~,~,~,~,~]=Assemble_M_K(N_TWR,N_TP_AW,N_TP_UW,N_pile_UNembed,N_pile_embed,...
    M_TWR_yz,K_TWR_yz,M_TP_AW_yz,K_TP_AW_yz,M_TP_UW_yz,K_TP_UW_yz,M_PILE_UNEMBED_yz,K_PILE_UNEMBED_yz,...
    M_PILE_EMBED_yz,K_PILE_EMBED_yz);

% segregating the nodes where hydro dynamic forces will act 

N_yz=N_xz+size(M_system_xz,1);

% Torsion matrix assembly
[M_system_tor,K_system_tor]=Assemble_tor(N_TWR,N_TP_AW,N_TP_UW,N_pile_UNembed,N_pile_embed,...
    M_TWR_tor,K_TWR_tor,M_TP_AW_tor,K_TP_AW_tor,M_TP_UW_tor,K_TP_UW_tor,M_PILE_UNEMBED_tor,K_PILE_UNEMBED_tor,...
    M_PILE_EMBED_tor,K_PILE_EMBED_tor);

%% Calculating natural frequency of the OWT structure

J_x=21775978.5;J_y=36748012;J_z=23079923.5; % Mass moments of Inertia including blade, hub and nacelle [Kgm^2] from NREL 5

% Adding the mass of rotor and Nacelle to the XZ Plane mass matrix

M_system_xz(1,1)=M_system_xz(1,1)+M_axial_TWR*1000;

M_system_xz(2,2)=M_system_xz(2,2)+J_y; % Adding rotary inertia to the mass matrix

% Adding the mass of rotor and Nacelle to the YZ Plane mass matrix

M_system_yz(1,1)=M_system_yz(1,1)+M_axial_TWR*1000;

M_system_yz(2,2)=M_system_yz(2,2)+J_x;

% Adding rotary inertia of the tower+nacelle to the global torsion matrix

M_system_tor(1,1)=M_system_tor(1,1)+J_z;

% Assembling all system matrices

[M_system,K_system] = system_assemble(M_system_xz,K_system_xz,M_system_yz,K_system_yz,M_system_tor,K_system_tor);

%% Natural Frequencies

% The last row and column have to be eliminated to assume non rotating conditions at the end

M_system_1=M_system(1:end-1,1:end-1);

K_system_1=K_system(1:end-1,1:end-1);

[hut,~,~]=nat_freq(M_system_1,K_system_1);

[~,~,v_xz]=nat_freq(M_system_xz,K_system_xz);

[~,~,v_yz]=nat_freq(M_system_yz,K_system_yz);

[~,~,v_tor]=nat_freq(M_system_tor(1:end-1,1:end-1),K_system_tor(1:end-1,1:end-1));

%% Modal parameters

n_modes=10; % number of participating modes

v=modal_vectors(v_xz,v_yz,v_tor,n_modes); % Reduced Modal matrix

M_mod=v'*M_system_1*v; % Reduced diagonal mass matrix

K_mod=v'*K_system_1*v; % Reduced diagonal stiffness matrix

%% Parameters for ODE solver

PRM.l_sub_tp=L_TP_UW/N_TP_UW; % length of one element for the submerged portion of the TP

PRM.Ntp_uw=N_TP_UW+1; % Number of nodes on submerged transition piece

PRM.l_sub_pile=L_pile_UNembed/N_pile_UNembed; % length of one element for the submerged portion of the Pile

PRM.Npile_uw=N_pile_UNembed+1; % Number of nodes on underwater pile

PRM.N_tp_uw=N_tp_uw; % Index numbers of submerged Transition piece nodes

PRM.N_pile_uw=N_pile_uw; % Index numbers of submerged Pile nodes

PRM.N_m=size(M_system_xz,1); % size of 1 plane mass or stiffness matrix

PRM.vel=vel; % total wave velocity in a direction theta 

PRM.F_mac=Force_MacCamy; % Maccamy force 

PRM.theta_waves=theta_waves; % angle of approach of waves

PRM.N=size(M_system_1,1); % total number of degrees of freedom after removing the last row and column
 
PRM.A=M_mod\K_mod;

PRM.B=M_mod\eye(size(M_mod,1),size(M_mod,1));

PRM.V_transp=v'; % for multiplication with force vector in ODE solver

PRM.V=v; % to find nodal displacements in terms of natural coordinates

PRM.F=zeros(PRM.N,1);

PRM.t=t1; % time for simulation

PRM.N_xz=N_xz; % index number of the submerged portions in global matrix for xz plane

PRM.N_yz=N_yz; % index number of the submerged portions in global matrix for yz plane

PRM.D_out=D_out_TP; % outer dia of the submerged portion

PRM.Cdw=Cdw; % Drag coefficient for waves

PRM.Cda=Cd_a; % Drag coefficient for air

PRM.n_modes=n_modes; % Number of participating modes

PRM.U_wind=red_fac*U_wind; % Wind velocties at the nodes exposed to air

PRM.t_wind=t_wind; % time series of wind loads

PRM.Acc_wind=red_fac*A_hub; % time series of wind accelerations

PRM.l_t=L_TWR_total/N_TWR; % length of one element for the tower

PRM.D_t=D_out_TWR; % Diameter of the tower

PRM.Nt=N_TWR+1; % Number of nodes on tower section

PRM.N_t=N_t; % Positions in the main vector where wind loads on the tower will act

PRM.N_tp_aw=N_tp_aw; % Positions in the main vector where wind loads on the above water transition piece will act

PRM.Ntp_aw=N_TP_AW+1; % Number of nodes on above water transition piece section

PRM.l_tp=L_TP_AW/N_TP_AW; % length of one element for the above water transition piece

% PRM.N_wind=N_wind; % start and end positions of nodes in air

PRM.theta_wind=theta_wind; % angle of approach for wind

PRM.rpm=rpm; % rotation speed of the turbine rotor

PRM.pitch=pitch; % pitch angle of the blades

PRM.overhang=overhang; % distance the rotor overhangs the tower axis

PRM.blade_prop=blade_prop; % matrix containing blade properties

PRM.Cd_blades=0.01; % drag coefficients for blade sections

% Initial conditions

q0=zeros(size(v,2),1);

q_dot_0=zeros(size(v,2),1);

IC=[q0;q_dot_0]; % defining initial conditions for natural coordinates

tspan=t1;% simulation period

h = waitbar(0,'Please wait...work in progress');

[t_out,q_out] = ode23(@(t,q)ode_func_lin(t,q,PRM),tspan,IC);

%% Getting responses

Q=(q_out(:,1:3*n_modes))';

X_l_out=v*Q;

X_l_out_dot=v*(q_out(:,3*n_modes+1:end))';

save(sprintf('x_l_out_%02d',ii),'X_l_out');

save(sprintf('x_l_out_dot_%02d',ii),'X_l_out_dot');

save(sprintf('t_%02d',ii),'t_out');

figure

plot(t_out,X_l_out(1,:));

title('Response of the rotor nacelle assembly in X-Z plane')

xlabel('t [s]')

ylabel('displacement [m]')

grid on

figure

plot(t_out,X_l_out(size(M_system_xz,1)+1,:));

title('Response of the rotor nacelle assembly in Y-Z plane')

xlabel('t [s]')

ylabel('displacement [m]')

grid on

figure

plot(t_out,X_l_out(141,:))


Y_data(:,1) = X_l_out(1,:);

Wind_data(:,1) = U_wind(1,:);

Wave_data(:,1) = Surf_el(1,:);



end
%% Animation

% H = anim(L_TWR_total,L_pile_total,L_TP_total,M_system_xz,rpm,L_blade,t1,X_l_out);
% 
% movie(H)
% 
% % close(h);
% 
% % Create AVI file.
% movie2avi(H, 'my_turb.avi', 'compression', 'None');

toc