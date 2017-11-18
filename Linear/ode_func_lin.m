function dq=ode_func_lin(t,q,PRM)

% q is a vector of size 6*n_modes.The 1st 3*n_modes rows contain the
% natural coordinate displacements and the next 3*n_modes rows contain the natural coordinate velocities

% dq is a vector of size 6*n_modes.The 1st 3*n_modes rows contain the
% natural coordinate velocities and the next 3*n_modes rows contain the natural coordinate accelerations

Q=q(1:3*PRM.n_modes); % Defining 1st 3*n_modes rows as natural coordinate displacements

Q_dot=q(3*PRM.n_modes+1:6*PRM.n_modes); % Defining next N rows as natural coordinate velocities

X=PRM.V*Q; % Nodal displacementvvs expressed in terms of natural coordinates

X_dot=PRM.V*Q_dot; % Nodal velocities expressed in terms of natural coordinates

%% Hydrodynamic Loading

n1=size(PRM.vel,1); % number of rows for wave velocities

for i=1:1:n1
    
wave_vel(i,1)=interp1(PRM.t,PRM.vel(i,:),t); % reading total wave velocities from table

f_mac(i,1)=interp1(PRM.t,PRM.F_mac(i,:),t); % reading MacCamy force from table

end

wave_vel_xz=wave_vel.*cosd(PRM.theta_waves); % wave velocities in xz plane

wave_vel_yz=wave_vel.*sind(PRM.theta_waves); % wave velocities in yz plane

f_mac_xz=f_mac.*cosd(PRM.theta_waves); % MacCamy force in xz plane

f_mac_yz=f_mac.*sind(PRM.theta_waves);  % MacCamy force in yz plane

%% Arranging MacCamy forces for combined underwater TP and submerged pile for final assembly in Force vector

F1=PRM.F;

F1(PRM.N_xz(1):2:PRM.N_xz(2),1)=f_mac_xz;

F2=PRM.F;

F2(PRM.N_yz(1):2:PRM.N_yz(2),1)=f_mac_yz;

%% Drag forces on under water TRANSITION PIECE in both XZ and YZ planes

% Forces in xz plane

f_dr_uwtp_xz=0.5*1025*PRM.Cdw*PRM.D_out.*abs(wave_vel(1:PRM.Ntp_uw,1)).*(wave_vel_xz(1:PRM.Ntp_uw,1))*PRM.l_sub_tp;

f_dr_uwtp_xz(1,1)=0.5*f_dr_uwtp_xz(1,1);f_dr_uwtp_xz(end,1)=f_dr_uwtp_xz(end,1)*0.5;

% assembling in force vector

F3=PRM.F;

F3(PRM.N_tp_uw(1):2:PRM.N_tp_uw(2),1)=f_dr_uwtp_xz;

% Forces in yz plane

f_dr_uwtp_yz=0.5*1025*PRM.Cdw*PRM.D_out*abs(wave_vel(1:PRM.Ntp_uw,1)).*(wave_vel_yz(1:PRM.Ntp_uw,1))*PRM.l_sub_tp;

f_dr_uwtp_yz(1,1)=0.5*f_dr_uwtp_yz(1,1);f_dr_uwtp_yz(end,1)=f_dr_uwtp_yz(end,1)*0.5;

% assembling in force vector

F4=PRM.F;

F4(PRM.N_tp_uw(1)+PRM.N_m:2:PRM.N_tp_uw(2)+PRM.N_m,1)=f_dr_uwtp_yz;

%% Drag forces on submerged PILE in both XZ and YZ planes

% drag forces in xz plane

f_dr_pileuw_xz=0.5*1025*PRM.Cdw*PRM.D_out*abs(wave_vel(PRM.Ntp_uw:end,1)).*(wave_vel_xz(PRM.Ntp_uw:end,1))*PRM.l_sub_pile;

f_dr_pileuw_xz(1,1)=0.5*f_dr_pileuw_xz(1,1);f_dr_pileuw_xz(end,1)=0.5*f_dr_pileuw_xz(end,1);

% assembling in force vector

F5=PRM.F;

F5(PRM.N_pile_uw(1):2:PRM.N_pile_uw(2),1)=f_dr_pileuw_xz;

% drag forces in yz plane

f_dr_pileuw_yz=0.5*1025*PRM.Cdw*PRM.D_out*abs(wave_vel(PRM.Ntp_uw:end,1)).*(wave_vel_yz(PRM.Ntp_uw:end,1))*PRM.l_sub_pile;

f_dr_pileuw_yz(1,1)=0.5*f_dr_pileuw_yz(1,1);f_dr_pileuw_yz(end,1)=0.5*f_dr_pileuw_yz(end,1);

% assembling in force vector

F6=PRM.F;

F6(PRM.N_pile_uw(1)+PRM.N_m:2:PRM.N_pile_uw(2)+PRM.N_m,1)=f_dr_pileuw_yz;

%% Aerodynamic loading

n2=size(PRM.U_wind,1); % number of rows for wind velocities

for i=1:1:n2

    wind_vel(i,1)=interp1(PRM.t_wind,PRM.U_wind(i,:),t); % reading values of wind velocties from table 
end

wind_vel_xz=wind_vel.*cosd(PRM.theta_wind); % wave velocities in xz plane

wind_vel_yz=wind_vel.*sind(PRM.theta_wind); % wave velocities in yz plane

% Calculating wind drag forces on tower XZ plane

f_wind_t_xz=0.5*1.225*PRM.Cda*PRM.D_t*wind_vel(1:PRM.Nt,1).*(wind_vel_xz(1:PRM.Nt,1))*PRM.l_t;

f_wind_t_xz(1,1)=0.5*f_wind_t_xz(1,1);f_wind_t_xz(end,1)=0.5*f_wind_t_xz(end,1);

% assembling in force vector

F7=PRM.F;

F7(PRM.N_t(1):2:PRM.N_t(2))=f_wind_t_xz;

% YZ plane wind drag on tower

f_wind_t_yz=0.5*1.225*PRM.Cda*PRM.D_t*wind_vel(1:PRM.Nt,1).*(wind_vel_yz(1:PRM.Nt,1))*PRM.l_t;

f_wind_t_yz(1,1)=0.5*f_wind_t_yz(1,1);f_wind_t_yz(end,1)=0.5*f_wind_t_yz(end,1);

% assembling in force vector

F8=PRM.F;

F8(PRM.N_t(1)+PRM.N_m:2:PRM.N_t(2)+PRM.N_m)=f_wind_t_yz;

% Calculating wind drag forces on above water transition piece in xz plane

f_wind_tp_xz=0.5*1.225*PRM.Cda*PRM.D_out*wind_vel(PRM.Nt:end,1).*(wind_vel_xz(PRM.Nt:end,1))*PRM.l_tp;

f_wind_tp_xz(1,1)=0.5*f_wind_tp_xz(1,1);f_wind_tp_xz(end,1)=0.5*f_wind_tp_xz(end,1);

% Assembling the wind drag forces from transition piece in xz plane in force vector

F9=PRM.F; % in this, the tower wind forces will be substituted

F9(PRM.N_tp_aw(1):2:PRM.N_tp_aw(2))=f_wind_tp_xz;

% Calculating wind drag forces on above water transition piece in YZ 

f_wind_tp_yz=0.5*1.225*PRM.Cda*PRM.D_out*wind_vel(PRM.Nt:end,1).*(wind_vel_yz(PRM.Nt:end,1))*PRM.l_tp;

f_wind_tp_yz(1,1)=0.5*f_wind_tp_yz(1,1);f_wind_tp_yz(end,1)=0.5*f_wind_tp_yz(end,1);

% Assembling the wind drag forces from transition piece in yz plane in force vector

F10=PRM.F; % in this, the tower wind forces will be substituted

F10(PRM.N_tp_aw(1)+PRM.N_m:2:PRM.N_tp_aw(2)+PRM.N_m)=f_wind_tp_yz;

%% Blade loading

U=[X(1,1);X(PRM.N_m+1,1);0;X(PRM.N_m+2,1);X(2,1);-X(2*PRM.N_m+1,1)];% nodal displacement vector of the top most node

U_dot=[X_dot(1,1);X_dot(PRM.N_m+1,1);0;X_dot(PRM.N_m+2,1);X_dot(2,1);-X_dot(2*PRM.N_m+1,1)]; % nodal velocity vector of the top most node

% U_dd=[PRM.V(1,:);PRM.V(PRM.N_m+1,:);zeros(1,3*PRM.n_modes);PRM.V(PRM.N_m+2,:);PRM.V(2,:);PRM.V(2*PRM.N_m+1,:)]; % translatory terms for inertia 

W=[wind_vel_xz(1,1);wind_vel_yz(1,1);0]; % global wind velocity vector 
    
w_dot=interp1(PRM.t_wind,PRM.Acc_wind(1,:),t);

W_dot=[w_dot*cosd(PRM.theta_wind);w_dot*sind(PRM.theta_wind);0]; % global wind acceleration vector
    
[F_tot,T_tot] = aero_force(W,W_dot,PRM.rpm,t,U,U_dot,PRM.blade_prop,PRM.pitch,PRM.Cd_blades,PRM.overhang);

% F_blade_vec=Mf*U_dd;
% 
% T_blade_vec=Mt*U_dd;
% 
% F_in=zeros(PRM.N,3*PRM.n_modes); % air inertia matrix (nxr) will be subtracted from reduced mass matrix
% 
% F_in(1,:)=F_blade_vec(1,:); % inertia terms in X direction 
% 
% F_in(2,:)=T_blade_vec(2,:);% inertia terms bout X axis
% 
% F_in(PRM.N_m+1,:)=F_blade_vec(2,:); % inertia term in Y direction 
% 
% F_in(PRM.N_m+2,:)=T_blade_vec(1,:) ;% inertia terms about Y axis
% 
% F_in(2*PRM.N_m+1,:)=T_blade_vec(3,:); % inertia term about Z axis
% 
% M_in=PRM.V_transp*F_in; % aerodynamic inertia matrix that is to be subtracted from main reduced mass matrix

% Assembling all the force vectors

F=F1+F2+F3+F4+F5+F6+F7+F8+F9+F10;

% Adding blade aerodynamic forces to the the force vector

F(1,1)=F(1,1)+F_tot(1,1); % Blade Force in XZ plane

F(2,1)=F(2,1)+T_tot(2,1); % Torque in XZ plane

F(PRM.N_m+1,1)=F(PRM.N_m+1,1)+F_tot(2,1); % Blade force in YZ plane

F(PRM.N_m+2,1)=F(PRM.N_m+2,1)+T_tot(1,1); % torque in YZ plane

F(2*PRM.N_m+1,1)=F(2*PRM.N_m+1,1)-T_tot(3,1); % torque about Z axis

% calculating natural coordinate accelerations 

% PRM.A=PRM.A-M_in;

Q_dot_dot=-PRM.A*Q+PRM.B*(PRM.V_transp*F);

dq=[Q_dot;Q_dot_dot];

waitbar(t/PRM.t(1,end));

end

