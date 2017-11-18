function [F,T] = aero_force(W_global,W_dot_global,rpm,t,U,U_dot,blade_prop,pitch,Cd,d)

                         % Function calculates aerodynamic forces and torques resulting from wind on
                         % the turbine blades. It assumes that the blades rortate in anticlockwise direction 
                         % when facing them from fore to end direction. The positive direction for rotational velocity is considered clockwise. 
                         % It takes the follwing inputs
                         % W_global = wind velocity vector (3x1) as observed from earth frame of reference
                         % W_dot_global = wind accelaration vector (3x1) as observed from earth frame of reference
                         % rpm = rotational speed of the wind turbine [revolutions per minute]
                         % t = time [s]
                         % U = (6x1) vector of nodal displacements containing x,y,z displacements and rotations of the top node of the tower
                         % U_dot = (6x1) vector of nodal velocities containing x,y,z translational and rotational velocities of the top node of the tower
                         % pitch = blade pitch angle [degrees]
                         % blade_prop = (17x4) vector containing blade element length [m], twist angle [degrees], chord length [m] &
                         % element radius [m] in its 1st, 2nd, 3rd and 4th column respectively
                         % Cd = wind drag coefficient 0.01
                         % d = overhang of the rotor from tower centre line

 omega=rpm*2*pi/60; % rotational speed of the wind turbine [rad/s] (+ for clockwise and - for anti clockwise)
 
 azimuth_1=-omega*t; % position of Blade 1 at any moment of time 't'
 
 azimuth_2=-omega*t-2*pi/3; % position of Blade 2 at any moment of time 't'
 
 azimuth_3=-omega*t-4*pi/3; % position of Blade 3 at any moment of time 't'                        
                         
 %% Forces calculation
 
 F_drag_glob_1=zeros(3,1);F_drag_glob_2=zeros(3,1);F_drag_glob_3=zeros(3,1); % initialising drag force vectors for 3 blades
 
 T_drag_glob_1=zeros(3,1);T_drag_glob_2=zeros(3,1);T_drag_glob_3=zeros(3,1); % initialising drag torque vectors for 3 blades
 
 F_lift_glob_1=zeros(3,1);F_lift_glob_2=zeros(3,1);F_lift_glob_3=zeros(3,1); % initialising lift force vectors for 3 blades
 
 T_lift_glob_1=zeros(3,1);T_lift_glob_2=zeros(3,1);T_lift_glob_3=zeros(3,1); % initialising lift torque vectors for 3 blades
 
 F_inertia_glob_1=zeros(3,1);F_inertia_glob_2=zeros(3,1);F_inertia_glob_3=zeros(3,1); % initialising inertia force vectors for 3 blades
 
 T_inertia_glob_1=zeros(3,1);T_inertia_glob_2=zeros(3,1);T_inertia_glob_3=zeros(3,1); % initialising inertia torque vectors for 3 blades
 
%  Mf_1=zeros(3,6);Mf_2=zeros(3,6);Mf_3=zeros(3,6);
%  
%  Mt_1=zeros(3,6);Mt_2=zeros(3,6);Mt_3=zeros(3,6);
  
 for i=1:1:17
     
 % 1st blade
 
[R_Inv_1,K_fd_1,C_fd_1,K_td_1,C_td_1,K_fl_1,C_fl_1,K_tl_1,C_tl_1,R_xro_1,CC_1,DD_1,q3_1,...
    K_fi_1,C_fi_1,K_ti_1,C_ti_1,C_1] = rot_mat(azimuth_1,W_global,W_dot_global,blade_prop(i,4),omega,d,pitch,blade_prop(i,2));
 
F_drag_glob_1=F_drag_glob_1+0.5*1.225*blade_prop(i,1)*blade_prop(i,3)*Cd*(R_Inv_1*CC_1+K_fd_1*U+C_fd_1*U_dot); % Forces in global coordinate system
 
T_drag_glob_1=T_drag_glob_1+0.5*1.225*blade_prop(i,1)*blade_prop(i,3)*Cd*(R_xro_1*R_Inv_1*CC_1+K_td_1*U+C_td_1*U_dot); % Torque due to drag in global coordinate system
 
F_lift_glob_1= F_lift_glob_1+1.225*pi*blade_prop(i,1)*blade_prop(i,3)*(R_Inv_1*DD_1*q3_1+K_fl_1*U+C_fl_1*U_dot); % Lift force in global coordinate system
 
T_lift_glob_1=T_lift_glob_1+1.225*pi*blade_prop(i,1)*blade_prop(i,3)*(R_xro_1*R_Inv_1*DD_1*q3_1+K_tl_1*U+C_tl_1*U_dot); % Torque due to lift force in global coordinate system

F_inertia_glob_1= F_inertia_glob_1+1.225*pi*blade_prop(i,1)*(blade_prop(i,3))^2*...
    (C_1+K_fi_1*U+C_fi_1*U_dot); % Inertia force in global coordinate system

T_inertia_glob_1= T_inertia_glob_1+1.225*pi*blade_prop(i,1)*(blade_prop(i,3))^2*...
    (R_xro_1*C_1+K_ti_1*U+C_ti_1*U_dot); % torque due to inertia force in global coordinates
    
% Mf_1=Mf_1+1.225*pi*blade_prop(i,1)*(blade_prop(i,3))^2*M_fi_1;

% Mt_1=Mt_1+1.225*pi*blade_prop(i,1)*(blade_prop(i,3))^2*M_ti_1;

% 2nd blade
 
[R_Inv_2,K_fd_2,C_fd_2,K_td_2,C_td_2,K_fl_2,C_fl_2,K_tl_2,C_tl_2,R_xro_2,CC_2,DD_2,q3_2,...
    K_fi_2,C_fi_2,K_ti_2,C_ti_2,C_2] = rot_mat(azimuth_2,W_global,W_dot_global,blade_prop(i,4),omega,d,pitch,blade_prop(i,2));
 
F_drag_glob_2=F_drag_glob_2+0.5*1.225*blade_prop(i,1)*blade_prop(i,3)*Cd*(R_Inv_2*CC_2+K_fd_2*U+C_fd_2*U_dot); % Forces in global coordinate system
 
T_drag_glob_2=T_drag_glob_2+0.5*1.225*blade_prop(i,1)*blade_prop(i,3)*Cd*(R_xro_2*R_Inv_2*CC_2+K_td_2*U+C_td_2*U_dot); % Torque due to drag in global coordinate system
 
F_lift_glob_2= F_lift_glob_2+1.225*pi*blade_prop(i,1)*blade_prop(i,3)*(R_Inv_2*DD_2*q3_2+K_fl_2*U+C_fl_2*U_dot); % Lift force in global coordinate system
 
T_lift_glob_2=T_lift_glob_2+1.225*pi*blade_prop(i,1)*blade_prop(i,3)*(R_xro_2*R_Inv_2*DD_2*q3_2+K_tl_2*U+C_tl_2*U_dot); % Torque due to lift force in global coordinate system
 
F_inertia_glob_2= F_inertia_glob_2+1.225*pi*blade_prop(i,1)*(blade_prop(i,3))^2*...
    (C_2+K_fi_2*U+C_fi_2*U_dot); % Inertia force in global coordinate system

T_inertia_glob_2= T_inertia_glob_2+1.225*pi*blade_prop(i,1)*(blade_prop(i,3))^2*...
    (R_xro_2*C_2+K_ti_2*U+C_ti_2*U_dot); % torque due to inertia force in global coordinates
    
% Mf_2=Mf_2+1.225*pi*blade_prop(i,1)*(blade_prop(i,3))^2*M_fi_2;

% Mt_2=Mt_2+1.225*pi*blade_prop(i,1)*(blade_prop(i,3))^2*M_ti_2;

% 3rd blade
 
[R_Inv_3,K_fd_3,C_fd_3,K_td_3,C_td_3,K_fl_3,C_fl_3,K_tl_3,C_tl_3,R_xro_3,CC_3,DD_3,q3_3,...
     K_fi_3,C_fi_3,K_ti_3,C_ti_3,C_3] = rot_mat(azimuth_3,W_global,W_dot_global,blade_prop(i,4),omega,d,pitch,blade_prop(i,2));
 
F_drag_glob_3=F_drag_glob_3+0.5*1.225*blade_prop(i,1)*blade_prop(i,3)*Cd*(R_Inv_3*CC_3+K_fd_3*U+C_fd_3*U_dot); % Forces in global coordinate system
 
T_drag_glob_3=T_drag_glob_3+0.5*1.225*blade_prop(i,1)*blade_prop(i,3)*Cd*(R_xro_3*R_Inv_3*CC_3+K_td_3*U+C_td_3*U_dot); % Torque due to drag in global coordinate system
 
F_lift_glob_3= F_lift_glob_3+1.225*pi*blade_prop(i,1)*blade_prop(i,3)*(R_Inv_3*DD_3*q3_3+K_fl_3*U+C_fl_3*U_dot); % Lift force in global coordinate system
 
T_lift_glob_3=T_lift_glob_3+1.225*pi*blade_prop(i,1)*blade_prop(i,3)*(R_xro_3*R_Inv_3*DD_3*q3_3+K_tl_3*U+C_tl_3*U_dot); % Torque due to lift force in global coordinate system
 
F_inertia_glob_3= F_inertia_glob_3+1.225*pi*blade_prop(i,1)*(blade_prop(i,3))^2*...
    (C_3+K_fi_3*U+C_fi_3*U_dot); % Inertia force in global coordinate system

T_inertia_glob_3= T_inertia_glob_3+1.225*pi*blade_prop(i,1)*(blade_prop(i,3))^2*...
    (R_xro_3*C_3+K_ti_3*U+C_ti_3*U_dot); % torque due to inertia force in global coordinates
    
% Mf_3=Mf_3+1.225*pi*blade_prop(i,1)*(blade_prop(i,3))^2*M_fi_3;

% Mt_3=Mt_3+1.225*pi*blade_prop(i,1)*(blade_prop(i,3))^2*M_ti_3; 
 
 end
 
F_drag= F_drag_glob_1+F_drag_glob_2+F_drag_glob_3;
 
T_drag=T_drag_glob_1+T_drag_glob_2+T_drag_glob_3;
 
F_lift=F_lift_glob_1+F_lift_glob_2+F_lift_glob_3;

T_lift=T_lift_glob_1+T_lift_glob_2+T_lift_glob_3;

F_inertia=F_inertia_glob_1+F_inertia_glob_2+F_inertia_glob_3;

T_inertia=T_inertia_glob_1+T_inertia_glob_2+T_inertia_glob_3;

F=F_lift+F_drag+F_inertia;

T=T_lift+T_drag+T_inertia;

% Mf=Mf_1+Mf_2+Mf_3;
% 
% Mt=Mt_1+Mt_2+Mt_3;
                         
end

