function [F_tot,T_tot] = aero_force(W_global,W_dot_global,rpm,t,U,U_dot,theta,theta_dot,blade_prop,pitch,Cd,d)

                         % Function calculates aerodynamic forces and torques resulting from wind on
                         % the turbine blades. It assumes that the blades rotate in anticlockwise direction 
                         % when facing them from fore to end direction. The positive direction for rotational velocity is considered clockwise. 
                         % It takes the follwing inputs
                         % W_global = wind velocity vector (3x1) as observed from earth frame of reference
                         % W_dot_global = wind accelaration vector (3x1) as observed from earth frame of reference
                         % rpm = rotational speed of the wind turbine [revolutions per minute]
                         % t = time [s]
                         % U = (3x1) vector of nodal displacements containing x,y,z displacements of the top node of the tower
                         % U_dot = (3x1) vector of nodal velocities containing x,y,z velocities of the top node of the tower
                         % theta = (3x1) vector of nodal rotations of the top node of the tower
                         % theta_dot = (3x1) vector of nodal rotation velocities of the top node of the tower
                         % pitch = blade pitch angle [degrees]
                         % blade_prop = (17x4) vector containing blade element length [m], twist angle [degrees], chord length [m] &
                         % Cd= blade drag coefficient =[0.01]
                         % element radius [m] of the 17 sections of the blade in its 1st, 2nd, 3rd and 4th column respectively
                         % d = overhang of the rotor from tower centre line [m]

 omega=rpm*2*pi/60; % rotational speed of the wind turbine [rad/s] 
 
 azimuth_1=-omega*t; % position of Blade 1 at any moment of time 't'(+ for clockwise and - for anti clockwise)
 
 azimuth_2=-omega*t-2*pi/3; % position of Blade 2 at any moment of time 't'
 
 azimuth_3=-omega*t-4*pi/3; % position of Blade 3 at any moment of time 't'                        
                         
 RR=inv([1 -theta(3) theta(2,1);theta(3,1) 1 -theta(1,1);-theta(2,1) theta(1,1) 1]);
 
  %% Forces calculation
 
 F_drag_glob_1=zeros(3,1);F_drag_glob_2=zeros(3,1);F_drag_glob_3=zeros(3,1); % initialising drag force vectors for 3 blades
 
 T_drag_glob_1=zeros(3,1);T_drag_glob_2=zeros(3,1);T_drag_glob_3=zeros(3,1); % initialising drag torque vectors for 3 blades
 
 F_lift_glob_1=zeros(3,1);F_lift_glob_2=zeros(3,1);F_lift_glob_3=zeros(3,1);  % initialising lift force vectors for 3 blades
 
 T_lift_glob_1=zeros(3,1);T_lift_glob_2=zeros(3,1);T_lift_glob_3=zeros(3,1);  % initialising lift torque vectors for 3 blades
 
 F_inertia_glob_1=zeros(3,1);F_inertia_glob_2=zeros(3,1);F_inertia_glob_3=zeros(3,1);  % initialising inertia force vectors for 3 blades
 
 T_inertia_glob_1=zeros(3,1);T_inertia_glob_2=zeros(3,1);T_inertia_glob_3=zeros(3,1);  % initialising inertia torque vectors for 3 blades
 
%  Ft2_1=zeros(3,3);Ft2_2=zeros(3,3);Ft2_3=zeros(3,3);
% 
%  Tt1_1=zeros(3,1);Tt1_2=zeros(3,1);Tt1_3=zeros(3,1);
%  
%  Tt2_1=zeros(3,3);Tt2_2=zeros(3,3);Tt2_3=zeros(3,3);
 
 p=[1 0 0;0 1 0;0 0 0]; % a matrix for making the last component of a vector zero
 
 % Calculating forces and torques and summing over the length of each blade (each blade divided into 17 parts)
 
 for i=1:1:17 
     
 beta_tp=[-sind(pitch+blade_prop(i,2)) cosd(pitch+blade_prop(i,2)) 0];
     
  % 1st blade
 
 [R_vec1,U_dot_local_1,U_dot_dot_local_1,R_glob2loc_1,R_dot_glob2loc_1,R1] = rot_mat(azimuth_1,U,U_dot,theta,theta_dot,d,blade_prop(i,4),omega); % obtaining local velocities and rotation vector
 
 W_local_1=R_glob2loc_1*W_global+[omega*blade_prop(i,4);0;0]; % relative velocity of wind w.r.t the ith blade section in local frame of reference of 1st blade
 
 W_dot_local_1=R_dot_glob2loc_1*W_global+R_glob2loc_1*W_dot_global; % acceleration of wind from blade reference frame
 
 V_local_1=p*(W_local_1-U_dot_local_1); % Net relative velocity vector including structural velocity
 
 V_dot_local_1=p*(W_dot_local_1-U_dot_dot_local_1); % Net relative acceleration vector including structural acceleration
 
 f_drag_glob_1=0.5*1.225*blade_prop(i,1)*blade_prop(i,3)*Cd*(R_glob2loc_1\V_local_1)*sqrt(sum(V_local_1.*V_local_1));
 
 F_drag_glob_1=F_drag_glob_1+f_drag_glob_1;
 
 T_drag_glob_1=T_drag_glob_1+cross(R_vec1,f_drag_glob_1); % Drag Torque in global coordinate system
 
 f_lift_glob_1=1.225*pi*blade_prop(i,1)*blade_prop(i,3)*beta_tp*...
     V_local_1*(R_glob2loc_1\cross([0;0;1],V_local_1)); % Lift force in global coordinate system
 
 F_lift_glob_1= F_lift_glob_1+f_lift_glob_1;
 
 T_lift_glob_1= T_lift_glob_1+cross(R_vec1,f_lift_glob_1); % Lift Torque in global coordinate system
 
 f_inertia_glob_1=1.225*pi*blade_prop(i,1)*(blade_prop(i,3))^2*(R_glob2loc_1\V_dot_local_1); % Inertia force in global coordinate system
 
 F_inertia_glob_1=F_inertia_glob_1+f_inertia_glob_1;
 
 T_inertia_glob_1=T_inertia_glob_1+cross(R_vec1,f_inertia_glob_1); % Inertia Torque in global coordinate system
 
%  Ft2_1=Ft2_1+1.225*pi*blade_prop(i,1)*(blade_prop(i,3))^2*R1;
%  
%  Tt1_1=Tt1_1+1.225*pi*blade_prop(i,1)*(blade_prop(i,3))^2*R_vec1;
%  
%  R11=RR*R1;
%  
%  Tt2_1=Tt2_1+1.225*pi*blade_prop(i,1)*(blade_prop(i,3))^2*...
%      [(R11(3,1)*R_vec1(2)-R11(2,1)*R_vec1(3)) (R11(3,2)*R_vec1(2)-R11(2,2)*R_vec1(3)) (R11(3,3)*R_vec1(2)-R11(2,3)*R_vec1(3));...
%      (R11(1,1)*R_vec1(3)-R11(3,1)*R_vec1(1)) (R11(1,2)*R_vec1(3)-R11(3,2)*R_vec1(1)) (R11(1,3)*R_vec1(3)-R11(3,3)*R_vec1(1));...
%      (R11(2,1)*R_vec1(1)-R11(1,1)*R_vec1(2)) (R11(2,2)*R_vec1(1)-R11(1,2)*R_vec1(2)) (R11(2,3)*R_vec1(1)-R11(1,3)*R_vec1(2))];
 
 
 % 2nd blade
 
 [R_vec2,U_dot_local_2,U_dot_dot_local_2,R_glob2loc_2,R_dot_glob2loc_2,R2] = rot_mat(azimuth_2,U,U_dot,theta,theta_dot,d,blade_prop(i,4),omega); % obtaining local velocities and rotation vector
 
 W_local_2=R_glob2loc_2*W_global+[omega*blade_prop(i,4);0;0]; % relative velocity of wind w.r.t the ith blade section in local frame of reference of 2nd blade
 
 W_dot_local_2=R_dot_glob2loc_2*W_global+R_glob2loc_2*W_dot_global; % acceleration of wind from blade
 
 V_local_2=p*(W_local_2-U_dot_local_2); % Net relative velocity vector
 
 V_dot_local_2=p*(W_dot_local_2-U_dot_dot_local_2); % Net relative acceleration vector including structural acceleration
 
 f_drag_glob_2=0.5*1.225*blade_prop(i,1)*blade_prop(i,3)*Cd*(R_glob2loc_2\V_local_2)*sqrt(sum(V_local_2.*V_local_2)); % Forces in global coordinate system
  
 F_drag_glob_2=F_drag_glob_2+f_drag_glob_2;
 
 T_drag_glob_2=T_drag_glob_2+cross(R_vec2,f_drag_glob_2); % Torque in global coordinate system
 
 f_lift_glob_2=1.225*pi*blade_prop(i,1)*blade_prop(i,3)*beta_tp*...
     V_local_2*(R_glob2loc_2\cross([0;0;1],V_local_2)); % Lift force in global coordinate system
 
 F_lift_glob_2=F_lift_glob_2+f_lift_glob_2;
 
 T_lift_glob_2= T_lift_glob_2+cross(R_vec2,f_lift_glob_2); % Torque in global coordinate system
 
 f_inertia_glob_2=1.225*pi*blade_prop(i,1)*(blade_prop(i,3))^2*(R_glob2loc_2\V_dot_local_2); % Inertia force in global coordinate system
 
 F_inertia_glob_2=F_inertia_glob_2+f_inertia_glob_2;
 
 T_inertia_glob_2=T_inertia_glob_2+cross(R_vec2,f_inertia_glob_2); % Torque in global coordinate system
 
%  Ft2_2=Ft2_2+1.225*pi*blade_prop(i,1)*(blade_prop(i,3))^2*R2;
%  
%  Tt1_2=Tt1_2+1.225*pi*blade_prop(i,1)*(blade_prop(i,3))^2*R_vec2;
%  
%  R22=RR*R2;
%  
%  Tt2_2=Tt2_2+1.225*pi*blade_prop(i,1)*(blade_prop(i,3))^2*...
%      [(R22(3,1)*R_vec2(2)-R22(2,1)*R_vec2(3)) (R22(3,2)*R_vec2(2)-R22(2,2)*R_vec2(3)) (R22(3,3)*R_vec2(2)-R22(2,3)*R_vec2(3));...
%      (R22(1,1)*R_vec2(3)-R22(3,1)*R_vec2(1)) (R22(1,2)*R_vec2(3)-R22(3,2)*R_vec2(1)) (R22(1,3)*R_vec2(3)-R22(3,3)*R_vec2(1));...
%      (R22(2,1)*R_vec2(1)-R22(1,1)*R_vec2(2)) (R22(2,2)*R_vec2(1)-R22(1,2)*R_vec2(2)) (R22(2,3)*R_vec2(1)-R22(1,3)*R_vec2(2))];
 
 % 3rd blade
 
 [R_vec3,U_dot_local_3,U_dot_dot_local_3,R_glob2loc_3,R_dot_glob2loc_3,R3] = rot_mat(azimuth_3,U,U_dot,theta,theta_dot,d,blade_prop(i,4),omega); % obtaining local velocities and rotation vector
 
 W_local_3=R_glob2loc_3*W_global+[omega*blade_prop(i,4);0;0]; % relative velocity of wind w.r.t the ith blade section in local frame of reference of 3rd blade
 
 W_dot_local_3=R_dot_glob2loc_3*W_global+R_glob2loc_3*W_dot_global; % acceleration of wind from blade
 
 V_local_3=p*(W_local_3-U_dot_local_3); % Net relative velocity vector
 
 V_dot_local_3=p*(W_dot_local_3-U_dot_dot_local_3); % Net relative acceleration vector including structural acceleration
 
 f_drag_glob_3=0.5*1.225*blade_prop(i,1)*blade_prop(i,3)*Cd*(R_glob2loc_3\V_local_3)*sqrt(sum(V_local_3.*V_local_3)); % Forces in global coordinate system
 
 F_drag_glob_3=F_drag_glob_3+f_drag_glob_3;
 
 T_drag_glob_3=T_drag_glob_3+cross(R_vec3,f_drag_glob_3); % Torque in global coordinate system
 
 f_lift_glob_3=1.225*pi*blade_prop(i,1)*blade_prop(i,3)*beta_tp*...
     V_local_3*(R_glob2loc_3\cross([0;0;1],V_local_3)); % Lift force in global coordinate system
 
 F_lift_glob_3= F_lift_glob_3+f_lift_glob_3;
 
 T_lift_glob_3= T_lift_glob_3+cross(R_vec3,f_lift_glob_3); % Lift Torque in global coordinate system
 
 f_inertia_glob_3=1.225*pi*blade_prop(i,1)*(blade_prop(i,3))^2*(R_glob2loc_3\V_dot_local_3); % Inertia force in global coordinate system
 
 F_inertia_glob_3=F_inertia_glob_3+f_inertia_glob_3;
 
 T_inertia_glob_3=T_inertia_glob_3+cross(R_vec3,f_inertia_glob_3); % Inertia Torque in global coordinate system
 
%  Ft2_3=Ft2_3+1.225*pi*blade_prop(i,1)*(blade_prop(i,3))^2*R3;
%   
%  Tt1_3=Tt1_3+1.225*pi*blade_prop(i,1)*(blade_prop(i,3))^2*R_vec3;
%  
%  R33=RR*R3;
%  
%  Tt2_3=Tt2_3+1.225*pi*blade_prop(i,1)*(blade_prop(i,3))^2*...
%      [(R33(3,1)*R_vec3(2)-R33(2,1)*R_vec3(3)) (R33(3,2)*R_vec3(2)-R33(2,2)*R_vec3(3)) (R33(3,3)*R_vec3(2)-R33(2,3)*R_vec3(3));...
%      (R33(1,1)*R_vec3(3)-R33(3,1)*R_vec3(1)) (R33(1,2)*R_vec3(3)-R33(3,2)*R_vec3(1)) (R33(1,3)*R_vec3(3)-R33(3,3)*R_vec3(1));...
%      (R33(2,1)*R_vec3(1)-R33(1,1)*R_vec3(2)) (R33(2,2)*R_vec3(1)-R33(1,2)*R_vec3(2)) (R33(2,3)*R_vec3(1)-R33(1,3)*R_vec3(2))];
 
 
 end
 
 F_drag= F_drag_glob_1+ F_drag_glob_2+ F_drag_glob_3;
 
 T_drag=T_drag_glob_1+T_drag_glob_2+T_drag_glob_3;
 
 F_lift=F_lift_glob_1+F_lift_glob_2+F_lift_glob_3;
 
 T_lift=T_lift_glob_1+T_lift_glob_2+T_lift_glob_3;
 
 F_inertia=F_inertia_glob_1+F_inertia_glob_2+F_inertia_glob_3;
 
 T_inertia=T_inertia_glob_1+T_inertia_glob_2+T_inertia_glob_3;
 
 F_tot=F_drag+F_lift+F_inertia;
 
 T_tot=T_drag+T_lift+T_inertia;
 
%  Ft1=-3*sum(1.225*pi.*blade_prop(i,1).*(blade_prop(i,3)).^2)*RR;
%  
%  Ft2=RR*(Ft2_1+Ft2_2+Ft2_3);
%  
%  tt1=Tt1_1+Tt1_2+Tt1_3;
%  
%  TT1=-[(RR(3,1)*tt1(2)-RR(2,1)*tt1(3)) (RR(3,2)*tt1(2)-RR(2,2)*tt1(3)) (RR(3,3)*tt1(2)-RR(2,3)*tt1(3));...
%      (RR(1,1)*tt1(3)-RR(3,1)*tt1(1)) (RR(1,2)*tt1(3)-RR(3,2)*tt1(1)) (RR(1,3)*tt1(3)-RR(3,3)*tt1(1));...
%      (RR(2,1)*tt1(1)-RR(1,1)*tt1(2)) (RR(2,2)*tt1(1)-RR(1,2)*tt1(2)) (RR(2,3)*tt1(1)-RR(1,3)*tt1(2))];
%  
%  TT2=Tt2_1+Tt2_2+Tt2_3;
end

