function [R_Inv,K_fd,C_fd,K_td,C_td,K_fl,C_fl,K_tl,C_tl,D,CC,DD,q0,...
    K_fi,C_fi,K_ti,C_ti,q3] = rot_mat(azimuth,W_glob,W_glob_dot,r,omega,d,pitch,twist)

                            % rot_mat is a function that returns U_dot_local which is (3x1) structural
                            % velocities in local frame of reference,
                            % a matrix R_glob2loc for transforming global variables to
                            % locals and vice versa and R_vec that gives
                            % position vector of a blade segment in global
                            % reference frame for calculating torque
                            % It takes the following inputs
                            % azimuths = spatial position of the turbine blade at any time "t"
                            % U =(3x1) vector of nodal displacements containing x,y,z displacements of the top node of the tower
                            % U_dot=(3x1) vector of nodal velocities containing x,y,z velocities of the top node of the tower
                            % theta = (3x1) vector containing rotations of the top most node of the tower 
                            % theta_dot=(3x1) vector of nodal rotation velocities of the top node of the tower
                            % d = overhang of the rotor [m] 
                            % r = section radius [m]
                            % omega = rotational velocity [rad/s]
                            
%% Creating matrix for transforming global frame of reference to local frame of reference

A=[-sin(azimuth) 0 -cos(azimuth);0 1 0;cos(azimuth) 0 -sin(azimuth)]; % R(X->Xjo(t))

% B=[0 -U(6,1) U(5,1);U(6,1) 0 -U(4,1);-U(5,1) U(4,1) 0]; % R_theta(t)

C=[omega*cos(azimuth) 0 -omega*sin(azimuth);0 0 0;omega*sin(azimuth) 0 omega*cos(azimuth)] ; % (R(X->Xjo(t)) time derivative)

C_dot=[omega^2*sin(azimuth) 0 omega^2*cos(azimuth);0 0 0;-omega^2*cos(azimuth) 0 omega^2*sin(azimuth)];% (R(X->Xjo(t)) double time derivative)

D=[0 r*sin(azimuth) -d;-r*sin(azimuth) 0 -r*cos(azimuth);d r*cos(azimuth) 0]; % R(X_ro(t))

E=[0 -omega*r*cos(azimuth) 0;omega*r*cos(azimuth) 0 -omega*r*sin(azimuth);0 omega*r*sin(azimuth) 0]; % R(X_ro(t)) time derivative

E_dot=[0 -r*omega^2*sin(azimuth) 0;r*omega^2*sin(azimuth) 0 r*omega^2*cos(azimuth);0 -r*omega^2*cos(azimuth) 0];% R(X_ro(t)) double time derivative

% R_glob2loc=A*(eye(3,3)+B); % R(X->Xj(t)) = R(X->Xjo(t))(eye(3,3)+R_theta(t)) rotational matrix

p=[1 0 0;0 1 0;0 0 0]; % a matrix for making the last component of a vector zero

W_loc=p*(A*W_glob+[omega*r;0;0]); % relative velocity of wind w.r.t the ith blade section in local frame of reference of blade including the effect of tangential velocity

W_loc_dot=p*(C*W_glob+A*W_glob_dot);
%% Drag parameters

R_Inv=A\eye(3,3); % Inverse of A (R(X->Xjo(t))) matrix

AA=sqrt(sum(W_loc.*W_loc));% |V| term

BB=(W_loc*W_loc'/AA);% VV'/|V| term

CC=AA*W_loc; % |V|V terms

Qw=[0 -W_glob(3,1) W_glob(2,1);W_glob(3,1) 0 -W_glob(1,1);-W_glob(2,1) W_glob(1,1) 0];

Qw_dot=[0 -W_glob_dot(3,1) W_glob_dot(2,1);W_glob_dot(3,1) 0 -W_glob_dot(1,1);-W_glob_dot(2,1) W_glob_dot(1,1) 0];

q1=R_Inv*W_loc;

Q1=[0 -q1(3,1) q1(2,1);q1(3,1) 0 -q1(1,1);-q1(2,1) q1(1,1) 0];

aa=-R_Inv*BB*C-R_Inv*AA*C;

bb=-R_Inv*BB*Qw-R_Inv*AA*(Qw-A*Q1)+...
    R_Inv*BB*(C*D+A*E)+AA*(R_Inv*C*D+E);

cc=-R_Inv*BB*A-AA*eye(3,3);

dd=R_Inv*BB*A*D+AA*D;

K_fd=[aa bb];

C_fd=[cc dd];

K_td=D*K_fd+AA*[-Q1 Q1*D];

C_td=D*C_fd;

%% Lift parameters

beta_tp=[-sind(pitch+twist) cosd(pitch+twist) 0];

DD=beta_tp*W_loc;

e_r=[0;0;1]; % unit vector in radial direction

Qr=[0 -1 0;1 0 0;0 0 0];

q0=cross(e_r,W_loc);

q2=R_Inv*q0;

Q2=[0 -q2(3,1) q2(2,1);q2(3,1) 0 -q2(1,1);-q2(2,1) q2(1,1) 0];

ee=-R_Inv*DD*Qr*C-R_Inv*(q0*beta_tp)*C;

ff=-R_Inv*(DD*Qr*Qw+(q0*beta_tp)*Qw-DD*A*Q2)+...
    R_Inv*(DD*Qr+q0*beta_tp)*(C*D+A*E);

gg=-R_Inv*DD*Qr*A-R_Inv*(q0*beta_tp)*A;

hh=R_Inv*DD*Qr*A*D+R_Inv*(q0*beta_tp)*A*D;

K_fl=[ee ff];

C_fl=[gg hh];

K_tl=D*K_fl+DD*[-Q2 Q2*D];

C_tl=D*C_fl;

%%  Inertia parameters

q3=R_Inv*W_loc_dot;

Q3=[0 -q3(3,1) q3(2,1);q3(3,1) 0 -q3(1,1);-q3(2,1) q3(1,1) 0];

ii=-R_Inv*C_dot;

jj=R_Inv*(C_dot*D+2*C*E+A*E_dot-(Qw_dot-A*Q3));

kk=-2*R_Inv*C;

ll=R_Inv*(2*C*D+2*A*E-Qw);

K_fi=[ii jj];

C_fi=[kk ll];

% M_fi=[eye(3,3);D]';

K_ti=D*K_fi+[-Q3 Q3*D];

C_ti=D*C_fi;

% M_ti=D*M_fi;
end

