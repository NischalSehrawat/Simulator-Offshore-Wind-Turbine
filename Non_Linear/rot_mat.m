function [R_vec,U_dot_local,U_dot_dot_local,R_glob2loc,R_dot_glob2loc,D] = rot_mat(azimuth,U,U_dot,theta,theta_dot,d,r,omega)

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

B=[0 -theta(3,1) theta(2,1);theta(3,1) 0 -theta(1,1);-theta(2,1) theta(1,1) 0]; % R_theta(t)

B_dot=[0 -theta_dot(3,1) theta_dot(2,1);theta_dot(3,1) 0 -theta_dot(1,1);-theta_dot(2,1) theta_dot(1,1) 0]; % R_theta(t) time derivative

C=[omega*cos(azimuth) 0 -omega*sin(azimuth);0 0 0;omega*sin(azimuth) 0 omega*cos(azimuth)] ; % (R(X->Xjo(t)) time derivative)

C_dot=[omega^2*sin(azimuth) 0 omega^2*cos(azimuth);0 0 0;-omega^2*cos(azimuth) 0 omega^2*sin(azimuth)];% (R(X->Xjo(t)) double time derivative)

D=[0 r*sin(azimuth) -d;-r*sin(azimuth) 0 -r*cos(azimuth);d r*cos(azimuth) 0]; % R(X_ro(t))

E=[0 -omega*r*cos(azimuth) 0;omega*r*cos(azimuth) 0 -omega*r*sin(azimuth);0 omega*r*sin(azimuth) 0]; % R(X_ro(t)) time derivative

E_dot=[0 -r*omega^2*sin(azimuth) 0;r*omega^2*sin(azimuth) 0 r*omega^2*cos(azimuth);0 -r*omega^2*cos(azimuth) 0];% R(X_ro(t)) double time derivative

R_glob2loc=A*(eye(3,3)+B); % R(X->Xj(t)) = R(X->Xjo(t))(eye(3,3)+R_theta(t)) rotational matrix

R_dot_glob2loc=C*(eye(3,3)+B)+A*B_dot; % Time derivative of R(X->Xj(t))

%% Structural Velocity calculations in local frame of reference

T1=C*U;

T2=A*U_dot;

T3=(C*D+A*E)*theta;

T4=A*D*theta_dot;

U_dot_local=T1+T2-T3-T4;

%% Structural Acceleration calculations in local frame of reference

U_dot_dot_local=C_dot*U+2*C*U_dot-...
    (C_dot*D+2*C*E+A*E_dot)*theta-...
    2*(C*D+A*E)*theta_dot;

%% Position vector for an element in global frame of reference

X_ro=[r*cos(azimuth);-d;-r*sin(azimuth)];

R_vec=U+X_ro+cross(theta,X_ro);

end

