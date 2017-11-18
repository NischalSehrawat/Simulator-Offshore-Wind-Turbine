function H = anim(L_TWR_total,L_pile_total,L_TP_total,M_system_xz,rpm,L_blade,t1,X_out)


L1=L_TWR_total+L_pile_total+L_TP_total;

n11=size(M_system_xz(1:2:end,:),1);

L=flipud((linspace(0,L1,n11))');

omega_rot=2*pi*rpm/60;

for i=1:1:length(t1)

X11=30*X_out(1:2:size(M_system_xz,1),i); % X coordinates

X22=30*X_out(size(M_system_xz,1)+1:2:2*size(M_system_xz,1),i); % Y coordinates

plot3(X11,X22,L,'LineWidth',7);hold on; % plotting OWT in 3D

title('Real time displacement of an OWT')

xlabel('X axis');

ylabel('Y axis');

zlabel('Height [m]');

grid on

% coordinates of 1st blade

X_b1=[X11(1,1) L_blade*cos(omega_rot*t1(i))+X11(1,1)];
Y_b1=[X22(1,1)-5 X22(1,1)-5];
Z_b1=[L1 L1+L_blade*sin(omega_rot*t1(i))];

plot3(X_b1,Y_b1,Z_b1,'r','LineWidth',5)

% coordinates of 2nd blade

X_b2=[X11(1,1) L_blade*cos(omega_rot*t1(i)+2*pi/3)+X11(1,1)];
Y_b2=[X22(1,1)-5 X22(1,1)-5];
Z_b2=[L1 L1+L_blade*sin(omega_rot*t1(i)+2*pi/3)];

plot3(X_b2,Y_b2,Z_b2,'r','LineWidth',5)

% coordinates of 3rd blade

X_b3=[X11(1,1) L_blade*cos(omega_rot*t1(i)+4*pi/3)+X11(1,1)];
Y_b3=[X22(1,1)-5 X22(1,1)-5];
Z_b3=[L1 L1+L_blade*sin(omega_rot*t1(i)+4*pi/3)];

plot3(X_b3,Y_b3,Z_b3,'r','LineWidth',5)

text(0.6*L_blade,1.1*(L1+L_blade),sprintf('t = %.2f s',t1(i)))

axis equal 

hold off

axis(1*[-100 100 -100 100 0 L1+L_blade]);

H(i) = getframe;

end



end

