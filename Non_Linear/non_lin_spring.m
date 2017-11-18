function K_non_lin = non_lin_spring(Num_el,depth_embed,k,D,c1,c2,c3,gamma,y)
                        % This function generates the non linear
                        % spring tables
                        % It takes the following inputs
                        % Num_el=number of elements for discretizing embedded pile
                        % depth_embed=embedment depth
                        % D=pile outer diameter
                        % c1,c2,c3 as obtained from DNV
                        % Gamma=soil specific weight in N/m^3
                        % k=soil stiffness parameter from DNV curves dependent on the angle of
                        % internal friction
n_nodes=Num_el+1;

l=depth_embed/Num_el;

% initial spring stiffness on 1st node is zero

k_1=0;

% initial spring stiffness on the last node 

k_last=Num_el*l*k*l/2;

% initial spring stiffness on intermediate nodes

k_int=zeros(n_nodes-2,1);

for i=2:1:n_nodes-1

    k_int(i-1)=(i-1)*l*k*l;
end

K_final=[k_1;k_int;k_last];

% Generating p-y curves

% Calculating ultimate strength of soil

x=0:l:depth_embed;

p_us=((c1.*x+c2.*D).*gamma.*x)'; % shallow depth ultimate strength

p_ud=(c3.*D.*gamma.*x)'; % deep depth ultimate strength

% Ultimate strength is minimum of shallow and deep depth strength

parfor i=1:1:n_nodes

p_u(i,1)=min(p_us(i,1),p_ud(i,1));

end

b=1;

for j=0:l:depth_embed

    P(b,:)=0.9.*p_u(b,1).*tanh(k.*j.*y/(0.9.*p_u(b,1)));

    b=b+1;
end

P(1,:)=0;

% Epy

parfor u=1:1:n_nodes

    E_py(u,:)=P(u,:)./y;
end
% non linear springs
K_non_lin=zeros(size(E_py,1),size(E_py,2));
K_non_lin(1,:)=0;
K_non_lin(2:size(E_py,1)-1,:)=E_py(2:size(E_py,1)-1,:).*l;
K_non_lin(end,:)=E_py(end,:).*l/2;

% Defining initial stiffness for the non linear springs

parfor i=1:1:n_nodes

    K_non_lin(i,1)=K_final(i,1);
end

end

