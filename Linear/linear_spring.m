function K = linear_spring(K1,Num_el,depth_embed,k)
                        % This function modifies the embedded pile stiffness matrices 
                        % and adds the soil spring stiffness to the respective place in the system stiffness matrix.
                        % It takes the following inputs
                        % K1=Embedded pile stiffness matrix
                        % Num_el=number of elements for discretizing embedded pile
                        % k=soil stiffness parameter from DNV curves dependent on the angle of
                        % internal friction
n_nodes =Num_el+1;
l=depth_embed/Num_el;
% spring constant on 1st node is zero
k_1=0;
% spring constant on the last node 
k_last=Num_el*l*k*l/2;
% spring constant on intermediate nodes
k_int=zeros(n_nodes-2,1);
for i=2:1:n_nodes-1
k_int(i-1)=(i-1)*l*k*l;
end
soil_spring_coefficients=[k_1;k_int;k_last];
j=1;
for i=1:1:n_nodes
    K1(j,j)=K1(j,j)+soil_spring_coefficients(i,1);
    j=j+2;
end
K=K1;
end

