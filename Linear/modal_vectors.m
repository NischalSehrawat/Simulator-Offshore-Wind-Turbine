function [v]=modal_vectors(v_xz,v_yz,v_tor,n_modes)
               
    %   MODAL_VECTORS gives the modal vector according to the number of modes selected
    %   It takes the modal vectors of xy,xz and tor plane as input and returns
    %   combined modal vectors for the entire model;

n1=size(v_xz,1)+size(v_yz,1)+size(v_tor,1);

n2=3*n_modes;

v=zeros(n1,n2);

V_xz=v_xz(:,1:n_modes);

V_yz=v_yz(:,1:n_modes);

V_tor=v_tor(:,1:n_modes);

v(1:size(v_xz,1),1:n_modes)=V_xz;

v(size(v_xz,1)+1:size(v_xz,1)+size(v_yz,1),n_modes+1:2*n_modes)=V_yz;

v(size(v_xz,1)+size(v_xz,1)+1:n1,2*n_modes+1:n2)=V_tor;


end

