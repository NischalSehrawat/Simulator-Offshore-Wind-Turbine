function [Nat_freq,Time_period,v]= nat_freq(M_system,K_system)
% nat_freq gives a vector consisiting of the natural frequencies of the OWT
% in Hz and time period in s
% It takes M_system and K_system as the input
% and gives a vector of natural frequencies as output
A=M_system\K_system;
%Obtain eigenvalues and eigenvectors of A
[V,D]=eig(A);
[omega_vec,i]=sort(sqrt(diag(D)));
Nat_freq=omega_vec./(2*pi);
% Nat_freq(2:2:end)=[];
Time_period=1./Nat_freq;
v=V(:,i);
end

