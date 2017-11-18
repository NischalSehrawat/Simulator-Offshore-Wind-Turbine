function [zeta_A,omega,phase,wavenum,t_sim]=create_waves(spectrum_type,hs,T0,depth,nfreq,freq_cutoff,rand_seed)
g=9.81;
omega_peak=2*pi/T0;
omega_max=omega_peak*freq_cutoff;
% freq step
delta_omega=omega_max/nfreq;
% simulation time
max_time_sim=2*pi/delta_omega;
% create frequency vector
omega_vec=delta_omega:delta_omega:omega_max;
% set random generator
rng(rand_seed,'v5uniform');
% create phase
phase=2*pi*rand(1,nfreq);
% spectrum
s_omega=A6_wave_spectrum(spectrum_type,[hs T0],omega_vec',0);
% zeta
parfor ii=1:nfreq
zeta_A(ii)=sqrt(2*s_omega(ii)*delta_omega);
% wave number
wavenum_0=omega_vec(ii)^2/g;
wavenum(ii)=fzero(@(k)k*tanh(k*depth)-omega_vec(ii)^2/g,[wavenum_0 1e10]);
end
omega=omega_vec;
t_sim=max_time_sim;
end

