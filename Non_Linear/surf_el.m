function zeta = surf_el( zeta_A,omega,phase,wavenum,t)
x=0;
for i=1:length(t)
    zeta(i)=sum(zeta_A.*cos(omega.*t(i)+phase-wavenum.*x));
end
end

