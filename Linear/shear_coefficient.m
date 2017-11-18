function   alfa =shear_coefficient(D_out,D_in,nu)
% this function calcultes the shear coefficient of any section.
% It takes D_out, D_in and poisson's ratio as input,
b=0.5*D_out;a=0.5*D_in;
c=6*(a^2+b^2)^2*(1+nu)^2;
d=7*a^4+34*a^2*b^2+7*b^4+nu*(12*a^4+48*a^2*b^2+12*b^4)+nu^2*(4*a^4+16*a^2*b^2+4*b^4);
alfa=d/c;
end

