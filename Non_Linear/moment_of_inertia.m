function [I,J] = moment_of_inertia( d_o,d_in)
I=(pi/64)*(d_o^4-d_in^4);
J=I*2;

end

