function y = BeamAssemble(M,m,i,j)
%BeamAssemble This function assembles the element mass
% matrix m of the beam element with nodes
% i and j into the global Mass matrix M.
% This function returns the global Mass
% matrix M after the element mass matrix
% m is assembled.
M(2*i-1,2*i-1) = M(2*i-1,2*i-1) + m(1,1);
M(2*i-1,2*i) = M(2*i-1,2*i) + m(1,2);
M(2*i-1,2*j-1) = M(2*i-1,2*j-1) + m(1,3);
M(2*i-1,2*j) = M(2*i-1,2*j) + m(1,4);
M(2*i,2*i-1) = M(2*i,2*i-1) + m(2,1);
M(2*i,2*i) = M(2*i,2*i) + m(2,2);
M(2*i,2*j-1) = M(2*i,2*j-1) + m(2,3);
M(2*i,2*j) = M(2*i,2*j) + m(2,4);
M(2*j-1,2*i-1) = M(2*j-1,2*i-1) + m(3,1);
M(2*j-1,2*i) = M(2*j-1,2*i) + m(3,2);
M(2*j-1,2*j-1) = M(2*j-1,2*j-1) + m(3,3);
M(2*j-1,2*j) = M(2*j-1,2*j) + m(3,4);
M(2*j,2*i-1) = M(2*j,2*i-1) + m(4,1);
M(2*j,2*i) = M(2*j,2*i) + m(4,2);
M(2*j,2*j-1) = M(2*j,2*j-1) + m(4,3);
M(2*j,2*j) = M(2*j,2*j) + m(4,4);
y = M;
end



