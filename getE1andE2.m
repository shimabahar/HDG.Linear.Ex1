% This function calculate matrices E1 and E2
function [E1 E2] = getE1andE2()
global ne tau sigma

E1 = zeros(ne, ne);
E1(1:ne-1, 1:ne-1) = eye(ne-1) * (-2*sigma);
E1(ne, ne) = 1 - sigma;

E2 = zeros(ne-1, ne-1);
E2 = (2*tau) * eye(ne-1);
