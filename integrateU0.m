% This function integrates u0*phi_j on the i'th element (i.e., [a,b])
function out = integrateU0(e, j)
global xL h

a = xL + (e-1)*h;
b = xL + e*h;

c0= xL - h/2;
c1= xL + h/2; % center of first element

if j == 1
  out = -cos(b) + cos(a);
elseif j == 2
  out = -(2/h)*cos(b)*c1 ...
       + (2/h)*cos(a)*c0 ...
       + (2/h)*(sin(b) - sin(a));
elseif j == 3
  out = (4/(h^2))*(-cos(b)*c1^2 + cos(a)*c0^2) ...
      + (8/(h^2))*(sin(b)*c1 - sin(a)*c0) ...
      - (8/(h^2))*(-cos(b) + cos(a));
end