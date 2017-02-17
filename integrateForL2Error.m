% This function integrates (U-Uexact)^2 on [a,b] using 
% Gauss Legendre integration method where U = sum(u_j*phi_j)
function out = integrateForL2Error(u, e, a, b, step)
global n dt

[x w np] = getIntegrationCoefficients();

% Gauss Legendre integration method on element 1 (on [0,h])
gl = 0;

for i = 1:np
  U = 0;
  for j = 1:n
    U = u(j) * phi(j, e, (a+b)/2+(b-a)/2*x(i)) + U;
  end
  gl = w(i) * (U - Uexact((a+b)/2+(b-a)/2*x(i), step*dt))^2 + gl;
end

out = ((b-a)/2)^2 * gl;


    