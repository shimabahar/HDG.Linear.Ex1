% This function integrates Dphi on [0,h] using Gauss Legendre integration method
function out = integrateDphi(i, j)
global n h

[x w np] = getIntegrationCoefficients();

% Gauss Legendre integration method on element 1 (on [0,h])
gl = 0;
for k = 1:np
  gl = w(k) * Dphi(i, 1, (h/2)*(x(k)+1)) * phi(j, 1, (h/2)*(x(k)+1)) + gl;
end
out = h/2 * gl;
