% This function calculates the L2 error for U - Uexact
function error = getL2Error(U, step)

global n ne xL dt h

error = 0;
z = 0;
k = 1;

for i = 1:ne  % differenet elements
  a = xL + (k-1)*h;
  b = xL + k*h;

  u = zeros(n, 1);
  for j = 1:n  % differenet base functions
    u(j) = U(j+z);
  end

  error = integrateForL2Error(u, i, a, b, step) + error;  
  k = k + 1;
  z = z + n;
end

error = error^(1/2);