% This function calculate matrices D1 and D2
function [D1 D2] = getD1andD2()
global n ne xL tau sigma h

D = zeros(n*ne, n*ne);
Ds = zeros(n, n); % sub-matrix D

for i = 1:n
  for j = 1:n
    Ds(i, j) = phi(j, 1, xL+h)*phi(i, 1, xL+h) + phi(j, 1, xL)*phi(i, 1, xL);
  end
end

% Building matrix D using sub-matrix Ds
k = n;
for i = 1:n:n*ne
  D(i:k, i:k) = Ds;
  k = k+n;
end

D1 = tau * D;
D2 = sigma * D;