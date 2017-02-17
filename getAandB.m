% This function calculate matrices A and B
function [A B] = getAandB()
global n ne

A = zeros(n*ne, n*ne);
B = zeros(n*ne, n*ne);
As = zeros(n, n); % sub-matrix A
Bs = zeros(n, n); % sub-matrix B

for i = 1:n
  for j = 1:n
    As(i, j) = integratePhi(i, j);
    Bs(i, j) = integrateDphi(i, j);
  end
end

% Building matrices A and B using sub-matrices As and Bs
k = n;
for i = 1:n:n*ne
  A(i:k, i:k) = As;
  B(i:k, i:k) = Bs;
  k = k+n;
end