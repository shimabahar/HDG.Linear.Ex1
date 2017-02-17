% This function calculate matrices C4T, C5T, and C6T
function [C4T C5T C6T] = getC4TandC5TandC6T()
global n ne xL xR tau sigma h

C4T = zeros(ne, n*ne);
C5T = zeros(ne-1, n*ne);
C6T = zeros(ne-1, n*ne);

c4t = zeros(1, 2*n);
c5t = zeros(1, 2*n);
c6t = zeros(1, 2*n);

for i = 1:n
  c4t(i) = phi(i, 1, xL+h);
  c4t(i+n) = phi(i, 2, xL+h);

  c5t(i) = c4t(1, i);
  c5t(i+n) = c4t(1, i+n);

  c6t(i) = c4t(1, i);
  c6t(i+n) = -c4t(1, i+n);
end

% Building matrix C4T (except the last row) using vector c4t
for i = 1:ne-1
  C4T(i, n*(i-1)+1:n*(i+1)) = sigma * c4t;
end

% Building the last row of matrix C4T
for i = 1:n
  C4T(ne, n*(ne-1)+i:n*ne) = sigma * c4t(i) - phi(i, ne, xR);
end

% Building matrix C5T using vector c5t
for i = 1:ne-1
  C5T(i, n*(i-1)+1:n*(i+1)) = tau * c5t;
end

% Building matrix C6T using vector c6t
for i = 1:ne-1
  C6T(i, n*(i-1)+1:n*(i+1)) = c6t;
end
