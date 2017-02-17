% This function calculates the right hand side vector R
function R = getRhsVector(A, U0, step)
global n ne xL xR dt tau sigma

% For the first step with the exact initial 
% condition, there is no need to multiply to A
if step == 0
  A = eye(n*ne);
end

%%%% calculate g(v) and K(z) %%%%
g1 = zeros(n, 1);
g2 = zeros(n, 1);
 
for i = 1:n
  g1(i) = guxL(step*dt) * phi(i, 1, xL);
  g2(i) = guxR(step*dt) * phi(i, ne, xR);
end

G = zeros(n*ne, 1);
K = zeros(n*ne, 1);

G(1:n) = tau * g1;
G(n*ne-n+1:n*ne) = tau * g2;

K(1:n) = -g1;
K(n*ne-n+1:n*ne) = g2;

%%%% calculate h(w) %%%%
h1 = zeros(n, 1);

for i = 1:n
  h1(i) = (-1-sigma)*gqxL(step*dt)*phi(i, 1, xL);
end

H = zeros(n*ne, 1);
H(1:n) = h1;

%%%% calculate s(nu) %%%%
S = zeros(ne-1, 1);

%%%% calculate the right hand side vector %%%%  
R = zeros(3*n*ne+2*ne-1, 1);

R(1:n*ne) = -G + (1/dt)*A*U0;
R(n*ne+1:2*n*ne) = H;
R(2*n*ne+1:3*n*ne) = K;
R(3*n*ne+ne+1:3*n*ne+2*ne-1) = S;