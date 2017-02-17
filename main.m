%{
Matlab code to solve KdV equation using HDG method.

More details can be found in the following paper:
A Hybridized Discontinuous Galerkin Method for the Nonlinear
Korteweg–de Vries Equation - Ali Samii et al.

Author(s):
Date:
%}

clc
clear all

tic

global n ne xL xR dt T tau sigma h X 

%%%%%% Problem configuration %%%%%%
degree = 2; % max degree of the base functions (phi)
n = degree + 1;
ne = 10; % number of elements in x domain
xL = 0;  % left hand side
xR = pi; % right hand side
dt = 10^-4;
T = 0.1;
tau = -10;
sigma = 10;
h = double((xR-xL)/ne);
X = xL:h:xR;

%%%% Build coefficient matrices and the left hand side matrix %%%%
[A B] = getAandB();
[D1 D2] = getD1andD2();
[C1 C2 C3] = getC1andC2andC3();
[C4T C5T C6T] = getC4TandC5TandC6T();
[E1 E2] = getE1andE2();
M = getLhsMatrix(A, B, D1, D2, C1, C2, C3, C4T, C5T, C6T, E1, E2);

%%%% Initial condition U0  %%%%
U0 = zeros(n*ne, 1);
k = 0;
for i = 1:ne % element
  for j = 1:n % phi_j
    U0(j+k) = integrateU0(i, j); 
  end
  k = k + n;
end

%%%% Solve the system and calculate the error  %%%%
for step = 0:T/dt
  R = getRhsVector(A, U0, step);
  U = M \ R;

  U0 = U(1:n*ne);
  Q = U(n*ne+1:2*n*ne);
  P = U(2*n*ne+1:3*n*ne);
  lambda = U(3*n*ne+1:3*n*ne+ne-1);
  psi = U(3*n*ne+ne:3*n*ne+2*ne-1);

  error = getL2Error(U0, step);
  disp(error);
end

toc