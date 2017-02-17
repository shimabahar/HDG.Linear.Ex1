% This function calculate left hand side matrix M
function M = getLhsMatrix(A, B, D1, D2, C1, C2, C3, C4T, C5T, C6T, E1, E2)
global n ne dt

%     [(1/dt)*A-D1      0     -B^T     0      C1]
%     [     0         B-D2      A     -C2      0]
% M = [     B          A        0      0     -C3]
%     [     0         C4T       0      E1      0]
%     [    C5T         0       C6T     0     -E2]

m = 3*n*ne + 2*ne - 1; 
M = zeros(m, m);

M(1:n*ne, 1:n*ne) = (1/dt)*A - D1;
M(1:n*ne, 2*n*ne + 1:3*n*ne) = -B';
M(1:n*ne, 3*n*ne + ne + 1:m) = C1;

M(n*ne + 1:2*n*ne, n*ne + 1:2*n*ne) = B - D2;
M(n*ne + 1:2*n*ne, 2*n*ne + 1:3*n*ne) = A;
M(n*ne + 1:2*n*ne, 3*n*ne + 1:3*n*ne + ne) = -C2;

M(2*n*ne + 1:3*n*ne, 1:n*ne) = B;
M(2*n*ne + 1:3*n*ne, n*ne + 1:2*n*ne) = A;
M(2*n*ne + 1:3*n*ne, 3*n*ne + ne + 1:m) = -C3;

M(3*n*ne + 1:3*n*ne+ne, n*ne + 1:2*n*ne) = C4T;
M(3*n*ne + 1:3*n*ne+ne, 3*n*ne + 1:3*n*ne + ne) = E1;

M(3*n*ne + ne + 1:m, 1:n*ne) = C5T;
M(3*n*ne + ne + 1:m, 2*n*ne + 1:3*n*ne) = C6T;
M(3*n*ne + ne + 1:m, 3*n*ne + ne + 1:m) = -E2;