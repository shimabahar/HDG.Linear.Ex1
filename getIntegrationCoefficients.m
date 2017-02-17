% This function creates integration coefficients for Gauss Legendre integration method
function [x w np] = getIntegrationCoefficients()
global n

% The integration method using np points is accurate for polynomials 
% up to degree 2*np-1. We choose number of points (np) for Gauss Legendre 
% integration method based on polynomiel degree (n-1). Note that in our
% calculations we have phi^2 so when degree of phi is d, we have to 
% consider 2*d as degree in our calculation. 
% | d | 2*d | n (=d+1) | np | accuracy  | 
%  --------------------------------------
% | 0 |  0  |    1     | 1  | 2*1-1 = 1 |
% | 1 |  2  |    2     | 2  | 2*2-1 = 3 |
% | 2 |  4  |    3     | 3  | 2*3-1 = 5 |
% | 3 |  6  |    4     | 4  | 2*4-1 = 7 |
np = n;

% Gauss Legendre integration method coefficients based on number of points (np)
if np == 1
  x = [0];
  w = [2];
elseif np == 2
  x = [-1/(3)^(1/2), 1/(3)^(1/2)];
  w = [1, 1];
elseif np == 3
  x = [-(3/5)^(1/2), 0, (3/5)^(1/2)];
  w = [5/9, 8/9, 5/9];
elseif np == 4
  x = [0.8611363116, 0.3399810436, -0.3399810436, -0.8611363116];
  w = [0.3478548451, 0.6521451549, 0.6521451549, 0.3478548451];
end